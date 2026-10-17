"""Recheck approximate sulfur events with the thorough connected-root search.

This audit is intentionally separate from D_01_differentiation_s.ipynb. It reads
that notebook's exported baseline history, reconstructs the exact incoming
inventory for each approximate event, and writes a non-destructive audit table.
"""
from __future__ import annotations

import json
import multiprocessing as mp
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

WORKFLOW_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = next(p for p in WORKFLOW_DIR.parents
                    if (p / "src" / "accrediff").is_dir())
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import accrediff as ad

TABLES = WORKFLOW_DIR / "Tables"
PHASE_ITEMS = [
    "MgO", "Al2O3", "CaO", "FeO", "NiO", "SiO2", "S_sil",
    "Fe", "Ni", "Si", "O", "S_met",
]
MANTLE_ITEMS = ["MgO", "Al2O3", "CaO", "FeO", "NiO", "SiO2", "S_sil"]
CORE_ITEMS = ["Fe", "Ni", "Si", "O", "S_met"]


def _scale(values, factor):
    return {name: float(values[name]) * factor for name in PHASE_ITEMS}


def _add(*values):
    return {name: sum(float(item[name]) for item in values) for name in PHASE_ITEMS}


def _subset(values, names):
    names = set(names)
    return {name: float(values[name]) if name in names else 0.0 for name in PHASE_ITEMS}


def split_event(target, impactor, row):
    mt, mi = float(target["m_e"]), float(impactor["m_e"])
    target_total = _scale(target, mt)
    impactor_total = _scale(impactor, mi)
    target_mantle = _subset(target_total, MANTLE_ITEMS)
    impactor_mantle = _subset(impactor_total, MANTLE_ITEMS)
    impactor_core = _subset(impactor_total, CORE_ITEMS)
    event = row["events"]
    if event == "small":
        mantle_fraction, core_fraction = 0.0, 0.7
    elif event == "global":
        mantle_fraction, core_fraction = float(row["k_mantle"]) * 0.05, 0.7
    elif event == "partial":
        mantle_fraction, core_fraction = float(row["k_mantle"]), 0.3
    else:
        raise ValueError(event)
    return _add(
        impactor_mantle,
        _scale(impactor_core, core_fraction),
        _scale(target_mantle, mantle_fraction),
    )


def reconstruct_cases():
    early = pd.read_csv(TABLES / "03_embryo_bulk_composition.csv", index_col=0)
    impacts = pd.read_excel(TABLES / "05_impact_events_S.xlsx", sheet_name=None)
    histories = pd.read_excel(TABLES / "06_differentiation_S.xlsx", sheet_name=None)
    initial = early[["i", "m_e"] + PHASE_ITEMS + ["S"]].copy().set_index("i")
    initial.index = initial.index.astype(int)
    cases = []
    for planet, events in impacts.items():
        state = initial.copy()
        history = histories[planet]
        if len(events) != len(history):
            raise ValueError(f"{planet}: impact/history row mismatch")
        for ordinal, ((_, row), (_, outcome)) in enumerate(
            zip(events.iterrows(), history.iterrows()), start=1
        ):
            target_id, impactor_id = int(row["target_id"]), int(row["impactor_id"])
            re_eq = split_event(state.loc[target_id], state.loc[impactor_id], row)
            if bool(outcome["accepted_approximate"]):
                cases.append({
                    "case_id": len(cases),
                    "planet": planet,
                    "event_number": ordinal,
                    "Time": float(row["Time"]),
                    "event": row["events"],
                    "target_id": target_id,
                    "impactor_id": impactor_id,
                    "P_equil": float(row["P_equil"]),
                    "T_equil": float(row["T_equil"]),
                    "baseline_KD_O_target": float(outcome["KD_O_target"]),
                    "baseline_KD_O_model": float(outcome["KD_O_model"]),
                    "baseline_O_log_residual": float(outcome["O_KD_log_residual"]),
                    "re_eq": re_eq,
                })
            for name in ["m_e"] + PHASE_ITEMS + ["S"]:
                state.loc[target_id, name] = outcome[name]
    return cases


def audit_case(case):
    started = time.perf_counter()
    calc = ad.SulfurKDCalculator()
    scale = sum(float(case["re_eq"][name]) for name in PHASE_ITEMS)
    normalized = {name: float(case["re_eq"][name]) / scale for name in PHASE_ITEMS}
    params = ad.sulfur_params_from_composition(
        normalized,
        P=case["P_equil"],
        T=case["T_equil"],
        major_kd_calculator=calc,
    )
    target = calc.get_KD("O", case["P_equil"], case["T_equil"])
    solver = ad.ForwardKDOSolverSulfur(
        params,
        P=case["P_equil"],
        T=case["T_equil"],
        sulfur_model="boujibar",
        sulfur_kd_calculator=calc,
        root_search="thorough",
        allow_nearest=True,
        sulfur_tol=1e-9,
        coupled_tol=1e-9,
    )
    try:
        result = solver.solve_x_for_KD_O(target, tol=1e-8, grid_N=200)
        status = "fast_missed_exact_root" if result.oxygen_converged else "no_exact_root_found_thorough"
        output = {
            key: value for key, value in case.items() if key != "re_eq"
        }
        output.update({
            "audit_status": status,
            "thorough_oxygen_converged": bool(result.oxygen_converged),
            "thorough_overall_converged": bool(result.overall_converged),
            "thorough_KD_O_model": float(result.KD_O_),
            "thorough_O_log_residual": float(result.residuals.get("O_KD_log", result.residual)),
            "thorough_sulfur_residual": float(result.sulfur_residual),
            "thorough_D_S_mole": float(result.D_S_mole),
            "thorough_FeO": float(result.x_ * scale),
            "thorough_NiO": float(result.y_ * scale),
            "thorough_SiO2": float(result.z_ * scale),
            "thorough_Fe": float(result.a_ * scale),
            "thorough_Ni": float(result.b_ * scale),
            "thorough_Si": float(result.c_ * scale),
            "thorough_O": float(result.d_ * scale),
            "thorough_S_sil": float(result.S_sil_ * scale),
            "thorough_S_met": float(result.S_met_ * scale),
            "audit_seconds": time.perf_counter() - started,
            "audit_error": "",
        })
        return output
    except Exception as exc:
        output = {key: value for key, value in case.items() if key != "re_eq"}
        output.update({
            "audit_status": "audit_failed",
            "thorough_oxygen_converged": False,
            "thorough_overall_converged": False,
            "audit_seconds": time.perf_counter() - started,
            "audit_error": f"{type(exc).__name__}: {exc}",
        })
        return output


def main():
    cases = reconstruct_cases()
    print(f"Reconstructed {len(cases)} approximate events", flush=True)
    results = []
    workers = min(4, max(1, mp.cpu_count() - 1))
    with mp.get_context("spawn").Pool(workers) as pool:
        for completed, result in enumerate(pool.imap_unordered(audit_case, cases), start=1):
            results.append(result)
            print(
                f"[{completed}/{len(cases)}] {result['planet']} event "
                f"{result['event_number']}: {result['audit_status']}, "
                f"residual={result.get('thorough_O_log_residual', np.nan):.6g}",
                flush=True,
            )
    audit = pd.DataFrame(results).sort_values("case_id").reset_index(drop=True)
    csv_path = TABLES / "06_thorough_root_audit.csv"
    xlsx_path = TABLES / "06_thorough_root_audit.xlsx"
    audit.to_csv(csv_path, index=False)
    with pd.ExcelWriter(xlsx_path) as writer:
        audit.to_excel(writer, sheet_name="event_audit", index=False)
        audit.groupby(["planet", "audit_status"]).size().unstack(fill_value=0).to_excel(
            writer, sheet_name="summary"
        )
    print(audit["audit_status"].value_counts(dropna=False), flush=True)
    print(f"Saved {csv_path}", flush=True)
    print(f"Saved {xlsx_path}", flush=True)


if __name__ == "__main__":
    main()
