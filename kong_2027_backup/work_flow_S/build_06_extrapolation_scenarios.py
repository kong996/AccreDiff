"""Build independent event-level extrapolation end members for audited no-root events."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from audit_06_thorough import PHASE_ITEMS, TABLES, reconstruct_cases


def elemental(values):
    return {
        "Fe": values["FeO"] + values["Fe"],
        "Ni": values["NiO"] + values["Ni"],
        "Si": values["SiO2"] + values["Si"],
        "O": values["FeO"] + values["NiO"] + 2.0 * values["SiO2"] + values["O"],
        "S": values["S_sil"] + values["S_met"],
        "MgO": values["MgO"],
        "Al2O3": values["Al2O3"],
        "CaO": values["CaO"],
    }


def conservation_error(before, after):
    b, a = elemental(before), elemental(after)
    scale = max(max(abs(v) for v in b.values()), 1e-30)
    return max(abs(a[k] - b[k]) for k in b) / scale


def main():
    audit = pd.read_csv(TABLES / "06_thorough_root_audit.csv")
    cases = {case["case_id"]: case for case in reconstruct_cases()}
    if set(audit["audit_status"]) != {"no_exact_root_found_thorough"}:
        raise ValueError("Scenario construction requires completed no-root audit results.")

    nearest_rows, preserved_rows, comparison_rows = [], [], []
    for _, audited in audit.iterrows():
        case = cases[int(audited["case_id"])]
        incoming = case["re_eq"]
        nearest = {
            "MgO": incoming["MgO"],
            "Al2O3": incoming["Al2O3"],
            "CaO": incoming["CaO"],
            "FeO": audited["thorough_FeO"],
            "NiO": audited["thorough_NiO"],
            "SiO2": audited["thorough_SiO2"],
            "Fe": audited["thorough_Fe"],
            "Ni": audited["thorough_Ni"],
            "Si": audited["thorough_Si"],
            "O": audited["thorough_O"],
            "S_sil": audited["thorough_S_sil"],
            "S_met": audited["thorough_S_met"],
        }
        common = {
            "case_id": int(audited["case_id"]),
            "planet": audited["planet"],
            "event_number": int(audited["event_number"]),
            "Time": audited["Time"],
            "event": audited["event"],
            "target_id": int(audited["target_id"]),
            "impactor_id": int(audited["impactor_id"]),
            "P_equil": audited["P_equil"],
            "T_equil": audited["T_equil"],
            "KD_O_target": audited["baseline_KD_O_target"],
        }
        nearest_row = {
            **common,
            "scenario": "thorough_nearest_physical",
            "oxygen_constraint_satisfied": False,
            "Ni_Si_S_equilibrium_satisfied": True,
            "O_KD_log_residual": audited["thorough_O_log_residual"],
            "max_conservation_error": conservation_error(incoming, nearest),
            **nearest,
        }
        preserved_row = {
            **common,
            "scenario": "phase_preserving_no_repartition",
            "oxygen_constraint_satisfied": False,
            "Ni_Si_S_equilibrium_satisfied": False,
            "O_KD_log_residual": np.nan,
            "max_conservation_error": conservation_error(incoming, incoming),
            **incoming,
        }
        nearest_rows.append(nearest_row)
        preserved_rows.append(preserved_row)
        comparison_rows.append({
            **common,
            "thorough_O_log_residual": audited["thorough_O_log_residual"],
            "nearest_S_sil_fraction": nearest["S_sil"] / (nearest["S_sil"] + nearest["S_met"]),
            "preserved_S_sil_fraction": incoming["S_sil"] / (incoming["S_sil"] + incoming["S_met"])
            if incoming["S_sil"] + incoming["S_met"] > 0 else np.nan,
            "delta_S_sil_fraction": (
                nearest["S_sil"] / (nearest["S_sil"] + nearest["S_met"])
                - incoming["S_sil"] / (incoming["S_sil"] + incoming["S_met"])
            ) if incoming["S_sil"] + incoming["S_met"] > 0 else np.nan,
        })

    nearest_df = pd.DataFrame(nearest_rows)
    preserved_df = pd.DataFrame(preserved_rows)
    comparison_df = pd.DataFrame(comparison_rows)
    if nearest_df["max_conservation_error"].max() > 1e-8:
        raise RuntimeError("Thorough-nearest scenario violates elemental conservation.")
    if preserved_df["max_conservation_error"].max() > 1e-12:
        raise RuntimeError("Phase-preserving scenario violates elemental conservation.")

    descriptions = pd.DataFrame([
        {
            "scenario": "thorough_nearest_physical",
            "purpose": "Model-constrained extrapolation end member",
            "definition": "Closest state found by thorough connected-root search; Ni/Si/S equilibrium and physical audit pass, oxygen KD does not.",
            "use": "Sensitivity branch only; never label as full equilibrium.",
        },
        {
            "scenario": "phase_preserving_no_repartition",
            "purpose": "No-chemical-redistribution end member",
            "definition": "Carry the re-equilibrating inventory forward with its incoming silicate/metal allocation unchanged.",
            "use": "Sensitivity branch only; bounds the effect of forcing repartition where the model has no full root.",
        },
    ])
    output = TABLES / "06_extrapolation_scenarios.xlsx"
    with pd.ExcelWriter(output) as writer:
        descriptions.to_excel(writer, sheet_name="scenario_definitions", index=False)
        comparison_df.to_excel(writer, sheet_name="event_comparison", index=False)
        nearest_df.to_excel(writer, sheet_name="thorough_nearest", index=False)
        preserved_df.to_excel(writer, sheet_name="phase_preserving", index=False)
    print(f"Saved {output}")
    print(f"Events: {len(comparison_df)}")
    print(f"Nearest max conservation error: {nearest_df['max_conservation_error'].max():.3g}")
    print(f"Maximum |oxygen KD log residual|: {nearest_df['O_KD_log_residual'].abs().max():.6g} dex")
    print(comparison_df["delta_S_sil_fraction"].describe().to_string())


if __name__ == "__main__":
    main()
