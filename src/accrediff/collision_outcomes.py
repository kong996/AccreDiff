"""Collision outcome classification utilities.

This module implements an EDACM-style analytic collision-outcome
classifier based on Leinhardt and Stewart (2012). It is intended as an
optional preprocessing layer for enriching AccreDiff impact-event tables
with non-perfect-merger diagnostics.

Developer:
    @author: Kang Shuai

Organizer:
    @author: Zhihui Kong
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Dict, Iterable, Optional, Tuple

import numpy as np  # type: ignore
import pandas as pd # type: ignore


class CollisionOutcome(IntEnum):
    """Integer collision-outcome labels used by the EDACM classifier."""

    PERFECT_MERGER = 1
    HIT_AND_RUN = 2
    EROSION = 3
    SUPER_CATASTROPHIC = 4
    PARTIAL_ACCRETION = 5
    HIT_AND_RUN_SUPER_CATASTROPHIC = 6


OUTCOME_LABELS = {
    CollisionOutcome.PERFECT_MERGER: "perfect_merger",
    CollisionOutcome.HIT_AND_RUN: "hit_and_run",
    CollisionOutcome.EROSION: "erosion",
    CollisionOutcome.SUPER_CATASTROPHIC: "super_catastrophic",
    CollisionOutcome.PARTIAL_ACCRETION: "partial_accretion",
    CollisionOutcome.HIT_AND_RUN_SUPER_CATASTROPHIC: (
        "hit_and_run_super_catastrophic"
    ),
}


@dataclass(frozen=True)
class EDACMConfig:
    """Configuration parameters for the EDACM-style classifier.

    Parameters are expressed in SI units. The two disruption-law parameter
    sets follow the small- and large-body branches used in the original
    prototype script.
    """

    rho: float = 3000.0
    rho_ref: float = 1000.0
    G: float = 6.67408e-11
    c_star_small: float = 5.0
    mu_small: float = 0.37
    c_star_large: float = 1.9
    mu_large: float = 0.36
    transition_radius: float = 1.0e6
    supercatastrophic_fraction: float = 0.1
    rocky_supercatastrophic_prefactor: float = 0.457
    rocky_supercatastrophic_exponent: float = -1.24


def impact_parameter_from_angle(theta_rad: np.ndarray | float) -> np.ndarray:
    """Convert impact angle in radians to impact parameter.

    This helper is provided for workflows that store collision geometry as an
    impact angle rather than the dimensionless impact parameter ``b``.
    """

    return np.sin(theta_rad)


def _broadcast_inputs(
    *values: Iterable[float] | float,
) -> Tuple[np.ndarray, ...]:
    arrays = [np.asarray(value, dtype=float) for value in values]
    return tuple(np.ravel(array) for array in np.broadcast_arrays(*arrays))


def _validate_inputs(
    m_target: np.ndarray,
    m_projectile: np.ndarray,
    r_target: np.ndarray,
    r_projectile: np.ndarray,
    v_impact: np.ndarray,
    impact_parameter: np.ndarray,
    require_target_larger: bool,
) -> None:
    if np.any(m_target <= 0.0) or np.any(m_projectile <= 0.0):
        raise ValueError("Collision masses must be positive.")
    if np.any(r_target <= 0.0) or np.any(r_projectile <= 0.0):
        raise ValueError("Collision radii must be positive.")
    if np.any(v_impact < 0.0):
        raise ValueError("Impact velocities must be non-negative.")
    if np.any((impact_parameter < 0.0) | (impact_parameter >= 1.0)):
        raise ValueError("Impact parameters must satisfy 0 <= b < 1.")
    if require_target_larger:
        if np.any(m_target < m_projectile):
            raise ValueError("Expected m_target >= m_projectile.")
        if np.any(r_target < r_projectile):
            raise ValueError("Expected r_target >= r_projectile.")


def _select_disruption_parameters(
    r_c1: np.ndarray,
    config: EDACMConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    c_star = np.where(
        r_c1 < config.transition_radius,
        config.c_star_small,
        config.c_star_large,
    )
    mu_param = np.where(
        r_c1 < config.transition_radius,
        config.mu_small,
        config.mu_large,
    )
    return c_star, mu_param


def _labels_from_outcomes(outcome_id: np.ndarray) -> np.ndarray:
    return np.asarray(
        [OUTCOME_LABELS[CollisionOutcome(value)] for value in outcome_id],
        dtype=object,
    )


def classify_edacm_outcomes(
    m_target: Iterable[float] | float,
    m_projectile: Iterable[float] | float,
    r_target: Iterable[float] | float,
    r_projectile: Iterable[float] | float,
    v_impact: Iterable[float] | float,
    impact_parameter: Iterable[float] | float,
    config: Optional[EDACMConfig] = None,
    require_target_larger: bool = True,
) -> Dict[str, np.ndarray]:
    """Classify collision outcomes using an EDACM-style analytic model.

    Parameters
    ----------
    m_target, m_projectile
        Target and projectile masses in kg.
    r_target, r_projectile
        Target and projectile radii in m.
    v_impact
        Impact velocity in m s-1.
    impact_parameter
        Dimensionless impact parameter ``b`` with ``0 <= b < 1``.
    config
        Optional :class:`EDACMConfig` instance.
    require_target_larger
        If ``True``, validate that the target is at least as massive and large
        as the projectile, matching the GENGA-oriented prototype convention.

    Returns
    -------
    dict
        A dictionary of NumPy arrays containing outcome IDs, readable labels,
        grazing and catastrophic flags, disruption thresholds, impact energies,
        and largest-remnant masses.
    """

    config = config or EDACMConfig()
    (
        m_target,
        m_projectile,
        r_target,
        r_projectile,
        v_impact,
        impact_parameter,
    ) = _broadcast_inputs(
        m_target,
        m_projectile,
        r_target,
        r_projectile,
        v_impact,
        impact_parameter,
    )
    _validate_inputs(
        m_target,
        m_projectile,
        r_target,
        r_projectile,
        v_impact,
        impact_parameter,
        require_target_larger=require_target_larger,
    )

    n_events = m_target.size
    outcome_id = np.zeros(n_events, dtype=int)

    # Interacting projectile mass fraction for oblique impacts.
    b_length = (r_target + r_projectile) * impact_parameter
    l_projectile = r_target + r_projectile - b_length
    alpha = (
        (3.0 * r_projectile * l_projectile**2 - l_projectile**3)
        / (4.0 * r_projectile**3)
    )
    full_projectile_overlap = b_length + r_projectile < r_target
    alpha[full_projectile_overlap] = 1.0
    l_projectile[full_projectile_overlap] = 2.0 * r_projectile[
        full_projectile_overlap
    ]
    alpha = np.clip(alpha, 0.0, 1.0)
    m_projectile_interact = m_projectile * alpha

    # Perfect merging threshold.
    m_oblique = m_target + m_projectile_interact
    r_oblique = (3.0 * m_oblique / (4.0 * np.pi * config.rho)) ** (1.0 / 3.0)
    v_escape_oblique = np.sqrt(2.0 * config.G * m_oblique / r_oblique)
    outcome_id[v_impact < v_escape_oblique] = CollisionOutcome.PERFECT_MERGER

    # Grazing geometry.
    b_crit = r_projectile / (r_target + r_projectile)
    grazing = impact_parameter >= b_crit

    # Catastrophic disruption threshold.
    m_total = m_target + m_projectile
    r_c1 = (3.0 * m_total / (4.0 * np.pi * config.rho_ref)) ** (1.0 / 3.0)
    c_star, mu_param = _select_disruption_parameters(r_c1, config)

    q_rd_gamma_1 = c_star * 4.0 / 5.0 * np.pi * config.rho_ref * config.G * r_c1**2
    v_threshold_gamma_1 = (
        (32.0 / 5.0 * np.pi * c_star * config.rho_ref * config.G) ** 0.5
        * r_c1
    )

    reduced_mass = m_target * m_projectile / m_total
    reduced_mass_alpha = alpha * m_target * m_projectile / (
        alpha * m_projectile + m_target
    )

    gamma = m_projectile / m_target
    gamma_base = (gamma + 1.0) ** 2 / (4.0 * gamma)
    q_rd = q_rd_gamma_1 * gamma_base ** (2.0 / (3.0 * mu_param) - 1.0)
    _ = v_threshold_gamma_1 * gamma_base ** (1.0 / (3.0 * mu_param))

    q_rd_star = q_rd * (reduced_mass / reduced_mass_alpha) ** (
        2.0 - 3.0 * mu_param / 2.0
    )
    v_threshold = np.sqrt(2.0 * q_rd_star * m_total / reduced_mass)

    q_erosion = ((m_target / m_total - 0.5) / (-0.5) + 1.0) * q_rd_star
    v_erosion = np.sqrt(q_erosion * m_total / (0.5 * reduced_mass))

    q_supercat = (
        (config.supercatastrophic_fraction - 0.5) / (-0.5) + 1.0
    ) * q_rd_star
    v_supercat = np.sqrt(q_supercat * m_total / (0.5 * reduced_mass))

    hit_and_run = grazing & (v_impact > v_escape_oblique) & (v_impact < v_erosion)
    outcome_id[hit_and_run] = CollisionOutcome.HIT_AND_RUN

    outcome_id[v_impact > v_erosion] = CollisionOutcome.EROSION
    outcome_id[v_impact > v_supercat] = CollisionOutcome.SUPER_CATASTROPHIC

    partial_accretion = (
        (~grazing) & (v_impact > v_escape_oblique) & (v_impact < v_erosion)
    )
    outcome_id[partial_accretion] = CollisionOutcome.PARTIAL_ACCRETION

    q_r = 0.5 * v_impact**2 * m_target * m_projectile / m_total**2
    m_largest_remnant_fraction = -0.5 * (q_r / q_rd_star - 1.0) + 0.5
    m_largest_remnant = m_total * m_largest_remnant_fraction

    q_rd_out = q_rd_star.copy()
    v_threshold_out = v_threshold.copy()
    q_r_out = q_r.copy()
    m_largest_remnant_fraction_out = m_largest_remnant_fraction.copy()
    m_largest_remnant_out = m_largest_remnant.copy()

    # Reverse-impact treatment for hit-and-run events.
    hit_and_run_idx = np.flatnonzero(outcome_id == CollisionOutcome.HIT_AND_RUN)
    if hit_and_run_idx.size:
        lj_r = l_projectile[hit_and_run_idx]
        r_target_r = r_target[hit_and_run_idx]
        r_projectile_r = r_projectile[hit_and_run_idx]
        c_star_r = c_star[hit_and_run_idx]
        mu_param_r = mu_param[hit_and_run_idx]
        v_impact_r = v_impact[hit_and_run_idx]

        arccos_arg = np.clip((lj_r - r_projectile_r) / r_projectile_r, -1.0, 1.0)
        phi = 2.0 * np.arccos(arccos_arg)
        area_interact = r_projectile_r**2 * (
            np.pi - (phi - np.sin(phi)) / 2.0
        )
        l_target = 2.0 * np.sqrt(
            np.maximum(r_target_r**2 - (r_target_r - lj_r / 2.0) ** 2, 0.0)
        )
        m_target_interact_r = area_interact * l_target * config.rho

        m_target_r = m_projectile[hit_and_run_idx]
        m_projectile_r = m_target_interact_r
        m_total_r = m_target_r + m_projectile_r
        r_c1_r = (
            3.0 * m_total_r / (4.0 * np.pi * config.rho_ref)
        ) ** (1.0 / 3.0)
        q_rd_gamma_1_r = (
            c_star_r * 4.0 / 5.0 * np.pi * config.rho_ref * config.G * r_c1_r**2
        )
        v_threshold_gamma_1_r = (
            (32.0 / 5.0 * np.pi * c_star_r * config.rho_ref * config.G) ** 0.5
            * r_c1_r
        )

        reduced_mass_r = m_target_r * m_projectile_r / m_total_r
        gamma_r = m_projectile_r / m_target_r
        gamma_base_r = (gamma_r + 1.0) ** 2 / (4.0 * gamma_r)
        q_rd_r = q_rd_gamma_1_r * gamma_base_r ** (
            2.0 / (3.0 * mu_param_r) - 1.0
        )
        v_threshold_r = v_threshold_gamma_1_r * gamma_base_r ** (
            1.0 / (3.0 * mu_param_r)
        )

        q_r_reverse = 0.5 * v_impact_r**2 * m_target_r * m_projectile_r / m_total_r**2
        m_largest_remnant_fraction_r = -0.5 * (q_r_reverse / q_rd_r - 1.0) + 0.5

        is_hit_and_run_supercat = (
            m_largest_remnant_fraction_r < config.supercatastrophic_fraction
        )
        outcome_id[
            hit_and_run_idx[is_hit_and_run_supercat]
        ] = CollisionOutcome.HIT_AND_RUN_SUPER_CATASTROPHIC

        q_rd_out[hit_and_run_idx] = q_rd_r
        v_threshold_out[hit_and_run_idx] = v_threshold_r
        q_r_out[hit_and_run_idx] = q_r_reverse

    supercat_idx = outcome_id == CollisionOutcome.SUPER_CATASTROPHIC
    if np.any(supercat_idx):
        m_largest_remnant_fraction_out[supercat_idx] = (
            config.rocky_supercatastrophic_prefactor
            * (q_r_out[supercat_idx] / q_rd_out[supercat_idx])
            ** config.rocky_supercatastrophic_exponent
        )
        m_largest_remnant_out[supercat_idx] = (
            m_total[supercat_idx] * m_largest_remnant_fraction_out[supercat_idx]
        )

    catastrophic = q_r_out > q_rd_out

    return {
        "outcome_id": outcome_id,
        "outcome": _labels_from_outcomes(outcome_id),
        "grazing": grazing,
        "catastrophic": catastrophic,
        "alpha_interact": alpha,
        "m_projectile_interact": m_projectile_interact,
        "q_r": q_r_out,
        "q_rd_star": q_rd_out,
        "v_threshold": v_threshold_out,
        "m_largest_remnant_fraction": m_largest_remnant_fraction_out,
        "m_largest_remnant": m_largest_remnant_out,
        "v_escape_oblique": v_escape_oblique,
        "v_erosion": v_erosion,
        "v_supercatastrophic": v_supercat,
    }


def add_collision_outcomes(
    df: pd.DataFrame,
    mass_target_col: str = "m_target",
    mass_projectile_col: str = "m_impactor",
    radius_target_col: str = "r_target",
    radius_projectile_col: str = "r_impactor",
    velocity_col: str = "v_impact",
    impact_parameter_col: str = "b",
    config: Optional[EDACMConfig] = None,
    require_target_larger: bool = True,
    prefix: str = "",
) -> pd.DataFrame:
    """Add EDACM collision-outcome diagnostics to an event table.

    The input table is not modified in place. Required columns are target and
    projectile mass, target and projectile radius, impact velocity, and impact
    parameter. All quantities are expected to be in SI units.
    """

    required_columns = [
        mass_target_col,
        mass_projectile_col,
        radius_target_col,
        radius_projectile_col,
        velocity_col,
        impact_parameter_col,
    ]
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        raise KeyError(f"Missing required collision-outcome columns: {missing}")

    outcomes = classify_edacm_outcomes(
        df[mass_target_col].to_numpy(),
        df[mass_projectile_col].to_numpy(),
        df[radius_target_col].to_numpy(),
        df[radius_projectile_col].to_numpy(),
        df[velocity_col].to_numpy(),
        df[impact_parameter_col].to_numpy(),
        config=config,
        require_target_larger=require_target_larger,
    )

    out = df.copy()
    for key, values in outcomes.items():
        out[f"{prefix}{key}"] = values
    return out

