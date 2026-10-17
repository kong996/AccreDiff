"""Orbit average of the same instantaneous drag used by future N-body work."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from paper_fig2_disk import DiskProfile
from paper_fig2_drag import MYR_PER_CODE, AerodynamicDrag


@dataclass(frozen=True)
class AveragedRates:
    da_dt_au_per_myr: float
    de_dt_per_myr: float


def orbit_averaged_rates(a_au: float, e: float, t_myr: float,
                         disk: DiskProfile, drag: AerodynamicDrag,
                         n_samples: int = 128) -> AveragedRates:
    """Mean-anomaly average of planar Kepler elements under instantaneous drag.

    The eccentricity-vector derivative avoids division by e near a circle.
    The result is a secular diagnostic; it does not model planetary resonances.
    """
    if not np.isfinite(a_au) or a_au <= 0 or not np.isfinite(e) or not 0 <= e < 1:
        raise ValueError("require a>0 and 0<=e<1")
    if not isinstance(n_samples, int) or n_samples < 16:
        raise ValueError("n_samples must be an integer >=16")
    mu = drag.flow.mu_code
    E = 2*np.pi*(np.arange(n_samples) + 0.5)/n_samples
    ce, se = np.cos(E), np.sin(E)
    root = np.sqrt(1 - e*e)
    weight = 1 - e*ce  # dM/dE
    vscale = np.sqrt(mu/a_au)/weight
    x = np.stack((a_au*(ce-e), a_au*root*se), axis=-1)
    v = np.stack((-vscale*se, vscale*root*ce), axis=-1)
    f = drag.acceleration(x, v, t_myr, disk)
    h = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
    h_dot = x[:, 0]*f[:, 1] - x[:, 1]*f[:, 0]
    energy_dot = np.sum(v*f, axis=-1)
    da_code = 2*a_au*a_au*energy_dot/mu
    de_code = (f[:, 1]*h + v[:, 1]*h_dot)/mu
    return AveragedRates(
        float(np.mean(da_code*weight)/MYR_PER_CODE),
        float(np.mean(de_code*weight)/MYR_PER_CODE),
    )
