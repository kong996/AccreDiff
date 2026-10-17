"""Instantaneous aerodynamic drag for the coplanar Fig. 2 baseline.

Position is AU, velocity is AU per GENGA-like code time, and the returned
acceleration is AU per code time squared. One code time is 1/(2*pi) year.
The paper and GENGA coefficient conventions are intentionally distinct.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from paper_fig2_disk import DiskProfile

AU_CM = 1.495978707e13
MSUN_G = 1.98847e33
MEARTH_TO_MSUN = 3.003489e-6
MEARTH_G = MEARTH_TO_MSUN * MSUN_G
YEAR_SECONDS = 365.25 * 86400.0
CODE_TIME_SECONDS = YEAR_SECONDS / (2 * np.pi)
MYR_PER_CODE = 1 / (2 * np.pi * 1e6)


@dataclass(frozen=True)
class Planetesimal:
    radius_km: float = 10.0
    solid_density_cgs: float = 1.0

    def __post_init__(self):
        if not (np.isfinite(self.radius_km) and self.radius_km > 0 and
                np.isfinite(self.solid_density_cgs) and self.solid_density_cgs > 0):
            raise ValueError("radius and solid density must be finite and positive")

    @property
    def radius_cm(self):
        return self.radius_km * 1e5

    @property
    def mass_g(self):
        return (4 * np.pi / 3) * self.solid_density_cgs * self.radius_cm**3

    @property
    def mass_mearth(self):
        return self.mass_g / MEARTH_G


@dataclass(frozen=True)
class GasFlowPolicy:
    """Midplane density and headwind choices, independent of disk structure.

    The paper's eta=(c_s/4v_K)^2 uses the closure c_s=H*Omega_K.
    genga_style uses eta=0.5*(p+2-beta)*(H/r)^2. For the gas.cu
    default beta=0.25, this is its p+1.75 expression.
    """

    convention: Literal["paper_2016", "genga_style"] = "paper_2016"
    beta: float = 0.25
    mu_code: float = 1.0

    def __post_init__(self):
        if self.convention not in ("paper_2016", "genga_style"):
            raise ValueError("unknown gas-flow convention")
        if not np.isfinite(self.beta) or not np.isfinite(self.mu_code) or self.mu_code <= 0:
            raise ValueError("beta must be finite and mu_code positive")

    def density_cgs(self, disk: DiskProfile, r_au, t_myr, z_au=0.0):
        r, z = np.asarray(r_au, dtype=float), np.asarray(z_au, dtype=float)
        sigma = np.asarray(disk.sigma_cgs(r, t_myr), dtype=float)
        height = np.asarray(disk.scale_height_au(r, t_myr), dtype=float)
        if self.convention == "paper_2016":
            if np.any(z != 0):
                raise ValueError("paper_2016 flow is defined for coplanar z=0 only")
            rho = sigma / (height * AU_CM)
        else:
            rho = sigma / (np.sqrt(2*np.pi) * height * AU_CM)
            rho = rho * np.exp(-0.5 * (z / height)**2)
        return float(rho) if rho.ndim == 0 else rho

    def eta(self, disk: DiskProfile, r_au, t_myr):
        r = np.asarray(r_au, dtype=float)
        height = np.asarray(disk.scale_height_au(r, t_myr), dtype=float)
        h = height / r
        if self.convention == "paper_2016":
            eta = (h / 4)**2
        else:
            p = np.asarray(disk.sigma_slope_p(r), dtype=float)
            eta = 0.5 * (p + 2 - self.beta) * h**2
        return float(eta) if np.ndim(eta) == 0 else eta

    def velocity_code(self, position_au, t_myr, disk: DiskProfile):
        """Return azimuthal gas velocity in AU per code time."""
        x = np.asarray(position_au, dtype=float)
        if x.shape[-1] not in (2, 3):
            raise ValueError("position must end in 2 or 3 Cartesian components")
        r = np.hypot(x[..., 0], x[..., 1])
        if np.any(r <= 0) or not np.all(np.isfinite(r)):
            raise ValueError("cylindrical radius must be finite and positive")
        vk = np.sqrt(self.mu_code / r)
        vphi = vk * (1 - self.eta(disk, r, t_myr))
        out = np.zeros_like(x)
        out[..., 0] = -vphi * x[..., 1] / r
        out[..., 1] = vphi * x[..., 0] / r
        return out


@dataclass(frozen=True)
class AerodynamicDrag:
    planetesimal: Planetesimal
    flow: GasFlowPolicy
    law: Literal["paper_2016", "genga_style"] = "paper_2016"
    drag_coefficient: float | None = None

    def __post_init__(self):
        if self.law not in ("paper_2016", "genga_style"):
            raise ValueError("unknown drag law")
        if self.drag_coefficient is not None and (
            not np.isfinite(self.drag_coefficient) or self.drag_coefficient <= 0
        ):
            raise ValueError("drag coefficient must be finite and positive")

    @property
    def cd(self):
        if self.drag_coefficient is not None:
            return self.drag_coefficient
        return 0.165 if self.law == "paper_2016" else 0.5

    def acceleration(self, position_au, velocity_code, t_myr, disk: DiskProfile):
        """Return drag acceleration in AU/(code time)**2.

        Paper: -Cd*rho/(rho_p*R)*|u|u (Zheng et al. eqs. 5-6).
        GENGA: -3Cd*rho/(8rho_p*R)*|u|u (gas.cu).
        These differ by a factor 3/8 for the same Cd, rho and velocity.
        """
        x = np.asarray(position_au, dtype=float)
        v = np.asarray(velocity_code, dtype=float)
        if x.shape != v.shape or x.shape[-1] not in (2, 3):
            raise ValueError("position and velocity shapes must match with 2 or 3 components")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(v)):
            raise ValueError("position and velocity must be finite")
        r = np.hypot(x[..., 0], x[..., 1])
        z = 0.0 if x.shape[-1] == 2 else x[..., 2]
        rho = np.asarray(self.flow.density_cgs(disk, r, t_myr, z), dtype=float)
        u = v - self.flow.velocity_code(x, t_myr, disk)
        u_cgs = u * (AU_CM / CODE_TIME_SECONDS)
        factor = 1.0 if self.law == "paper_2016" else 3.0 / 8.0
        coeff = factor * self.cd * rho / (
            self.planetesimal.solid_density_cgs * self.planetesimal.radius_cm)
        a_cgs = -coeff[..., None] * np.linalg.norm(u_cgs, axis=-1)[..., None] * u_cgs
        return a_cgs * CODE_TIME_SECONDS**2 / AU_CM
