"""Softened axisymmetric gapped-disk gravity for Fig. 2 experiments.

This is a numerical ring-potential approximation, NOT the thin-disk Ward
series in Zheng et al. (2016), eqs. (3)-(4). Its softening and outer edge are
explicit configuration parameters whose effects on resonance positions must
be checked before claiming a quantitative Fig. 2 reproduction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from paper_fig2_disk import DiskProfile
from paper_fig2_drag import AU_CM, MSUN_G
from scipy.interpolate import CubicSpline
from scipy.special import ellipk


@dataclass(frozen=True)
class GapGravityConfig:
    inner_edge_au: float = 0.05
    gap_inner_au: float = 4.5
    gap_outer_au: float = 11.0
    outer_edge_au: float = 50.0
    softening_h: float = 0.6
    n_ring_per_region: int = 320
    table_min_au: float = 0.3
    table_max_au: float = 13.0
    n_table: int = 768

    def __post_init__(self):
        radii = (self.inner_edge_au, self.gap_inner_au,
                 self.gap_outer_au, self.outer_edge_au)
        if not all(np.isfinite(r) for r in radii) or not (
            0 < radii[0] < radii[1] < radii[2] < radii[3]
        ):
            raise ValueError("require ordered positive disk edges")
        if not (np.isfinite(self.softening_h) and self.softening_h > 0):
            raise ValueError("softening_h must be positive")
        if not (0 < self.table_min_au < self.table_max_au):
            raise ValueError("invalid force-table interval")
        if self.n_ring_per_region < 16 or self.n_table < 32:
            raise ValueError("force-table and ring resolutions are too small")


class PaperGapDiskGravity:
    """Gapped disk potential and radial force in G=Msun=AU=1 units.

    At t=0, each ring has surface density from ``disk``. The entire potential
    and force scale by ``disk.sigma(1,t)/disk.sigma(1,0)``. This is valid for
    fixed-shape exponential depletion and is checked at construction time.
    The potential is tabulated and differentiated with a cubic spline.
    """

    def __init__(self, disk: DiskProfile, config: GapGravityConfig | None = None):
        self.disk = disk
        self.config = config or GapGravityConfig()
        cfg = self.config
        self._ring_r, self._ring_weight = self._build_rings()
        self._softening = cfg.softening_h * np.asarray(
            disk.scale_height_au(self._ring_r, 0.0), dtype=float)
        self.table_r = np.linspace(cfg.table_min_au, cfg.table_max_au,
                                   cfg.n_table)
        self.table_phi = self._potential_reference(self.table_r)
        self._spline = CubicSpline(self.table_r, self.table_phi,
                                   bc_type="natural")
        self._sigma0 = float(disk.sigma_cgs(1.0, 0.0))
        if self._sigma0 <= 0:
            raise ValueError("disk surface density at 1 AU must be positive")
        sigma_initial = np.asarray(disk.sigma_cgs(self._ring_r, 0.0), dtype=float)
        if not np.all(np.isfinite(sigma_initial)) or np.any(sigma_initial <= 0):
            raise ValueError("disk ring surface density must be finite and positive")
        for t_probe in (0.5, 2.0):
            ratio = np.asarray(disk.sigma_cgs(self._ring_r, t_probe), dtype=float) / sigma_initial
            if not np.all(np.isfinite(ratio)) or not np.allclose(ratio, ratio[0], rtol=1e-10):
                raise ValueError("disk depletion must be separable in radius and time")

    def _build_rings(self):
        cfg = self.config
        nodes, weights = np.polynomial.legendre.leggauss(cfg.n_ring_per_region)
        pieces = []
        for lo, hi in ((cfg.inner_edge_au, cfg.gap_inner_au),
                       (cfg.gap_outer_au, cfg.outer_edge_au)):
            radius = 0.5*(hi-lo)*nodes + 0.5*(hi+lo)
            dr = 0.5*(hi-lo)*weights
            sigma_code = (np.asarray(self.disk.sigma_cgs(radius, 0), dtype=float)
                          * AU_CM**2 / MSUN_G)
            pieces.append((radius, sigma_code*radius*dr))
        return (np.concatenate([p[0] for p in pieces]),
                np.concatenate([p[1] for p in pieces]))

    def _potential_reference(self, radii):
        r = np.atleast_1d(np.asarray(radii, dtype=float))[:, None]
        s = self._ring_r[None, :]
        d2 = (r+s)**2 + self._softening[None, :]**2
        m = np.minimum(4*r*s/d2, 1-1e-14)
        # Angular integral of the softened ring potential is 4K(m)/sqrt(d2).
        out = -4 * np.sum(self._ring_weight[None, :]
                          * ellipk(m) / np.sqrt(d2), axis=1)
        return out

    def _checked_radius(self, r_au):
        r = np.asarray(r_au, dtype=float)
        if not np.all(np.isfinite(r)) or np.any(r < self.table_r[0]) or np.any(r > self.table_r[-1]):
            raise ValueError("radius outside calibrated disk-force table")
        return r

    def depletion_factor(self, t_myr):
        t = np.asarray(t_myr, dtype=float)
        if not np.all(np.isfinite(t)) or np.any(t < 0):
            raise ValueError("time must be finite and nonnegative")
        return np.asarray(self.disk.sigma_cgs(1.0, t), dtype=float) / self._sigma0

    def potential(self, r_au, t_myr=0.0):
        r = self._checked_radius(r_au)
        return self._spline(r) * self.depletion_factor(t_myr)

    def radial_acceleration(self, r_au, t_myr=0.0):
        """Outward-positive radial acceleration in AU/code_time**2."""
        r = self._checked_radius(r_au)
        return -self._spline(r, 1) * self.depletion_factor(t_myr)

    def radial_derivative(self, r_au, t_myr=0.0):
        """d(radial acceleration)/dr in 1/code_time**2."""
        r = self._checked_radius(r_au)
        return -self._spline(r, 2) * self.depletion_factor(t_myr)

    def acceleration(self, position_au, t_myr=0.0):
        """Cartesian disk acceleration in AU/code_time**2."""
        x = np.asarray(position_au, dtype=float)
        if x.shape[-1] not in (2, 3):
            raise ValueError("position must have 2 or 3 components")
        r = self._checked_radius(np.hypot(x[..., 0], x[..., 1]))
        fr = self.radial_acceleration(r, t_myr)
        out = np.zeros_like(x)
        out[..., 0] = fr*x[..., 0]/r
        out[..., 1] = fr*x[..., 1]/r
        return out
