"""Interchangeable local gas-disk profiles for the Fig. 2 drag baseline.

Public radii and heights are AU, times are Myr, and surface densities are
g/cm**2. This module supplies local gas properties only; it has no disk
gravity, gap, Type-I torque, or resonance model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

import numpy as np

if TYPE_CHECKING:
    from accrediff.gas_migration import Gas_DiskModel, Gas_UnitSystem


def _result(value):
    arr = np.asarray(value, dtype=float)
    return float(arr) if arr.ndim == 0 else arr


def _radius(r_au):
    r = np.asarray(r_au, dtype=float)
    if not np.all(np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("r_au must be finite and positive")
    return r


def _time(t_myr):
    t = np.asarray(t_myr, dtype=float)
    if not np.all(np.isfinite(t)) or np.any(t < 0):
        raise ValueError("t_myr must be finite and nonnegative")
    return t


class DiskProfile(Protocol):
    """Minimum local disk interface used by gas flow and drag."""

    def sigma_cgs(self, r_au, t_myr): ...

    def scale_height_au(self, r_au, t_myr): ...

    def sigma_slope_p(self, r_au): ...


@dataclass(frozen=True)
class PaperFig2Disk:
    """Zheng et al. (2016), eqs. (1)-(2), default model without gap gravity."""

    sigma0_cgs: float = 1700.0
    surface_density_exponent: float = 1.5
    h1_au: float = 0.025
    height_exponent: float = 1.25
    tau_decay_myr: float = 1.0

    def __post_init__(self):
        vals = (self.sigma0_cgs, self.surface_density_exponent,
                self.h1_au, self.height_exponent, self.tau_decay_myr)
        if not all(np.isfinite(v) for v in vals):
            raise ValueError("disk parameters must be finite")
        if self.sigma0_cgs <= 0 or self.h1_au <= 0 or self.tau_decay_myr <= 0:
            raise ValueError("sigma0, h1 and tau_decay must be positive")

    def sigma_cgs(self, r_au, t_myr):
        r, t = _radius(r_au), _time(t_myr)
        out = self.sigma0_cgs * r**(-self.surface_density_exponent)
        out = out * np.exp(-t / self.tau_decay_myr)
        return _result(out)

    def scale_height_au(self, r_au, t_myr):
        r = _radius(r_au)
        _time(t_myr)
        return _result(self.h1_au * r**self.height_exponent)

    def sigma_slope_p(self, r_au):
        r = _radius(r_au)
        return _result(np.full_like(r, self.surface_density_exponent))


@dataclass(frozen=True)
class GasMigrationDiskAdapter:
    """Expose an existing Gas_DiskModel with the same physical-unit interface."""

    disk: Gas_DiskModel
    units: Gas_UnitSystem

    def sigma_cgs(self, r_au, t_myr):
        r, t = _radius(r_au), _time(t_myr)
        return _result(self.disk.sigma_cgs(r, self.units.myr_to_code(t)))

    def scale_height_au(self, r_au, t_myr):
        r = _radius(r_au)
        _time(t_myr)
        return _result(self.disk.scale_height(r))

    def sigma_slope_p(self, r_au):
        r = _radius(r_au)
        return _result(self.disk.sigma_slope_p(r))
