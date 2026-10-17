"""Explicit initial conditions and force switches for Fig. 2 experiments."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Fig2SystemConfig:
    star_mass_msun: float = 1.0
    jupiter_mass_msun: float = 9.5458e-4
    saturn_mass_msun: float = 2.858e-4
    jupiter_a_au: float = 5.2
    saturn_a_au: float = 9.58
    jupiter_e: float = 0.05
    saturn_e: float = 0.0
    jupiter_omega_rad: float = 0.0
    saturn_omega_rad: float = 0.0
    jupiter_mean_anomaly_rad: float = 0.0
    saturn_mean_anomaly_rad: float = 0.0
    particle_a_au: tuple[float, ...] = (3.5, 3.2, 2.7, 2.3)
    particle_e: float | tuple[float, ...] = 0.01
    particle_omega_rad: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0)
    particle_mean_anomaly_rad: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0)
    use_planets: bool = True
    use_disk_gravity: bool = True
    use_drag: bool = True

    def __post_init__(self):
        positives = (self.star_mass_msun, self.jupiter_mass_msun,
                     self.saturn_mass_msun, self.jupiter_a_au,
                     self.saturn_a_au, *self.particle_a_au)
        if not all(np.isfinite(x) and x > 0 for x in positives):
            raise ValueError("masses and semimajor axes must be positive and finite")
        if not (0 <= self.jupiter_e < 1 and 0 <= self.saturn_e < 1):
            raise ValueError("eccentricities must be in [0,1)")
        n = len(self.particle_a_au)
        eccentricities = self.particle_e_values
        if len(eccentricities) != n or not all(0 <= e < 1 for e in eccentricities):
            raise ValueError("particle eccentricities must have length N and lie in [0,1)")
        if n == 0 or len(self.particle_omega_rad) != n or len(self.particle_mean_anomaly_rad) != n:
            raise ValueError("one omega and mean anomaly required per test particle")
        angles = (self.jupiter_omega_rad, self.saturn_omega_rad,
                  self.jupiter_mean_anomaly_rad, self.saturn_mean_anomaly_rad,
                  *self.particle_omega_rad, *self.particle_mean_anomaly_rad)
        if not all(np.isfinite(x) for x in angles):
            raise ValueError("orbital angles must be finite")

    @property
    def particle_e_values(self) -> tuple[float, ...]:
        """One eccentricity per test particle, expanding a scalar if needed."""
        if np.isscalar(self.particle_e):
            return (float(self.particle_e),) * len(self.particle_a_au)
        return tuple(float(e) for e in self.particle_e)


@dataclass(frozen=True)
class IntegrationConfig:
    integrator: str = "ias15"
    initial_step_years: float = 0.1
    force_backend: str = "python"

    def __post_init__(self):
        if self.integrator not in ("ias15", "whfast"):
            raise ValueError("integrator must be ias15 or whfast")
        if self.force_backend not in ("python", "compiled"):
            raise ValueError("force_backend must be python or compiled")
        if not np.isfinite(self.initial_step_years) or self.initial_step_years <= 0:
            raise ValueError("initial_step_years must be positive")

    def estimated_fixed_steps(self, t_end_myr: float):
        return int(np.ceil(t_end_myr * 1e6 / self.initial_step_years))
