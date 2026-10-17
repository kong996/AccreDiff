"""Per-particle absorbing boundary events for Fig. 2 experiments."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

INNER_DOMAIN_EXIT = 1
OUTER_DOMAIN_EXIT_BOUND = 2
EJECTED = 3
NONFINITE = 4

EVENT_NAMES = {
    INNER_DOMAIN_EXIT: "inner_domain_exit",
    OUTER_DOMAIN_EXIT_BOUND: "outer_domain_exit_bound",
    EJECTED: "ejected",
    NONFINITE: "nonfinite",
}


@dataclass(frozen=True)
class ParticleBoundaryConfig:
    inner_termination_au: float = 0.3
    outer_termination_au: float = 13.0
    event_check_interval_years: float = 100.0
    require_positive_energy_for_ejection: bool = True

    def __post_init__(self):
        values = (self.inner_termination_au, self.outer_termination_au,
                  self.event_check_interval_years)
        if not all(np.isfinite(x) and x > 0 for x in values):
            raise ValueError("particle boundaries and check interval must be finite and positive")
        if self.inner_termination_au >= self.outer_termination_au:
            raise ValueError("inner particle boundary must be below outer boundary")

    def validate_force_table(self, table_min_au: float, table_max_au: float):
        if not (table_min_au < self.inner_termination_au <
                self.outer_termination_au < table_max_au):
            raise ValueError(
                "particle termination boundaries must lie strictly inside the disk-force table"
            )


@dataclass(frozen=True)
class ParticleEvent:
    particle_index: int
    particle_name: str
    event_code: int
    event_type: str
    time_myr: float
    radius_au: float
    position_au: tuple[float, float, float]
    velocity_code: tuple[float, float, float]
    radial_velocity_code: float
    kepler_energy_code: float

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, value):
        item = dict(value)
        item["position_au"] = tuple(item["position_au"])
        item["velocity_code"] = tuple(item["velocity_code"])
        return cls(**item)


def classify_state(position_au, velocity_code, star_mass_msun: float,
                   config: ParticleBoundaryConfig):
    """Return event diagnostics, or ``None`` while the test particle is active."""
    position = np.asarray(position_au, dtype=float)
    velocity = np.asarray(velocity_code, dtype=float)
    if position.shape != (3,) or velocity.shape != (3,):
        raise ValueError("event state vectors must have three components")
    if not np.all(np.isfinite(position)) or not np.all(np.isfinite(velocity)):
        return NONFINITE, np.nan, np.nan, np.nan
    radius = float(np.linalg.norm(position))
    if not np.isfinite(radius) or radius <= 0:
        return NONFINITE, radius, np.nan, np.nan
    radial_velocity = float(np.dot(position, velocity) / radius)
    energy = float(0.5*np.dot(velocity, velocity) - star_mass_msun/radius)
    if radius <= config.inner_termination_au:
        return INNER_DOMAIN_EXIT, radius, radial_velocity, energy
    if radius >= config.outer_termination_au:
        escaped = radial_velocity > 0 and (
            energy >= 0 or not config.require_positive_energy_for_ejection)
        code = EJECTED if escaped else OUTER_DOMAIN_EXIT_BOUND
        return code, radius, radial_velocity, energy
    return None
