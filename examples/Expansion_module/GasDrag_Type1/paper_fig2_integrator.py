"""REBOUND experimental backend for gravity + disk force + aerodynamic drag.

Short-run validation backend. A 10 Myr run with Python additional-force calls
can be prohibitively slow; use benchmark results before attempting it.
"""

from __future__ import annotations

import hashlib
import inspect
import json
from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path

import numpy as np
from paper_fig2_disk import DiskProfile
from paper_fig2_drag import MYR_PER_CODE, AerodynamicDrag
from paper_fig2_events import (
    EVENT_NAMES,
    ParticleBoundaryConfig,
    ParticleEvent,
    classify_state,
)
from paper_fig2_gravity import PaperGapDiskGravity
from paper_fig2_system import Fig2SystemConfig, IntegrationConfig


@dataclass
class Fig2Run:
    time_myr: np.ndarray
    names: tuple[str, ...]
    position_au: np.ndarray
    velocity_code: np.ndarray
    a_au: np.ndarray
    e: np.ndarray
    omega_rad: np.ndarray
    active: np.ndarray
    events: tuple[ParticleEvent, ...]
    metadata: dict

    @property
    def q_au(self):
        return self.a_au * (1 - self.e)

    @property
    def Q_au(self):
        return self.a_au * (1 + self.e)

    def save_npz(self, path):
        """Save all numerical tracks and a JSON metadata string."""
        dest = Path(path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(dest, time_myr=self.time_myr,
                            names=np.asarray(self.names),
                            position_au=self.position_au,
                            velocity_code=self.velocity_code,
                            a_au=self.a_au, e=self.e,
                            omega_rad=self.omega_rad,
                            active=self.active,
                            events_json=json.dumps(
                                [event.to_dict() for event in self.events], sort_keys=True),
                            metadata_json=json.dumps(self.metadata, sort_keys=True))

    @classmethod
    def load_npz(cls, path):
        """Reload a saved track without executing an N-body integration."""
        with np.load(path, allow_pickle=False) as saved:
            a_au = saved["a_au"]
            active = (saved["active"] if "active" in saved
                      else np.isfinite(a_au))
            events = (tuple(ParticleEvent.from_dict(item) for item in
                            json.loads(str(saved["events_json"])))
                      if "events_json" in saved else ())
            metadata = json.loads(str(saved["metadata_json"]))
            if "active" not in saved:
                metadata = {**metadata, "event_schema_version": 0}
            return cls(
                time_myr=saved["time_myr"],
                names=tuple(str(x) for x in saved["names"]),
                position_au=saved["position_au"],
                velocity_code=saved["velocity_code"],
                a_au=a_au,
                e=saved["e"],
                omega_rad=saved["omega_rad"],
                active=active, events=events, metadata=metadata,
            )


def integrate_fig2(system: Fig2SystemConfig, disk: DiskProfile,
                   gravity: PaperGapDiskGravity | None,
                   drag: AerodynamicDrag | None,
                   output_times_myr,
                   integration: IntegrationConfig | None = None,
                   boundaries: ParticleBoundaryConfig | None = None) -> Fig2Run:
    """Integrate the specified experiment at requested times.

    REBOUND integration time is code time, not Myr. All saved positions are
    heliocentric; velocities are relative to the Sun. Test particles have
    zero gravitational mass. This function does not call orbit averaging.
    """
    try:
        import rebound
    except ImportError as exc:
        raise ImportError("REBOUND is needed for direct Fig. 2 experiments") from exc

    if system.use_disk_gravity and gravity is None:
        raise ValueError("disk gravity is enabled but no gravity model was given")
    if system.use_drag and drag is None:
        raise ValueError("drag is enabled but no drag model was given")
    cfg = integration or IntegrationConfig()
    times = np.asarray(output_times_myr, dtype=float)
    if times.ndim != 1 or len(times) == 0 or not np.all(np.isfinite(times)) or np.any(times < 0):
        raise ValueError("output times must be a finite, nonnegative 1D array")
    if np.any(np.diff(times) <= 0):
        raise ValueError("output times must be strictly increasing")

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = cfg.integrator
    sim.dt = cfg.initial_step_years / (1e6 * MYR_PER_CODE)
    sim.add(m=system.star_mass_msun, hash="Sun")
    primary = sim.particles[0]
    names = []
    if system.use_planets:
        sim.add(m=system.jupiter_mass_msun, a=system.jupiter_a_au,
                e=system.jupiter_e, omega=system.jupiter_omega_rad,
                M=system.jupiter_mean_anomaly_rad, hash="Jupiter",
                primary=primary)
        sim.add(m=system.saturn_mass_msun, a=system.saturn_a_au,
                e=system.saturn_e, omega=system.saturn_omega_rad,
                M=system.saturn_mean_anomaly_rad, hash="Saturn",
                primary=primary)
        names += ["Jupiter", "Saturn"]
    sim.N_active = sim.N
    for j, (a, eccentricity, omega, mean_anomaly) in enumerate(zip(
            system.particle_a_au, system.particle_e_values, system.particle_omega_rad,
            system.particle_mean_anomaly_rad)):
        name = f"particle_{j}"
        sim.add(m=0, a=a, e=eccentricity,
                omega=omega, M=mean_anomaly, hash=name,
                primary=primary)
        names.append(name)
    sim.testparticle_type = 0
    sim.move_to_com()
    n_planets = 2 if system.use_planets else 0
    first_test = 1 + n_planets
    if boundaries is not None and gravity is not None and system.use_disk_gravity:
        boundaries.validate_force_table(gravity.table_r[0], gravity.table_r[-1])

    original_names = tuple(names)
    current_names = list(names)
    recorded_events = []
    pending_python_events = []
    inactive_python_names = set()
    compiled_force = None

    def install_forces():
        nonlocal compiled_force, pending_python_events, inactive_python_names
        pending_python_events = []
        inactive_python_names = set()
        uses_callback = system.use_disk_gravity or system.use_drag or boundaries is not None
        if not uses_callback:
            compiled_force = None
            return
        if cfg.force_backend == "compiled":
            from paper_fig2_compiled import CompiledFig2Force

            compiled_force = CompiledFig2Force(
                disk, gravity, drag, first_test, system.use_disk_gravity,
                system.use_drag, sim.N, boundaries)
            sim.additional_forces = compiled_force.function
        else:
            compiled_force = None

            def additional_forces(reb_sim):
                current = reb_sim.contents
                particles = current.particles
                sun = particles[0]
                t_myr = current.t * MYR_PER_CODE
                for particle_index in range(1, current.N):
                    p = particles[particle_index]
                    name = current_names[particle_index-1]
                    xyz = np.array([p.x-sun.x, p.y-sun.y, p.z-sun.z])
                    velocity_xyz = np.array([p.vx-sun.vx, p.vy-sun.vy, p.vz-sun.vz])
                    if boundaries is not None and particle_index >= first_test:
                        if name in inactive_python_names:
                            continue
                        result = classify_state(
                            xyz, velocity_xyz, system.star_mass_msun, boundaries)
                        if result is not None:
                            code, radius, radial_velocity, energy = result
                            pending_python_events.append(ParticleEvent(
                                particle_index=original_names.index(name),
                                particle_name=name, event_code=code,
                                event_type=EVENT_NAMES[code], time_myr=t_myr,
                                radius_au=radius,
                                position_au=tuple(float(x) for x in xyz),
                                velocity_code=tuple(float(x) for x in velocity_xyz),
                                radial_velocity_code=radial_velocity,
                                kepler_energy_code=energy,
                            ))
                            inactive_python_names.add(name)
                            continue
                    acc = np.zeros(2)
                    if system.use_disk_gravity:
                        acc += gravity.acceleration(xyz[:2], t_myr)
                    if system.use_drag and particle_index >= first_test:
                        acc += drag.acceleration(xyz[:2], velocity_xyz[:2], t_myr, disk)
                    p.ax += acc[0]
                    p.ay += acc[1]

            sim.additional_forces = additional_forces
        if system.use_drag:
            sim.force_is_velocity_dependent = 1

    def collect_events_and_remove():
        nonlocal current_names
        if compiled_force is not None:
            compiled_force.check_domain(current_names)
            new_events = compiled_force.events(current_names, original_names)
        else:
            new_events = list(pending_python_events)
        if not new_events:
            return
        known = {event.particle_name for event in recorded_events}
        new_events = [event for event in new_events if event.particle_name not in known]
        if not new_events:
            return
        recorded_events.extend(new_events)
        removal = sorted(
            (current_names.index(event.particle_name)+1, event.particle_name)
            for event in new_events)
        for particle_index, name in reversed(removal):
            sim.remove(index=particle_index, keep_sorted=True)
            current_names.pop(particle_index-1)
        install_forces()

    install_forces()

    count = len(original_names)
    position = np.full((len(times), count, 3), np.nan)
    velocity = np.empty_like(position)
    velocity.fill(np.nan)
    a_arr = np.full((len(times), count), np.nan)
    e_arr = np.full_like(a_arr, np.nan)
    omega_arr = np.full_like(a_arr, np.nan)
    active_arr = np.zeros((len(times), count), dtype=bool)
    actual_times = np.empty(len(times))

    for k, t_myr in enumerate(times):
        target = t_myr / MYR_PER_CODE
        interval = ((boundaries.event_check_interval_years /
                     (1e6*MYR_PER_CODE)) if boundaries is not None else np.inf)
        while sim.t < target:
            segment_end = min(target, sim.t+interval)
            sim.integrate(segment_end,
                          exact_finish_time=0 if cfg.integrator == "whfast" else 1)
            collect_events_and_remove()
        collect_events_and_remove()
        actual_times[k] = sim.t * MYR_PER_CODE
        sun = sim.particles[0]
        for j, name in enumerate(original_names):
            if name not in current_names:
                continue
            p = sim.particles[name]
            active_arr[k, j] = True
            position[k, j] = (p.x-sun.x, p.y-sun.y, p.z-sun.z)
            velocity[k, j] = (p.vx-sun.vx, p.vy-sun.vy, p.vz-sun.vz)
            orbit = p.orbit(primary=sun)
            a_arr[k, j], e_arr[k, j], omega_arr[k, j] = (
                orbit.a, orbit.e, orbit.omega)

    disk_description = {"type": type(disk).__name__}
    if is_dataclass(disk):
        # The existing-disk adapter contains a model object; serialize its
        # configuration instead of a non-JSON-safe model instance.
        if hasattr(disk, "disk") and hasattr(disk.disk, "cfg"):
            disk_description["config"] = asdict(disk.disk.cfg)
        else:
            disk_description["config"] = asdict(disk)
    settings = {
        "system": asdict(system), "integration": asdict(cfg),
        "units": "AU, Msun, G=1, code_time=year/(2pi)",
        "disk": disk_description,
        "gravity": asdict(gravity.config) if gravity else None,
        "boundaries": asdict(boundaries) if boundaries else None,
        "event_schema_version": 1,
        "drag": ({"body": asdict(drag.planetesimal),
                  "flow": asdict(drag.flow), "law": drag.law,
                  "drag_coefficient": drag.cd} if drag else None),
        "rebound_version": rebound.__version__,
    }
    config_bytes = json.dumps(settings, sort_keys=True).encode("utf-8")
    settings["config_sha256"] = hashlib.sha256(config_bytes).hexdigest()
    source_hash = hashlib.sha256()
    for filename in ("paper_fig2_disk.py", "paper_fig2_drag.py",
                     "paper_fig2_gravity.py", "paper_fig2_system.py",
                     "paper_fig2_events.py",
                     "paper_fig2_integrator.py", "paper_fig2_experiment.py",
                     "paper_fig2_compiled.py",
                     "paper_fig2_force.c"):
        source_hash.update(filename.encode("utf-8"))
        source_hash.update(Path(__file__).with_name(filename).read_bytes())
    if hasattr(disk, "disk"):
        disk_source = Path(inspect.getfile(type(disk.disk)))
        source_hash.update(disk_source.read_bytes())
    settings["source_sha256"] = source_hash.hexdigest()
    settings["event_counts"] = {
        event_type: sum(event.event_type == event_type for event in recorded_events)
        for event_type in EVENT_NAMES.values()
    }
    return Fig2Run(actual_times, original_names, position, velocity,
                   a_arr, e_arr, omega_arr, active_arr,
                   tuple(recorded_events), settings)
