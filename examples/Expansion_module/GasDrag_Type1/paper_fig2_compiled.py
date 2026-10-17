"""Build and bind the compiled REBOUND force for two gas-disk providers."""

from __future__ import annotations

import ctypes
import hashlib
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from paper_fig2_disk import GasMigrationDiskAdapter, PaperFig2Disk
from paper_fig2_drag import AU_CM, MSUN_G, MYR_PER_CODE, AerodynamicDrag
from paper_fig2_events import EVENT_NAMES, ParticleBoundaryConfig, ParticleEvent
from paper_fig2_gravity import PaperGapDiskGravity

DoublePtr = ctypes.POINTER(ctypes.c_double)
IntPtr = ctypes.POINTER(ctypes.c_int)


class DiskForceDomainError(RuntimeError):
    """A body left the radial domain covered by the compiled force table."""

    def __init__(self, *, count, particle_index, particle_name, kind,
                 time_myr, radius_au, kepler_energy, table_bounds_au,
                 offending_range_au):
        self.count = count
        self.particle_index = particle_index
        self.particle_name = particle_name
        self.kind = kind
        self.time_myr = time_myr
        self.radius_au = radius_au
        self.kepler_energy = kepler_energy
        self.table_bounds_au = table_bounds_au
        self.offending_range_au = offending_range_au
        reason = "radial table boundary" if kind == 1 else "nonfinite/nonpositive radius"
        orbit_state = ("unbound in the star-only two-body diagnostic" if kepler_energy >= 0
                       else "bound in the star-only two-body diagnostic")
        super().__init__(
            f"compiled disk force failed: first {reason} violation by {particle_name} "
            f"at t={time_myr:.9g} Myr, r={radius_au:.9g} AU; "
            f"Kepler energy={kepler_energy:.9g} ({orbit_state}); "
            f"valid table=[{table_bounds_au[0]:.9g}, {table_bounds_au[1]:.9g}] AU; "
            f"{count} force evaluations violated the domain "
            f"(finite offending range=[{offending_range_au[0]:.9g}, "
            f"{offending_range_au[1]:.9g}] AU)"
        )


class _Params(ctypes.Structure):
    _fields_ = [
        ("n_table", ctypes.c_int),
        ("table_r", DoublePtr),
        ("table_fr", DoublePtr),
        ("table_dfr", DoublePtr),
        ("tau_code", ctypes.c_double),
        ("sigma0_cgs", ctypes.c_double),
        ("alpha1", ctypes.c_double),
        ("alpha2", ctypes.c_double),
        ("h1_au", ctypes.c_double),
        ("height_exponent", ctypes.c_double),
        ("mu_code", ctypes.c_double),
        ("drag_coefficient", ctypes.c_double),
        ("solid_density_cgs", ctypes.c_double),
        ("radius_cm", ctypes.c_double),
        ("au_cm", ctypes.c_double),
        ("first_test", ctypes.c_int),
        ("use_disk_gravity", ctypes.c_int),
        ("use_drag", ctypes.c_int),
        ("use_boundaries", ctypes.c_int),
        ("inner_boundary", ctypes.c_double),
        ("outer_boundary", ctypes.c_double),
        ("require_positive_energy", ctypes.c_int),
        ("active", IntPtr),
        ("event_code", IntPtr),
        ("event_time", DoublePtr),
        ("event_radius", DoublePtr),
        ("event_position", DoublePtr),
        ("event_velocity", DoublePtr),
        ("event_radial_velocity", DoublePtr),
        ("event_kepler_energy", DoublePtr),
    ]


class CompiledFig2Force:
    """Own NumPy force tables and C library for one active simulation."""

    def __init__(self, disk: PaperFig2Disk | GasMigrationDiskAdapter,
                 gravity: PaperGapDiskGravity | None,
                 drag: AerodynamicDrag | None, first_test: int,
                 use_disk_gravity: bool, use_drag: bool, n_particles: int,
                 boundaries: ParticleBoundaryConfig | None = None):
        if isinstance(disk, PaperFig2Disk):
            tau_myr = disk.tau_decay_myr
            sigma0 = disk.sigma0_cgs
            alpha1, alpha2 = 0.0, disk.surface_density_exponent
            h1, height_exponent = disk.h1_au, disk.height_exponent
        elif isinstance(disk, GasMigrationDiskAdapter):
            config = disk.disk.cfg
            if config.G != 1 or config.Mstar != 1:
                raise ValueError("compiled force requires G=Mstar=1 code units")
            if (not np.isclose(config.AU_in_cm, AU_CM, rtol=1e-12) or
                    not np.isclose(config.Msun_in_g, MSUN_G, rtol=1e-12)):
                raise ValueError("disk physical units differ from drag/gravity units")
            if disk.units is not disk.disk.units:
                raise ValueError("adapter and Gas_DiskModel must use the same time units")
            tau_myr = config.tau_decay_myr
            sigma0 = config.Sigma0_cgs
            alpha1, alpha2 = config.a_alpha1, config.a_alpha2
            h1, height_exponent = config.def_h_1, 1 + config.beta
        else:
            raise TypeError("compiled force supports paper or GasMigrationDiskAdapter")
        if not np.all(np.isfinite((tau_myr, sigma0, alpha1, alpha2,
                                   h1, height_exponent))) or min(tau_myr, sigma0, h1) <= 0:
            raise ValueError("compiled disk parameters must be finite and positive where required")
        if use_disk_gravity and gravity is None:
            raise ValueError("disk gravity model required")
        if use_drag and (
            drag is None or drag.law != "paper_2016" or
            drag.flow.convention != "paper_2016"
        ):
            raise ValueError("compiled force supports paper_2016 drag and flow only")
        if gravity is not None and gravity.disk is not disk:
            raise ValueError("gravity and drag must share the same disk")

        import rebound

        source = Path(__file__).with_name("paper_fig2_force.c")
        header_dir = Path(rebound.__file__).resolve().parent
        key = hashlib.sha256(source.read_bytes() +
                             (header_dir / "rebound.h").read_bytes()).hexdigest()[:16]
        build_dir = Path(tempfile.gettempdir()) / "accrediff_fig2_compiled"
        build_dir.mkdir(parents=True, exist_ok=True)
        lib_path = build_dir / f"fig2_force_{key}.dylib"
        if not lib_path.exists():
            temp_path = build_dir / f"fig2_force_{key}.tmp.dylib"
            subprocess.run(["cc", "-O3", "-std=c11", "-fPIC", "-shared",
                            "-I", str(header_dir), str(source),
                            "-o", str(temp_path)], check=True)
            temp_path.replace(lib_path)
        self.lib = ctypes.CDLL(str(lib_path))
        self.lib.fig2_set_params.argtypes = [ctypes.POINTER(_Params)]
        self.lib.fig2_set_params.restype = None
        self.lib.fig2_get_domain_violations.restype = ctypes.c_int
        self.lib.fig2_get_first_violation_particle.restype = ctypes.c_int
        self.lib.fig2_get_first_violation_kind.restype = ctypes.c_int
        for getter in ("fig2_get_first_violation_time",
                       "fig2_get_first_violation_radius",
                       "fig2_get_first_violation_kepler_energy",
                       "fig2_get_minimum_violation_radius",
                       "fig2_get_maximum_violation_radius"):
            getattr(self.lib, getter).restype = ctypes.c_double
        self.function = self.lib.fig2_additional_forces

        if gravity is not None:
            radii = gravity.table_r
            fr = gravity.radial_acceleration(radii, 0)
            dfr = gravity.radial_derivative(radii, 0)
        else:
            radii = np.array([0.1, 20.0])
            fr = dfr = np.zeros(2)
        self.radii = np.ascontiguousarray(radii, dtype=np.float64)
        self.fr = np.ascontiguousarray(fr, dtype=np.float64)
        self.dfr = np.ascontiguousarray(dfr, dtype=np.float64)
        self.active = np.ones(n_particles, dtype=np.int32)
        self.event_code = np.zeros(n_particles, dtype=np.int32)
        self.event_time = np.full(n_particles, np.nan)
        self.event_radius = np.full(n_particles, np.nan)
        self.event_position = np.full((n_particles, 3), np.nan)
        self.event_velocity = np.full((n_particles, 3), np.nan)
        self.event_radial_velocity = np.full(n_particles, np.nan)
        self.event_kepler_energy = np.full(n_particles, np.nan)
        inner = boundaries.inner_termination_au if boundaries else 0.0
        outer = boundaries.outer_termination_au if boundaries else 0.0
        self.params = _Params(
            len(self.radii), self.radii.ctypes.data_as(DoublePtr),
            self.fr.ctypes.data_as(DoublePtr), self.dfr.ctypes.data_as(DoublePtr),
            tau_myr / MYR_PER_CODE,
            sigma0, alpha1, alpha2,
            h1, height_exponent,
            drag.flow.mu_code if drag else 1.0,
            drag.cd if drag else 0.0,
            drag.planetesimal.solid_density_cgs if drag else 1.0,
            drag.planetesimal.radius_cm if drag else 1.0,
            AU_CM, first_test, int(use_disk_gravity), int(use_drag),
            int(boundaries is not None), inner, outer,
            int(boundaries.require_positive_energy_for_ejection) if boundaries else 1,
            self.active.ctypes.data_as(IntPtr),
            self.event_code.ctypes.data_as(IntPtr),
            self.event_time.ctypes.data_as(DoublePtr),
            self.event_radius.ctypes.data_as(DoublePtr),
            self.event_position.ctypes.data_as(DoublePtr),
            self.event_velocity.ctypes.data_as(DoublePtr),
            self.event_radial_velocity.ctypes.data_as(DoublePtr),
            self.event_kepler_energy.ctypes.data_as(DoublePtr),
        )
        self.lib.fig2_set_params(ctypes.byref(self.params))

    def events(self, current_names, original_names):
        """Copy newly recorded C events before affected particles are removed."""
        found = []
        for index in np.flatnonzero(self.event_code):
            if index == 0 or index > len(current_names):
                continue
            name = current_names[index-1]
            code = int(self.event_code[index])
            found.append(ParticleEvent(
                particle_index=original_names.index(name), particle_name=name,
                event_code=code, event_type=EVENT_NAMES[code],
                time_myr=float(self.event_time[index] * MYR_PER_CODE),
                radius_au=float(self.event_radius[index]),
                position_au=tuple(float(x) for x in self.event_position[index]),
                velocity_code=tuple(float(x) for x in self.event_velocity[index]),
                radial_velocity_code=float(self.event_radial_velocity[index]),
                kepler_energy_code=float(self.event_kepler_energy[index]),
            ))
        return found

    def check_domain(self, names=()):
        count = self.lib.fig2_get_domain_violations()
        if count:
            index = self.lib.fig2_get_first_violation_particle()
            kind = self.lib.fig2_get_first_violation_kind()
            time_myr = self.lib.fig2_get_first_violation_time() * MYR_PER_CODE
            radius = self.lib.fig2_get_first_violation_radius()
            energy = self.lib.fig2_get_first_violation_kepler_energy()
            minimum = self.lib.fig2_get_minimum_violation_radius()
            maximum = self.lib.fig2_get_maximum_violation_radius()
            name = names[index-1] if 1 <= index <= len(names) else f"simulation index {index}"
            raise DiskForceDomainError(
                count=count, particle_index=index, particle_name=name, kind=kind,
                time_myr=time_myr, radius_au=radius,
                kepler_energy=energy,
                table_bounds_au=(self.radii[0], self.radii[-1]),
                offending_range_au=(minimum, maximum),
            )
