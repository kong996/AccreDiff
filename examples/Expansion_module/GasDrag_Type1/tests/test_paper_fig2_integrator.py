import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_disk import GasMigrationDiskAdapter, PaperFig2Disk
from paper_fig2_drag import AerodynamicDrag, GasFlowPolicy, Planetesimal
from paper_fig2_gravity import GapGravityConfig, PaperGapDiskGravity
from paper_fig2_integrator import Fig2Run, integrate_fig2
from paper_fig2_system import Fig2SystemConfig, IntegrationConfig

from accrediff.gas_migration import Gas_DiskModel, Gas_ModelConfig, Gas_UnitSystem


def test_pure_kepler_track_and_force_switches():
    disk = PaperFig2Disk()
    base = Fig2SystemConfig(use_planets=False, use_disk_gravity=False,
                            use_drag=False, particle_a_au=(2.3,),
                            particle_omega_rad=(0.0,),
                            particle_mean_anomaly_rad=(0.0,))
    times = np.array([0.0, 2e-6])  # 2 years
    integ = IntegrationConfig("ias15", 0.1)
    kepler = integrate_fig2(base, disk, None, None, times, integ)
    assert kepler.names == ("particle_0",)
    np.testing.assert_allclose(kepler.a_au[:, 0], 2.3, atol=1e-10)
    np.testing.assert_allclose(kepler.e[:, 0], base.particle_e, atol=1e-10)

    drag = AerodynamicDrag(Planetesimal(), GasFlowPolicy())
    # An instantaneous osculating a can increase during one orbital phase;
    # compare over many orbits to test the secular inward drift.
    drag_run = integrate_fig2(replace(base, use_drag=True), disk, None, drag,
                              np.array([0.0, 1e-4]),
                              IntegrationConfig("whfast", 0.05))
    assert drag_run.a_au[-1, 0] < kepler.a_au[-1, 0]


def test_individual_initial_semimajor_axes_and_eccentricities():
    disk = PaperFig2Disk()
    system = Fig2SystemConfig(
        use_planets=False, use_disk_gravity=False, use_drag=False,
        particle_a_au=(0.8, 1.4, 2.6),
        particle_e=(0.01, 0.08, 0.2),
        particle_omega_rad=(0.0, 0.1, 0.2),
        particle_mean_anomaly_rad=(0.0, 0.0, 0.0),
    )
    run = integrate_fig2(system, disk, None, None, [0.0],
                         IntegrationConfig("ias15", 0.1))
    np.testing.assert_allclose(run.a_au[0], system.particle_a_au, rtol=1e-12)
    np.testing.assert_allclose(run.e[0], system.particle_e_values, rtol=1e-12)
    with pytest.raises(ValueError, match="length N"):
        Fig2SystemConfig(particle_e=(0.01, 0.02))


def test_planet_perturbations_and_saved_shapes(tmp_path):
    disk = PaperFig2Disk()
    system = Fig2SystemConfig()
    gravity = PaperGapDiskGravity(disk)
    drag = AerodynamicDrag(Planetesimal(), GasFlowPolicy())
    run = integrate_fig2(system, disk, gravity, drag, np.array([0, 1e-6]),
                         IntegrationConfig("ias15", 0.1))
    assert run.position_au.shape == (2, 6, 3)
    assert run.q_au.shape == run.a_au.shape
    assert "Jupiter" in run.names and "particle_0" in run.names
    np.testing.assert_allclose(run.a_au[0],
        [system.jupiter_a_au, system.saturn_a_au, *system.particle_a_au],
        rtol=1e-12)
    # Zero-mass test particles do not alter giant-planet tracks.
    one_particle = replace(system, particle_a_au=(3.5,),
                           particle_omega_rad=(0.0,),
                           particle_mean_anomaly_rad=(0.0,))
    run_one = integrate_fig2(one_particle, disk, gravity, drag,
                             np.array([0, 1e-6]), IntegrationConfig("ias15", 0.1))
    np.testing.assert_allclose(run.position_au[:, :2],
                               run_one.position_au[:, :2], atol=1e-10)
    dest = tmp_path / "short_run.npz"
    run.save_npz(dest)
    with np.load(dest) as saved:
        np.testing.assert_allclose(saved["a_au"], run.a_au)
    restored = Fig2Run.load_npz(dest)
    assert restored.names == run.names
    np.testing.assert_allclose(restored.q_au, run.q_au)


def test_compiled_force_matches_python_over_many_orbits():
    disk = PaperFig2Disk()
    system = Fig2SystemConfig()
    gravity = PaperGapDiskGravity(disk)
    drag = AerodynamicDrag(Planetesimal(), GasFlowPolicy())
    times = np.array([0.0, 1e-4])  # 100 years
    python_run = integrate_fig2(system, disk, gravity, drag, times,
                                IntegrationConfig("whfast", 0.025, "python"))
    compiled_run = integrate_fig2(system, disk, gravity, drag, times,
                                  IntegrationConfig("whfast", 0.025, "compiled"))
    np.testing.assert_allclose(compiled_run.a_au, python_run.a_au, atol=1e-10)
    np.testing.assert_allclose(compiled_run.e, python_run.e, atol=1e-10)


def test_gas_migration_disk_compiled_matches_python():
    units = Gas_UnitSystem()
    config = Gas_ModelConfig(Sigma0_cgs=1700, a_alpha1=1.0, a_alpha2=0.5,
                             tau_decay_myr=1.0)
    disk = GasMigrationDiskAdapter(Gas_DiskModel(config, units), units)
    gravity = PaperGapDiskGravity(disk)
    drag = AerodynamicDrag(Planetesimal(), GasFlowPolicy())
    system = Fig2SystemConfig()
    times = np.array([0.0, 1e-4])
    python_run = integrate_fig2(system, disk, gravity, drag, times,
                                IntegrationConfig("whfast", 0.025, "python"))
    compiled_run = integrate_fig2(system, disk, gravity, drag, times,
                                  IntegrationConfig("whfast", 0.025, "compiled"))
    np.testing.assert_allclose(compiled_run.a_au, python_run.a_au, atol=1e-10)
    np.testing.assert_allclose(compiled_run.e, python_run.e, atol=1e-10)


def test_matched_adapter_reproduces_paper_short_track():
    units = Gas_UnitSystem()
    config = Gas_ModelConfig(Sigma0_cgs=1700, a_alpha1=0, a_alpha2=1.5,
                             def_h_1=0.025, beta=0.25, tau_decay_myr=1)
    adapted = GasMigrationDiskAdapter(Gas_DiskModel(config, units), units)
    paper = PaperFig2Disk()
    drag = AerodynamicDrag(Planetesimal(), GasFlowPolicy())
    system = Fig2SystemConfig()
    times = np.array([0.0, 1e-4])
    cfg = IntegrationConfig("whfast", 0.025, "compiled")
    baseline = integrate_fig2(system, paper, PaperGapDiskGravity(paper), drag, times, cfg)
    matched = integrate_fig2(system, adapted, PaperGapDiskGravity(adapted), drag, times, cfg)
    np.testing.assert_allclose(matched.a_au, baseline.a_au, atol=1e-10)
    np.testing.assert_allclose(matched.e, baseline.e, atol=1e-10)


def test_compiled_domain_error_identifies_body_time_and_radius():
    disk = PaperFig2Disk()
    gravity = PaperGapDiskGravity(
        disk, GapGravityConfig(table_min_au=2.0, table_max_au=2.29))
    system = Fig2SystemConfig(
        use_planets=False, use_drag=False, particle_a_au=(2.3,),
        particle_omega_rad=(0.0,), particle_mean_anomaly_rad=(0.0,))
    with pytest.raises(RuntimeError, match=(
            r"first radial table boundary violation by particle_0 .*"
                r"valid table=\[2, 2\.29\] AU")):
        integrate_fig2(system, disk, gravity, None, [0.0, 5e-6],
                       IntegrationConfig("whfast", 0.025, "compiled"))
