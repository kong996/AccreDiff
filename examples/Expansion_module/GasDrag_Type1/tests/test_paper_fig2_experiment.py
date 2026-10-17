import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_drag import AerodynamicDrag, GasFlowPolicy, Planetesimal
from paper_fig2_experiment import make_disk, matches_configuration, output_path
from paper_fig2_gravity import PaperGapDiskGravity
from paper_fig2_integrator import integrate_fig2
from paper_fig2_system import Fig2SystemConfig, IntegrationConfig


def test_named_disks_match_or_depart_as_documented():
    paper = make_disk("paper_reference")
    matched = make_disk("gas_migration_matched")
    target = make_disk("gas_profile_baseline")
    radii = np.array([0.5, 1.0, 2.3, 3.5, 5.2, 9.58, 12.0])
    for t in (0.0, 1.0, 5.0):
        np.testing.assert_allclose(matched.sigma_cgs(radii, t),
                                   paper.sigma_cgs(radii, t), rtol=1e-13)
        np.testing.assert_allclose(matched.scale_height_au(radii, t),
                                   paper.scale_height_au(radii, t), rtol=1e-13)
    np.testing.assert_allclose(matched.sigma_slope_p(radii), 1.5)
    assert not np.isclose(target.sigma_cgs(3.5, 0), paper.sigma_cgs(3.5, 0))
    assert target.disk.cfg.a_alpha1 == 1.0
    assert target.disk.cfg.a_alpha2 == 0.5
    assert target.disk.cfg.tau_decay_myr == 1.0


def test_saved_run_match_rejects_changed_initial_conditions_and_radius(tmp_path):
    disk = make_disk("paper_reference")
    system = replace(Fig2SystemConfig(), use_planets=False, use_disk_gravity=False,
                     particle_a_au=(2.3,), particle_omega_rad=(0.0,),
                     particle_mean_anomaly_rad=(0.0,))
    drag = AerodynamicDrag(Planetesimal(radius_km=10), GasFlowPolicy())
    cfg = IntegrationConfig("whfast", 0.05, "compiled")
    run = integrate_fig2(system, disk, None, drag, [0.0, 1e-6], cfg)
    assert matches_configuration(run, system, disk, None, drag, cfg)
    big_body = AerodynamicDrag(Planetesimal(radius_km=100), GasFlowPolicy())
    assert not matches_configuration(run, system, disk, None, big_body, cfg)
    assert not matches_configuration(run, replace(system, particle_e=0.02),
                                     disk, None, drag, cfg)
    assert output_path(tmp_path, "gas_profile_baseline") != output_path(
        tmp_path, "paper_reference")


def test_target_disk_generates_finite_gap_force():
    disk = make_disk("gas_profile_baseline")
    gravity = PaperGapDiskGravity(disk)
    radii = np.array([1.0, 2.3, 3.5, 5.2, 9.58, 12.0])
    assert np.all(np.isfinite(gravity.radial_acceleration(radii, 0)))
    assert np.all(np.isfinite(gravity.radial_derivative(radii, 0)))
    np.testing.assert_allclose(gravity.radial_acceleration(radii, 1),
                               gravity.radial_acceleration(radii, 0) * np.exp(-1),
                               rtol=1e-12)
