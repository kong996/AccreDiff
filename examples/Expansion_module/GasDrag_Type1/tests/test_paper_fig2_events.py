import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_disk import PaperFig2Disk
from paper_fig2_events import (
    OUTER_DOMAIN_EXIT_BOUND,
    ParticleBoundaryConfig,
    classify_state,
)
from paper_fig2_integrator import Fig2Run, integrate_fig2
from paper_fig2_system import Fig2SystemConfig, IntegrationConfig


def _system(semimajor_axes):
    count = len(semimajor_axes)
    return Fig2SystemConfig(
        use_planets=False, use_disk_gravity=False, use_drag=False,
        particle_a_au=tuple(semimajor_axes), particle_e=0.5,
        particle_omega_rad=(0.0,)*count,
        particle_mean_anomaly_rad=(0.0,)*count)


def test_boundary_configuration_and_classifier():
    config = ParticleBoundaryConfig(0.2, 3.0, 10.0)
    config.validate_force_table(0.1, 10.0)
    with pytest.raises(ValueError, match="strictly inside"):
        config.validate_force_table(0.2, 10.0)
    result = classify_state([3.1, 0, 0], [0, 0.5, 0], 1.0, config)
    assert result[0] == OUTER_DOMAIN_EXIT_BOUND
    assert result[1] == 3.1


@pytest.mark.parametrize("backend", ["python", "compiled"])
def test_escaped_test_particle_does_not_change_survivor(backend, tmp_path):
    disk = PaperFig2Disk()
    boundaries = ParticleBoundaryConfig(0.1, 3.0, 0.25)
    integration = IntegrationConfig("whfast", 0.01, backend)
    times = np.linspace(0, 2e-5, 21)  # 20 years
    combined = integrate_fig2(_system((1.0, 2.3)), disk, None, None,
                              times, integration, boundaries)
    control = integrate_fig2(_system((1.0,)), disk, None, None,
                             times, integration, boundaries)

    assert len(combined.events) == 1
    event = combined.events[0]
    assert event.particle_name == "particle_1"
    assert event.event_type == "outer_domain_exit_bound"
    assert event.radius_au >= boundaries.outer_termination_au
    assert not combined.active[-1, 1]
    assert np.all(np.isnan(combined.a_au[~combined.active[:, 1], 1]))
    np.testing.assert_allclose(combined.position_au[:, 0],
                               control.position_au[:, 0], atol=2e-12)
    np.testing.assert_allclose(combined.a_au[:, 0], control.a_au[:, 0], atol=2e-12)
    np.testing.assert_allclose(combined.e[:, 0], control.e[:, 0], atol=2e-12)

    destination = tmp_path / f"events_{backend}.npz"
    combined.save_npz(destination)
    restored = Fig2Run.load_npz(destination)
    np.testing.assert_array_equal(restored.active, combined.active)
    assert restored.events == combined.events
    np.testing.assert_allclose(restored.a_au, combined.a_au, equal_nan=True)
