import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_diagnostics import mmr_location_au
from paper_fig2_disk import PaperFig2Disk
from paper_fig2_gravity import PaperGapDiskGravity
from paper_fig2_resonances import (
    laplace_coefficient,
    locate_secular_resonances,
    secular_matrix,
)
from paper_fig2_system import Fig2SystemConfig


def test_laplace_and_mmr_baselines():
    assert np.isclose(laplace_coefficient(0.0, 0), 2, rtol=1e-12)
    assert np.isclose(laplace_coefficient(0.0, 1), 0, atol=1e-12)
    assert np.isclose(mmr_location_au(5.2, 2, 1),
                      5.2 * 2**(-2/3), rtol=1e-12)


def test_disk_changes_secular_matrix_and_roots_sweep_inward():
    system = Fig2SystemConfig()
    gravity = PaperGapDiskGravity(PaperFig2Disk())
    assert not np.allclose(secular_matrix(system, gravity, 0),
                           secular_matrix(system, None, 0))
    t = np.array([1.0, 2.0])
    result = locate_secular_resonances(
        system, gravity, t, np.linspace(1.5, 4.3, 300), n_nodes=96)
    assert len(result.nu5_roots_au[0]) == 1
    assert len(result.nu5_roots_au[1]) == 1
    assert result.nu5_roots_au[1][0] < result.nu5_roots_au[0][0]
    assert result.g5_code.shape == t.shape
