import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_disk import GasMigrationDiskAdapter, PaperFig2Disk

from accrediff.gas_migration import Gas_build_model, Gas_ModelConfig


def test_paper_normalization_decay_and_density():
    disk = PaperFig2Disk()
    assert disk.sigma_cgs(1, 0) == 1700
    assert disk.scale_height_au(1, 0) == 0.025
    assert np.isclose(disk.sigma_cgs(1, 1) / disk.sigma_cgs(1, 0), np.exp(-1))
    assert np.isclose(disk.sigma_cgs(3, 0), 1700 * 3**-1.5)
    assert np.isclose(disk.sigma_slope_p(3), 1.5)


def test_existing_disk_adapter_reproduces_paper_profile():
    model = Gas_build_model(Gas_ModelConfig(
        Sigma0_cgs=1700, a_alpha1=0, a_alpha2=1.5,
        def_h_1=0.025, beta=0.25, tau_decay_myr=1,
    ))
    paper = PaperFig2Disk()
    adapted = GasMigrationDiskAdapter(model["disk"], model["units"])
    radius = np.array([1.0, 2.3, 3.5])
    for t in (0.0, 1.0, 5.0):
        np.testing.assert_allclose(adapted.sigma_cgs(radius, t),
                                   paper.sigma_cgs(radius, t), rtol=1e-13)
        np.testing.assert_allclose(adapted.scale_height_au(radius, t),
                                   paper.scale_height_au(radius, t), rtol=1e-13)
    np.testing.assert_allclose(adapted.sigma_slope_p(radius),
                               paper.sigma_slope_p(radius))
