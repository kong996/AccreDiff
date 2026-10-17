import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_disk import PaperFig2Disk
from paper_fig2_gravity import GapGravityConfig, PaperGapDiskGravity


def test_force_is_potential_gradient_and_depletes():
    gravity = PaperGapDiskGravity(PaperFig2Disk())
    for r in (2.3, 3.5, 5.2, 9.58):
        h = 1e-4
        finite_difference = -(gravity.potential(r+h) - gravity.potential(r-h))/(2*h)
        assert np.isclose(gravity.radial_acceleration(r), finite_difference,
                          rtol=2e-4)
        vec = gravity.acceleration(np.array([r, 0.0]), 0)
        np.testing.assert_allclose(vec, [gravity.radial_acceleration(r), 0])
        assert np.isclose(gravity.radial_acceleration(r, 1) /
                          gravity.radial_acceleration(r, 0), np.exp(-1), rtol=1e-12)


def test_ring_and_table_resolution_stability_away_from_gap_edges():
    disk = PaperFig2Disk()
    coarse = PaperGapDiskGravity(disk, GapGravityConfig(
        n_ring_per_region=160, n_table=384))
    fine = PaperGapDiskGravity(disk, GapGravityConfig(
        n_ring_per_region=320, n_table=768))
    for r in (2.3, 3.0, 5.2, 9.58):
        assert np.isclose(coarse.radial_acceleration(r),
                          fine.radial_acceleration(r), rtol=0.02)
