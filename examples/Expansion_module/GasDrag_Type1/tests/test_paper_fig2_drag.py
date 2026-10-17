import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from paper_fig2_disk import PaperFig2Disk
from paper_fig2_drag import (
    AU_CM,
    AerodynamicDrag,
    GasFlowPolicy,
    Planetesimal,
)
from paper_fig2_orbit_average import orbit_averaged_rates


def test_planetesimal_mass_and_midplane_conventions():
    disk = PaperFig2Disk()
    body = Planetesimal()
    assert np.isclose(body.mass_g, 4*np.pi/3 * (1e6)**3)
    assert np.isclose(body.mass_mearth, 7.013640807294941e-10)
    paper = GasFlowPolicy("paper_2016")
    genga = GasFlowPolicy("genga_style")
    assert np.isclose(paper.density_cgs(disk, 1, 0), 1700/(0.025*AU_CM))
    assert np.isclose(paper.density_cgs(disk, 1, 0) /
                      genga.density_cgs(disk, 1, 0), np.sqrt(2*np.pi))
    assert np.isclose(paper.eta(disk, 1, 0), (0.025/4)**2)


def test_drag_coefficient_energy_and_vectorization():
    disk = PaperFig2Disk()
    body = Planetesimal()
    flow = GasFlowPolicy("paper_2016")
    paper = AerodynamicDrag(body, flow, law="paper_2016", drag_coefficient=0.165)
    genga = AerodynamicDrag(body, flow, law="genga_style", drag_coefficient=0.165)
    x = np.array([[1.0, 0.0], [0.0, 2.3]])
    gas_v = flow.velocity_code(x, 0, disk)
    assert np.all(paper.acceleration(x, gas_v, 0, disk) == 0)
    v = gas_v + np.array([[0.01, 0.02], [-0.03, 0.01]])
    a_paper = paper.acceleration(x, v, 0, disk)
    a_genga = genga.acceleration(x, v, 0, disk)
    np.testing.assert_allclose(a_genga, 3/8 * a_paper, rtol=1e-14)
    assert np.all(np.sum((v-gas_v)*a_paper, axis=-1) < 0)


def test_orbit_average_converges_and_damps_eccentricity():
    disk = PaperFig2Disk()
    drag = AerodynamicDrag(Planetesimal(), GasFlowPolicy())
    coarse = orbit_averaged_rates(3.0, 0.01, 0.0, disk, drag, n_samples=64)
    fine = orbit_averaged_rates(3.0, 0.01, 0.0, disk, drag, n_samples=256)
    assert coarse.da_dt_au_per_myr < 0
    assert coarse.de_dt_per_myr < 0
    assert np.isclose(coarse.da_dt_au_per_myr, fine.da_dt_au_per_myr, rtol=1e-5)
    assert np.isclose(coarse.de_dt_per_myr, fine.de_dt_per_myr, rtol=1e-5)
