"""Linear secular diagnostics for the Fig. 2 Sun-Jupiter-Saturn model.

These curves identify candidate nu5/nu6 locations. They are diagnostics,
not forces, and linear secular theory is unreliable near strong MMRs.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations

import numpy as np
from paper_fig2_gravity import PaperGapDiskGravity
from paper_fig2_system import Fig2SystemConfig


def laplace_coefficient(alpha, harmonic: int, n_nodes: int = 256):
    """b_{3/2}^{(j)}(alpha), Gauss-Legendre quadrature over [0,2pi]."""
    a = np.asarray(alpha, dtype=float)
    if np.any(a < 0) or np.any(a >= 1) or not np.all(np.isfinite(a)):
        raise ValueError("Laplace coefficient requires 0 <= alpha < 1")
    if harmonic < 0 or n_nodes < 32:
        raise ValueError("invalid harmonic or quadrature resolution")
    x, w = np.polynomial.legendre.leggauss(n_nodes)
    theta = np.pi * (x + 1)
    denom = 1 - 2*a[..., None]*np.cos(theta) + a[..., None]**2
    value = np.sum(w * np.cos(harmonic*theta) / denom**1.5,
                   axis=-1)
    return float(value) if np.ndim(alpha) == 0 else value


def disk_apsidal_rate(r_au, t_myr: float, gravity: PaperGapDiskGravity,
                      star_mass_msun: float = 1.0):
    """Omega-kappa for nearly circular planar orbits, in rad/code time."""
    r = np.asarray(r_au, dtype=float)
    fr = gravity.radial_acceleration(r, t_myr)
    dfr = gravity.radial_derivative(r, t_myr)
    omega2 = star_mass_msun/r**3 - fr/r
    kappa2 = star_mass_msun/r**3 - dfr - 3*fr/r
    if np.any(omega2 <= 0) or np.any(kappa2 <= 0):
        raise ValueError("disk model has unstable or invalid circular frequency")
    out = np.sqrt(omega2) - np.sqrt(kappa2)
    return float(out) if out.ndim == 0 else out


def secular_matrix(system: Fig2SystemConfig, gravity: PaperGapDiskGravity | None,
                   t_myr: float, n_nodes: int = 256):
    """Planar Jupiter-Saturn Laplace-Lagrange matrix in rad/code time."""
    aj, ast = system.jupiter_a_au, system.saturn_a_au
    if not aj < ast:
        raise ValueError("Jupiter must be interior to Saturn")
    alpha = aj/ast
    b1 = laplace_coefficient(alpha, 1, n_nodes)
    b2 = laplace_coefficient(alpha, 2, n_nodes)
    nj = np.sqrt(system.star_mass_msun/aj**3)
    ns = np.sqrt(system.star_mass_msun/ast**3)
    mj, ms, mstar = (system.jupiter_mass_msun,
                      system.saturn_mass_msun, system.star_mass_msun)
    # The inner body's secular coefficient has alpha*baralpha=alpha**2;
    # the outer body's has alpha*baralpha=alpha.
    cj = nj*ms/(4*(mstar+mj))*alpha**2
    cs = ns*mj/(4*(mstar+ms))*alpha
    matrix = np.array([[cj*b1, -cj*b2], [-cs*b2, cs*b1]], dtype=float)
    if gravity is not None:
        matrix[0, 0] += disk_apsidal_rate(aj, t_myr, gravity, mstar)
        matrix[1, 1] += disk_apsidal_rate(ast, t_myr, gravity, mstar)
    return matrix


def particle_apsidal_rate(a_au, system: Fig2SystemConfig,
                          gravity: PaperGapDiskGravity | None,
                          t_myr: float, n_nodes: int = 256):
    """Free test-particle apsidal rate inside Jupiter, in rad/code time."""
    a = np.asarray(a_au, dtype=float)
    if np.any(a <= 0) or np.any(a >= system.jupiter_a_au):
        raise ValueError("particle radii must be positive and inside Jupiter")
    n = np.sqrt(system.star_mass_msun/a**3)
    rate = np.zeros_like(a)
    for mass, ap in ((system.jupiter_mass_msun, system.jupiter_a_au),
                     (system.saturn_mass_msun, system.saturn_a_au)):
        alpha = a/ap
        rate += n*mass/(4*system.star_mass_msun)*alpha**2 * laplace_coefficient(
            alpha, 1, n_nodes)
    if gravity is not None:
        rate += disk_apsidal_rate(a, t_myr, gravity, system.star_mass_msun)
    return float(rate) if rate.ndim == 0 else rate


@dataclass
class SecularResonanceTracks:
    time_myr: np.ndarray
    g5_code: np.ndarray
    g6_code: np.ndarray
    nu5_roots_au: tuple[np.ndarray, ...]
    nu6_roots_au: tuple[np.ndarray, ...]


def _mode_frequencies(matrix, previous_vectors=None):
    values, vectors = np.linalg.eig(matrix)
    if np.max(np.abs(values.imag)) > 1e-12:
        raise ValueError("complex secular eigenfrequencies")
    values, vectors = values.real, vectors.real
    vectors /= np.linalg.norm(vectors, axis=0)
    if previous_vectors is None:
        # Initially label the mode with greater Jupiter participation as g5.
        order = (0, 1) if abs(vectors[0, 0]) >= abs(vectors[0, 1]) else (1, 0)
    else:
        candidates = list(permutations((0, 1)))
        scores = [sum(abs(np.dot(previous_vectors[:, j], vectors[:, order[j]]))
                      for j in (0, 1)) for order in candidates]
        order = candidates[int(np.argmax(scores))]
    return values[list(order)], vectors[:, list(order)]


def locate_secular_resonances(system: Fig2SystemConfig,
                              gravity: PaperGapDiskGravity | None,
                              time_myr, radius_grid_au,
                              n_nodes: int = 256) -> SecularResonanceTracks:
    """Find all roots g_particle(a,t)=g5(t) or g6(t) on a radius grid."""
    times = np.asarray(time_myr, dtype=float)
    radii = np.asarray(radius_grid_au, dtype=float)
    if (times.ndim != 1 or radii.ndim != 1 or len(times) == 0 or
        len(radii) < 3 or np.any(np.diff(times) <= 0) or
        np.any(np.diff(radii) <= 0) or np.any(times < 0) or
        radii[0] <= 0 or radii[-1] >= system.jupiter_a_au):
        raise ValueError("times and radii must be ordered, valid 1D grids")
    g5, g6 = np.empty(len(times)), np.empty(len(times))
    roots5, roots6 = [], []
    previous = None
    for i, t in enumerate(times):
        modes, previous = _mode_frequencies(
            secular_matrix(system, gravity, float(t), n_nodes), previous)
        g5[i], g6[i] = modes
        gp = particle_apsidal_rate(radii, system, gravity, float(t), n_nodes)
        for mode, output in ((g5[i], roots5), (g6[i], roots6)):
            diff = gp-mode
            roots = []
            for j in np.where(diff[:-1]*diff[1:] < 0)[0]:
                # Linear interpolation between calculated secular rates.
                root = radii[j] - diff[j] * (radii[j+1]-radii[j]) / (
                    diff[j+1]-diff[j])
                roots.append(root)
            roots.extend(float(radii[j]) for j in np.where(diff == 0)[0])
            output.append(np.asarray(sorted(set(roots)), dtype=float))
    return SecularResonanceTracks(times, g5, g6,
                                  tuple(roots5), tuple(roots6))
