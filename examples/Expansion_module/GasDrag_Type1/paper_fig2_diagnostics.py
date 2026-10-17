"""Plotting and event helpers for Fig. 2-like orbital outputs."""

from __future__ import annotations

import numpy as np
from paper_fig2_integrator import Fig2Run
from paper_fig2_resonances import SecularResonanceTracks


def mmr_location_au(jupiter_a_au, numerator: int, denominator: int):
    """Interior p:q MMR radius for planet a_J, with p>q."""
    if not numerator > denominator > 0 or jupiter_a_au <= 0:
        raise ValueError("require a_J>0 and integer p>q>0")
    return jupiter_a_au * (denominator/numerator)**(2/3)


def default_mmr_locations(jupiter_a_au=5.2):
    return {"2:1": mmr_location_au(jupiter_a_au, 2, 1),
            "7:3": mmr_location_au(jupiter_a_au, 7, 3),
            "5:2": mmr_location_au(jupiter_a_au, 5, 2),
            "3:1": mmr_location_au(jupiter_a_au, 3, 1),
            "4:1": mmr_location_au(jupiter_a_au, 4, 1)}


def derive_orbital_elements(run: Fig2Run):
    """Return saved heliocentric osculating elements and peri/aphelion."""
    return {"t_myr": run.time_myr, "names": run.names,
            "a_au": run.a_au, "e": run.e,
            "omega_rad": run.omega_rad,
            "q_au": run.q_au, "Q_au": run.Q_au,
            "active": run.active, "events": run.events}


def plot_fig2_style(run: Fig2Run, resonances: SecularResonanceTracks | None = None,
                    tdep_myr: float = 1.0, ax=None):
    """Plot q/a/Q tracks, candidate secular roots and Jupiter MMRs.

    This figure does not imply quantitative reproduction of Zheng et al. Fig. 2.
    All roots at a given time are plotted; disconnected roots are not joined.
    """
    import matplotlib.pyplot as plt

    if tdep_myr <= 0:
        raise ValueError("tdep_myr must be positive")
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5.5))
    t = run.time_myr/tdep_myr
    good = t > 0
    colors = ("tab:red", "limegreen", "tab:blue", "magenta")
    test_names = [name for name in run.names if name.startswith("particle_")]
    for j, name in enumerate(test_names):
        idx = run.names.index(name)
        color = colors[j % len(colors)]
        ax.plot(t[good], run.a_au[good, idx], color=color, lw=1.3,
                label=f"{name}: a")
        ax.plot(t[good], run.q_au[good, idx], color=color, lw=0.6,
                alpha=0.65)
        ax.plot(t[good], run.Q_au[good, idx], color=color, lw=0.6,
                alpha=0.65)
    jupiter_a = (run.a_au[:, run.names.index("Jupiter")]
                 if "Jupiter" in run.names else np.full(len(t), 5.2))
    for label, loc in default_mmr_locations(float(jupiter_a[0])).items():
        ax.axhline(loc, color="gold", lw=0.8, alpha=0.8)
        ax.text(t[good][-1], loc, label, ha="right", va="bottom", fontsize=7)
    if resonances is not None:
        for label, roots_by_time, color in (
            ("nu5 candidate", resonances.nu5_roots_au, "black"),
            ("nu6 candidate", resonances.nu6_roots_au, "gray"),
        ):
            xs, ys = [], []
            for ti, roots in zip(resonances.time_myr/tdep_myr, roots_by_time):
                if ti <= 0:
                    continue
                for root in roots:
                    xs.append(ti)
                    ys.append(root)
            if xs:
                ax.scatter(xs, ys, c=color, s=4, label=label, zorder=4)
    ax.set(xscale="log", xlabel=r"$t/T_{dep}$", ylabel="q, a, Q (AU)",
           title="Fig. 2 dynamical experiment (model assumptions in metadata)")
    ymin, ymax = ax.get_ylim()
    for event in run.events:
        if not event.particle_name.startswith("particle_") or event.time_myr <= 0:
            continue
        color = colors[test_names.index(event.particle_name) % len(colors)]
        event_x = event.time_myr/tdep_myr
        if ymin <= event.radius_au <= ymax:
            ax.scatter(event_x, event.radius_au, marker="x", s=34,
                       linewidths=1.2, color=color, zorder=6)
        else:
            ax.plot(event_x, 0.97, marker="x", ms=6, color=color,
                    transform=ax.get_xaxis_transform(), zorder=6)
            ax.annotate(f"{event.particle_name}: {event.event_type}, "
                        f"r={event.radius_au:.2f} AU",
                        xy=(event_x, 0.97), xycoords=ax.get_xaxis_transform(),
                        xytext=(3, -4), textcoords="offset points",
                        fontsize=6, color=color, rotation=90, va="top")
    ax.set_ylim(ymin, ymax)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=7, ncol=2)
    return ax
