# Gas drag, Type-I migration, and secular resonance experiments

This directory contains research examples for comparing aerodynamic drag with
Type-I migration, integrating test-particle orbits in a dissipating gas disk,
and tracking secular resonances. These are experimental workflows; passing the
tests does not establish quantitative reproduction of every reference figure.

## Setup

From the repository root, using Python 3.10 or later:

```bash
python -m pip install -e ".[dev,examples]"
python -m pip install -r examples/Expansion_module/GasDrag_Type1/requirements.txt
```

The additional dependency is REBOUND 4.4.8, the version used for validation.
NumPy, SciPy, plotting, notebook, and test dependencies come from AccreDiff and
its extras. This installation step is not performed automatically by notebooks.

The `python` force backend requires no local C compilation. The `compiled`
backend, and tests that exercise it, require a C compiler exposed as `cc` and
the `rebound.h` header supplied by REBOUND. On macOS, install the Xcode Command
Line Tools if `cc` is unavailable. The bridge builds its library in the system
temporary directory under `accrediff_fig2_compiled`; generated libraries are
not committed. Validation was performed on macOS; other platforms have not
been validated here. Compiled force parameters are process-global: use
separate processes for concurrent experiments with different configurations.

## Notebook entry points

| Notebook | Purpose |
| --- | --- |
| `run_gas_drag_type1.ipynb` | Compare drag and Type-I migration/damping time scales |
| `gas_drag_dispersal_tracks.ipynb` | Compare independent drag and Type-I tracks during gas dispersal |
| `reproduce_thommes2008_e_damping.ipynb` | Explore the reference eccentricity-damping calculation |
| `run_paper_fig2_drag.ipynb` | Configure direct orbital integration with drag and disk forces |
| `run_paper_fig2_resonances.ipynb` | Combine orbital evolution and resonance-track diagnostics |

Launch notebooks from the repository root or this directory. Review the time
span, time step, force backend, and execution switches before running: long
direct integrations can be expensive. Saved notebook outputs may originate
from earlier runs; regenerate results for the selected configuration.

## Modules and units

- `src/accrediff/gas_drag_type1_model.py` provides the reusable comparison
  configuration and factory, built on `accrediff.gas_migration`.
- `paper_fig2_disk.py`, `paper_fig2_drag.py`, and `paper_fig2_orbit_average.py`
  implement the example disk adapters, drag conventions, and orbit averaging.
- `paper_fig2_gravity.py` supplies tabulated disk gravity;
  `paper_fig2_resonances.py` computes secular resonance tracks.
- `paper_fig2_system.py`, `paper_fig2_integrator.py`, and
  `paper_fig2_events.py` configure initial conditions, integration, and
  per-particle boundary events.
- `paper_fig2_compiled.py` and `paper_fig2_force.c` implement the optional
  compiled force callback. `paper_fig2_experiment.py` checks saved-run
  configuration compatibility; `paper_fig2_diagnostics.py` plots diagnostics.

Orbital integration uses AU, solar masses, and G=1; its internal time unit is
year/(2*pi). Saved evolution times use Myr. Disk and drag calculations also
use explicit CGS quantities at their interfaces; use the provided adapters
rather than passing Myr directly to code-time functions.

Generated NPZ tracks and plots belong in `outputs/`, which is ignored by Git.
The directory's development-plan documents record research decisions and
pending work; they are not a guarantee that every proposed feature is complete.

## Verification

Run from the repository root:

```bash
python -m pytest -q examples/Expansion_module/GasDrag_Type1/tests
```

The default project test command only collects `tests/`; explicitly include
this directory when validating these experiments. Tests cover disk adapters,
drag and orbit averaging, disk gravity, short integrations, Python/compiled
force agreement, boundary events, saved-run compatibility, and resonances.
They do not run the full 10 Myr notebook experiments.
