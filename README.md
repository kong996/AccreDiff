# 🪐 AccreDiff

**AccreDiff** is an open-source Python framework for coupling  
**N-body accretion histories** with **core–mantle differentiation models**.

It provides a modular, reproducible, and extensible toolkit for studying
the chemical evolution of terrestrial planets.

AccreDiff is designed for community use and collaborative development.

---

## Overview

AccreDiff enables researchers to:

- Load N-body simulation outputs (e.g., GENGA)
- Extract and process impact histories
- Apply metal–silicate equilibration models
- Track core–mantle compositional evolution
- Control redox evolution (ΔIW-based models)
- Analyze and compare planetary analogs
- Generate publication-ready visualizations

The framework separates physical models from workflow logic, allowing
new differentiation models or partitioning laws to be integrated easily.

---

## Design Principles

AccreDiff is built around:

- **Modularity** — Physical models are interchangeable.
- **Reproducibility** — All parameters are controlled via configuration files.
- **Transparency** — No hidden constants; all assumptions are explicit.
- **Extensibility** — Designed for future community contributions.

---

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/kong996/AccreDiff.git
cd AccreDiff
pip install -r requirements.txt
```

Or install in editable mode for development:

```bash
pip install -e .
```

### Requirements

- Python ≥ 3.8
- numpy ≥ 1.21.0
- scipy ≥ 1.7.0
- pandas ≥ 1.3.0
- matplotlib ≥ 3.4.0
- pyyaml ≥ 6.0
- jupyter ≥ 1.0.0
- ipykernel ≥ 6.0.0
- ipywidgets ≥ 7.6.0
- tqdm ≥ 4.62.0

---

## Project Structure

```
AccreDiff/
├── src/
│   └── accrediff/          # Core Python package
├── examples/               # Example workflows & tutorial notebooks
│   ├── quickstart/         # Step-by-step tutorial notebooks
│   ├── basics/             # Basic usage notebooks
│   └── dataset/            # Example datasets
├── notebooks/              # Additional Jupyter Notebooks
├── tests/                  # Unit tests
├── requirements.txt
├── pyproject.toml
└── README.md
```

---

## Quick Start

The recommended way to get started is through the step-by-step quickstart notebooks.  
Each notebook builds on the previous one and covers a complete accretion–differentiation workflow:

| Notebook | Description |
|----------|-------------|
| [`01_Constants_and_setup.ipynb`](examples/quickstart/01_Constants_and_setup.ipynb) | Load built-in constants, element tables, and meteorite / planet reference compositions |
| [`02_Meteorite_bulk_composition.ipynb`](examples/quickstart/02_Meteorite_bulk_composition.ipynb) | Convert wt% → molar, normalize to MgO & CI, export to CSV |
| [`03_Embryo_bulk_composition.ipynb`](examples/quickstart/03_Embryo_bulk_composition.ipynb) | Reconstruct embryo bulk compositions at ~5 Myr from N-body collision histories |
| [`04_Accrection_history.ipynb`](examples/quickstart/04_Accrection_history.ipynb) | Track mass evolution and reconstruct merger trees across 0–100 Myr |
| [`05_Impact_events.ipynb`](examples/quickstart/05_Impact_events.ipynb) | Classify impact events (small / global / partial melting) and assign equilibration pressure |
| [`06_Differentiation.ipynb`](examples/quickstart/06_Differentiation.ipynb) | Simulate multi-stage metal–silicate differentiation using K_D equilibrium solver |
| [`07_Comparison_bulk.ipynb`](examples/quickstart/07_Comparison_bulk.ipynb) | Compare modeled planet compositions against Earth & Mars reference values |

A minimal usage example:

```python
import accrediff as ad

# Load built-in constants and reference compositions
elements  = ad.constants.Elements
earth     = ad.constants.Earth
ci        = ad.constants.CI_bulk

# Load N-body snapshots (GENGA .aei format)
A_dict = ad.load_aei_snapshots("examples/dataset/N_body_dataset/")

# Trace collision history and compute embryo compositions at 5 Myr
tracer   = ad.CollisionTracer(A_dict)
df_embryo = tracer.get_composition_at(t_myr=5.0)

# Classify impact events and assign equilibration pressure
processor = ad.ImpactEventProcessor(df_embryo)
data_dict = processor.run()

# Run multi-stage differentiation
df_result = ad.run_differentiation(data_dict)
```

> 💡 For the full annotated walkthrough, open the notebooks in [`examples/quickstart/`](examples/quickstart/).

---

## Physical Models

### Metal–Silicate Equilibration

AccreDiff implements pressure- and temperature-dependent partition coefficients
following the approach of Fischer et al. (2015) and Rubie et al. (2015).
Equilibration depth is treated as a free parameter or derived from impact scaling laws.

### Redox Evolution (ΔIW)

Oxygen fugacity is tracked relative to the Iron–Wüstite buffer (ΔIW).
The framework supports both fixed-ΔIW and evolving-ΔIW scenarios.

### Accretion History Input

Currently supported formats:
- **GENGA** binary/text outputs
- Custom CSV impact logs (see `data/format_spec.md`)

---

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

This framework draws on methodologies from:

- Rubie et al. (2011), *Earth and Planetary Science Letters* — Heterogeneous accretion, composition and core–mantle differentiation of the Earth
- Rubie et al. (2015), *Icarus* — Accretion and differentiation of the terrestrial planets with  implications for the compositions of early-formed Solar  System bodies and accretion of water
- Fischer et al. (2017), *Earth and Planetary Science Letters* — Sensitivities of Earth’s core and mantle compositions to accretion and differentiation processes
- Kong et al. (2026),*Astronomy and Astrophysics* - Coupling Dynamical Accretion and Chemical Differentiation: A Unified Framework for Earth--Mars Diversity (submit)

---

*AccreDiff is actively developed. Feedback, issues, and pull requests are warmly welcomed.*