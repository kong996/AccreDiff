# AccreDiff

AccreDiff is an open-source Python framework for coupling N-body accretion
histories with multi-stage core-mantle differentiation calculations for
terrestrial planet formation.

The project is developed as a reproducible model-workflow package for the GMD
model-description paper:

> AccreDiff v1.0.0: a modular framework for coupling planetary accretion and
> core-mantle differentiation

## Main Capabilities

AccreDiff provides tools to:

- read and process GENGA-like N-body accretion outputs;
- reconstruct collision histories, parent-product relationships, and source
  reservoir fractions;
- initialize cosmochemical compositions from CI-based or meteoritic-reservoir
  prescriptions;
- assign impact regimes and pressure-temperature equilibration conditions;
- perform event-by-event metal-silicate differentiation;
- track mantle and core inventories, redox state, bulk composition, and core
  mass fraction;
- run extension examples for gas-disk migration, collision-outcome
  classification, and Fe-Ni-S sulfur-bearing differentiation.

AccreDiff is a post-processing and coupled geochemical evolution framework. It
does not recompute N-body dynamics and does not solve hydrodynamic impact
physics directly.

## Installation

AccreDiff requires Python 3.10 or later.

Clone the repository:

```bash
git clone https://github.com/kong996/AccreDiff.git
cd AccreDiff
```

Install the package in editable mode:

```bash
pip install -e .
```

For development and example notebooks, install the optional dependencies:

```bash
pip install -e ".[dev,examples]"
```

Alternatively, create a conda environment:

```bash
conda env create -f environment.yml
conda activate accrediff
```

## Quick Verification

Run the test suite from the repository root:

```bash
python -m pytest -q
```

Build the source and wheel distributions:

```bash
python -m build --sdist --wheel
```

Check the installed package version:

```bash
python -c "import accrediff as ad; print(ad.__version__)"
```

For the GMD release described in the paper, this should print:

```text
1.0.0
```

## Project Structure

```text
AccreDiff/
├── src/accrediff/                    # Core Python package
├── examples/                         # Reproducible example workflows
│   ├── quickstart/                   # Seven-step end-to-end workflow
│   ├── Cosmochemistry/               # Cosmochemical initialization examples
│   ├── Accrection/                   # Accretion reconstruction examples
│   ├── Impact/                       # Impact-regime examples
│   ├── Differentiation/              # Core-mantle differentiation examples
│   ├── Outcome_analysis/             # Diagnostic post-processing examples
│   ├── Expansion_module/             # Extension-module demonstrations
│   └── dataset/                      # Bundled demonstration datasets
├── tests/                            # Unit tests
├── docs/                             # Manuscript-support notes and documents
├── CITATION.cff                      # Software citation metadata
├── .zenodo.json                      # Zenodo archive metadata
├── CHANGELOG.md                      # Release notes
├── environment.yml                   # Conda environment for examples/tests
├── pyproject.toml                    # Package metadata
└── README.md
```

The directory name `examples/Accrection/` is retained for compatibility with
the distributed notebooks and manuscript paths.

## Example Workflows

The recommended starting point is the quickstart workflow. The notebooks are
designed to be run in order:

| Notebook | Purpose |
| --- | --- |
| `examples/quickstart/01_Constants_and_setup.ipynb` | Load constants and reference compositions |
| `examples/quickstart/02_Meteorite_bulk_composition.ipynb` | Build meteoritic endmember compositions |
| `examples/quickstart/03_Embryo_bulk_composition.ipynb` | Reconstruct embryo compositions |
| `examples/quickstart/04_Accrection_history.ipynb` | Reconstruct planetary accretion histories |
| `examples/quickstart/05_Impact_events.ipynb` | Assign impact regimes and equilibration pressure |
| `examples/quickstart/06_Differentiation.ipynb` | Run multi-stage core-mantle differentiation |
| `examples/quickstart/07_Comparison_bulk.ipynb` | Compare final model outputs with reference compositions |

Additional module-specific examples are provided in:

- `examples/Cosmochemistry/`
- `examples/Accrection/`
- `examples/Impact/`
- `examples/Differentiation/`
- `examples/Outcome_analysis/`
- `examples/Expansion_module/gas_disk/`
- `examples/Expansion_module/collision_detail/`
- `examples/Expansion_module/Fe_Ni_S_ternary/`

Each example folder includes a README describing its inputs, outputs, and main
AccreDiff APIs.

## Minimal Usage

```python
import accrediff as ad

print(ad.__version__)

elements = ad.constants.Elements
earth = ad.constants.Earth
ci = ad.constants.CI_bulk

mm = ad.MolarMassCalculator()
mg_o_mass = mm.molar_mass("MgO")
print(mg_o_mass)
```

For complete workflows, use the notebooks in `examples/quickstart/`.

## Extension Modules

AccreDiff v1.0.0 includes prototype extension modules that demonstrate the
intended modular interface:

- `src/accrediff/gas_migration.py`: Python-side gas-disk migration and torque
  diagnostics, with a GENGA-compatible CUDA configuration example in
  `examples/dataset/gas_config/gas.cu`.
- `src/accrediff/collision_outcomes.py`: EDACM-style collision-outcome
  diagnostics for GENGA-like collision records.
- `src/accrediff/differentiation_sulfur.py`: standalone Fe-Ni-S
  metal-silicate differentiation extension with Rose-Weston and Boujibar sulfur
  partitioning examples.

These modules are distributed as extension examples rather than mandatory parts
of every default workflow.

## Citation

If you use AccreDiff, please cite the archived software release and the
accompanying GMD model-description paper. Citation metadata are provided in
`CITATION.cff`.

After creating the GitHub release and Zenodo archive, cite the exact archived
version, for example:

```text
Kong, Z.: AccreDiff v1.0.0: a modular framework for coupling planetary
accretion and core-mantle differentiation, Zenodo [code],
https://doi.org/10.5281/zenodo.XXXXXXX, 2026.
```

## GMD Release Checklist

Before submitting the GMD discussion paper, the release should be frozen and
archived:

1. Confirm that `pyproject.toml` and `accrediff.__version__` both report
   `1.0.0`.
2. Run `python -m pytest -q`.
3. Run `python -m build --sdist --wheel`.
4. Create a GitHub release tagged `v1.0.0`.
5. Archive the release on Zenodo and obtain a DOI.
6. Add the Zenodo DOI to the manuscript code-availability section and, if
   desired, to `CITATION.cff`.

## License

AccreDiff is distributed under the MIT License. See `LICENSE` for details.

## Acknowledgements

AccreDiff draws on methods and parameterizations from planetary accretion,
metal-silicate partitioning, impact physics, and cosmochemical reservoir
studies. Key scientific references are discussed in the accompanying GMD
manuscript and in the example notebooks.
