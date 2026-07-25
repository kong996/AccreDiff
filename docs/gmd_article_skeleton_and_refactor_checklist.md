# GMD article skeleton and project refactor checklist

Working title:

> AccreDiff v0.4.0: a Python framework coupling N-body accretion histories with metal-silicate differentiation for terrestrial planet composition modelling

## 1. Positioning

The previous A&A research paper, "Coupling dynamical accretion and chemical differentiation: A unified framework for the diversity of Earth and Mars", should serve as the scientific demonstration case for the GMD manuscript. The GMD paper should not simply repeat the Earth-Mars science result. Its main contribution should be the model/framework description:

- what AccreDiff computes;
- which equations and numerical assumptions control the output;
- how N-body collision histories are transformed into geochemical evolution;
- how users can reproduce a benchmark workflow;
- how the model is verified, evaluated, and versioned.

Recommended GMD manuscript type: model description paper. This is appropriate because GMD accepts descriptions of full models, model components, frameworks, and geoscientific software toolboxes. The title should include the model name and version number.

## 2. Core message

AccreDiff is a modular Python framework for tracing the coupled dynamical and chemical evolution of terrestrial planets. It links N-body accretion histories to impact-regime-dependent metal-silicate equilibration and predicts final mantle and core compositions, redox-sensitive element partitioning, Mg#, and core mass fraction.

The Earth-Mars application from the previous paper can be framed as the first benchmark demonstration:

- Earth analogs preferentially accrete reduced inner-disk reservoirs and experience deeper equilibration.
- Mars analogs are scattered outward earlier, sample more oxidized material, and equilibrate at shallower conditions.
- The divergence in FeO, Mg#, and CMF emerges from accretion provenance, radial redox structure, and impact-controlled equilibration.

## 3. Proposed manuscript skeleton

### Abstract

State the modelling gap: N-body simulations usually provide realistic accretion histories but not self-consistent chemical differentiation, while geochemical differentiation models often prescribe simplified growth histories.

Introduce AccreDiff as a Python framework that couples collision histories, radial compositional reservoirs, impact regimes, pressure-temperature conditions, and multistage metal-silicate equilibration.

Report the benchmark: reproducing the workflow used for Earth-Mars analogs, including source fractions, mantle/core compositions, Mg#, and CMF.

End with availability: versioned source code, reproducible examples, notebooks, and benchmark outputs.

### 1. Introduction

- Scientific motivation: Earth and Mars differ in mass, mantle FeO, Mg#, and CMF.
- Methodological gap: dynamical simulations and geochemical differentiation are often separated.
- Need for an open, modular, reproducible tool.
- Relationship to the previous A&A paper: that paper demonstrates a scientific application; this paper documents the model, equations, software architecture, and reproducibility package.

### 2. Model overview

Describe AccreDiff as a pipeline:

```text
N-body snapshots and collision records
  -> collision genealogy and source reservoirs
  -> early embryo bulk compositions
  -> impact regime classification
  -> pressure-temperature equilibration conditions
  -> metal-silicate partitioning solver
  -> time-dependent mantle/core compositions
  -> final planetary diagnostics and figures
```

Suggested figure: a software/model architecture diagram linking modules:

- `constant.py`: reference data and endmember compositions
- `chemistry.py`: molar masses, KD, IW, CMF
- `accretion.py`: collision tracing, growth fitting, impact processing
- `differentiation.py`: KD solver and multistage core formation
- `melt_model.py`: impact-generated melt pool and magma-ocean scaling
- `gas_migration.py`: optional gas-disk migration utilities
- `examples/quickstart`: reproducible user workflows

### 3. Input data model

Define required inputs clearly:

- N-body snapshots: particle IDs, masses, semi-major axes, eccentricities, inclinations, times.
- Collision records: time, target/impactor IDs, target/impactor masses, product IDs.
- Initial compositional reservoirs: EF, EC, OC, CI or user-defined alternatives.
- Model parameters: gas disk lifetime, impact thresholds, equilibration fractions, P-T relations, KD coefficients.

This section should also define current limitations in supported formats. If GENGA `.aei` files are advertised, the loader must exist as a stable public API or the README must be corrected.

### 4. Accretion and source tracing

Describe:

- collision product ID resolution;
- full collision genealogy with `CollisionTracer`;
- mass-weighted compositional mixing;
- source fraction reconstruction from initial semi-major axis;
- growth curve fitting and timescale diagnostics.

Benchmark figure candidates:

- cumulative source fraction versus initial semi-major axis;
- EF/EC/OC/CI source fractions for Earth and Mars analogs;
- mass evolution through time.

### 5. Impact regime classification

Translate the A&A impact-regime framework into software documentation:

- early phase impacts, `t <= 5 Myr`;
- late small impacts;
- late large impacts and post-large-impact magma-ocean accretion;
- participating fractions of impactor core, impactor mantle, target core, target mantle;
- equilibration pressure definitions P1, P2, P3, P4;
- equilibration temperature from Rubie-style pressure-temperature relations.

Suggested figure: simplified version of the A&A Fig. 3, but redrawn as a model-process schematic rather than a science-result figure.

### 6. Metal-silicate differentiation solver

Document the mathematical core of AccreDiff:

- mass conservation for Fe, Ni, Si, and O;
- Mg, Al, and Ca as lithophile components in the silicate phase;
- KD formulation and P-T dependence;
- oxygen fugacity / IW calculation;
- `ForwardKDOSolver` algorithm and convergence behavior;
- how partial equilibration updates target and impactor reservoirs through time.

Verification tests should be referenced here:

- mass conservation residuals;
- KD residuals;
- non-negativity and physical bounds;
- regression against known single-event solutions.

### 7. Reference workflow and user interface

Use the quickstart notebooks as the user-facing workflow:

1. Constants and setup
2. Meteorite bulk composition
3. Embryo bulk composition
4. Accretion history
5. Impact events
6. Differentiation
7. Bulk composition comparison

For GMD, this must be made reproducible without hidden local paths. Ideally provide:

- a small packaged toy dataset;
- a complete benchmark dataset used in the manuscript;
- scripts or notebooks that regenerate every GMD figure;
- stable outputs in `examples/.../outputs` or `benchmarks/...`.

### 8. Benchmark demonstration: Earth-Mars analogs

Use the previous A&A paper as the scientific evaluation case, but keep the emphasis on model capability and reproducibility.

Recommended outputs:

- final orbital architectures or selected system summary;
- source fractions: Earth analogs enriched in EF/EC, Mars analogs enriched in OC/CI;
- mantle/core compositions normalized to Mg and references;
- mass, Mg#, and CMF comparison;
- a table comparing benchmark outputs against Earth/Mars reference values.

Keep the wording careful: this is an evaluation of model behavior and plausibility, not a fresh claim that needs to duplicate all of the A&A interpretation.

### 9. Verification and evaluation

Separate two questions explicitly:

- Verification: does the code solve the intended equations correctly?
- Evaluation: does the framework produce plausible planetary compositions in a realistic use case?

Minimum verification package:

- import test;
- molar-mass calculator tests;
- KD calculator tests at fixed P-T values;
- `ForwardKDOSolver` mass conservation tests;
- impact-regime classification tests;
- one miniature end-to-end workflow test.

Minimum evaluation package:

- reproduce the A&A Earth-Mars benchmark table;
- compare Mg#, CMF, mantle FeO, and selected core light elements against reference values;
- report known offsets: Earth core O, Mars mantle NiO, sulfur-free CMF interpretation.

### 10. Discussion

Suggested subsections:

- Relationship to existing coupled accretion-differentiation models.
- Why modularity matters: swapping reservoirs, impact physics, or N-body inputs.
- Scope of applicability: terrestrial planets, rocky exoplanets, Earth-Mars analogs.
- Limitations:
  - sulfur not included;
  - volatile loss not included;
  - fixed/parameterized equilibration fractions;
  - perfect merger assumption inherited from collision treatment;
  - only selected major elements in the current framework;
  - dependence on N-body input quality.

### 11. Conclusions

Focus on what the software enables:

- AccreDiff provides an open bridge between dynamical accretion histories and geochemical differentiation.
- The modular design exposes each modelling assumption.
- The Earth-Mars benchmark demonstrates the complete workflow.
- The code is ready for extension to sulfur, volatiles, trace elements, and alternative formation scenarios.

### Code and data availability

For GMD readiness, this section should include:

- exact model version;
- GitHub repository URL;
- Zenodo DOI or another persistent archive for the exact release;
- license;
- installation instructions;
- benchmark data DOI or supplement;
- notebooks or scripts needed to reproduce manuscript figures.

## 4. Figure plan

1. AccreDiff workflow diagram.
2. Module architecture and data flow.
3. Impact-regime schematic.
4. Solver verification: residuals / mass conservation / KD match.
5. Benchmark source fractions for Earth and Mars analogs.
6. Benchmark composition comparison: mantle/core spider or normalized ratios.
7. Bulk diagnostics: mass, Mg#, CMF.
8. Optional: runtime or scaling summary for workflow stages.

## 5. Project refactor checklist

### A. Public API and documentation

- [ ] Decide publication version: `v0.4.0`, `v0.5.0`, or `v1.0.0-gmd`.
- [ ] Align `pyproject.toml`, `src/accrediff/__init__.py`, package metadata, and release tags.
- [ ] Fix README examples: `load_aei_snapshots()` and `run_differentiation()` are advertised but not currently exported.
- [ ] Add a clear "Model equations" document matching the manuscript.
- [ ] Add a user manual or `docs/user_guide.md`.
- [ ] Add an API reference page for public classes and functions.

### B. Installation and environment

- [ ] Ensure `pip install -e .` works from a clean environment.
- [ ] Add `python >= 3.10` consistently in README and requirements.
- [ ] Remove or regenerate stale `src/accrediff.egg-info` metadata.
- [ ] Provide a minimal environment file, for example `environment.yml`.
- [ ] Avoid import-time matplotlib/font cache side effects from top-level imports where possible.

### C. Data and reproducibility

- [ ] Add a small toy dataset for CI and documentation.
- [ ] Add a benchmark dataset for the GMD paper, or publish it separately with DOI.
- [ ] Make all notebook inputs relative to the repository or documented data paths.
- [ ] Add scripts that regenerate manuscript figures from raw/intermediate data.
- [ ] Add a manifest listing all benchmark inputs, outputs, and checksums.

### D. Testing and verification

- [ ] Expand beyond `tests/test_import.py`.
- [ ] Test `MolarMassCalculator`.
- [ ] Test `KDCalculator` at fixed P-T values.
- [ ] Test `ForwardKDOSolver` mass conservation and residual tolerance.
- [ ] Test `ImpactEventProcessor` event classification.
- [ ] Test source-fraction reconstruction with a tiny collision tree.
- [ ] Add a smoke test for the full quickstart workflow on toy data.

### E. Scientific validation package

- [ ] Define the exact Earth-Mars benchmark case inherited from the A&A paper.
- [ ] Reproduce source fractions for seven selected systems or a documented subset.
- [ ] Reproduce mean mantle/core compositions.
- [ ] Reproduce mass, Mg#, and CMF diagnostics.
- [ ] Document known discrepancies and why they are expected: sulfur-free Mars, volatile-free model, core O uncertainty.

### F. Software engineering cleanup

- [ ] Standardize naming: `Accrection` should probably become `Accretion` in examples, or be clearly retained for backward compatibility.
- [ ] Decide whether notebooks under `notebooks/` are archival, active, or deprecated.
- [ ] Separate core package code from paper-specific analysis.
- [ ] Move reusable notebook helper code into package modules.
- [ ] Add type hints and docstrings for public APIs.
- [ ] Add logging or structured output for long workflows.

### G. GMD submission assets

- [ ] Create manuscript repository or `paper/` directory.
- [ ] Add Copernicus/GMD LaTeX or Word template.
- [ ] Prepare `Code availability`, `Data availability`, and `Author contributions`.
- [ ] Archive release on Zenodo and reserve DOI before submission.
- [ ] Prepare a short summary for GMD submission.
- [ ] Ensure all figures are 300 dpi or vector and colorblind-friendly.

## 6. Recommended next work sequence

1. Freeze the model scope for the GMD paper: major elements only, sulfur/volatiles as future work.
2. Fix package install/import and version metadata.
3. Create the toy dataset and miniature end-to-end workflow test.
4. Convert the A&A benchmark into reproducible scripts/notebooks.
5. Write the Methods and Verification sections before polishing the Introduction.
6. Archive code and benchmark data only after all figures can be regenerated cleanly.

