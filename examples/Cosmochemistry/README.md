# Cosmochemistry Examples

This folder contains example notebooks for computing bulk planetary compositions
under different cosmochemical assumptions, built on top of the `accrediff` package.

## Notebooks

### 1. `A_IW15IW35.ipynb` — Oxygen Fugacity Endmember Compositions

Computes bulk compositions for two oxygen fugacity endmembers (**IW-1.5** and
**IW-3.5**) via forward core–mantle differentiation modeling, starting from the
Bulk Earth Reference enriched in refractory elements (Rubie et al., 2011).

**Key steps:**
- Compute metal–silicate partition coefficients (Ni, Si, O) at a reference P–T
  condition using `KDCalculator`
- Convert the Bulk Earth Reference (wt%) to molar amounts using `MolarMassCalculator`
- Search for the oxygen loss parameter (O_L) that satisfies each target IW buffer
  via `OLSolver`
- Forward-solve core–mantle element partitioning with `ForwardKDOSolver`
- Validate calculated compositions against Fischer et al. (2017) reference models
  with a reusable two-panel plot (bar comparison + residuals with MAE/RMSE/MAPE)
- Export MgO-normalized molar compositions to CSV

**Output:**
- `Tables/A1_IW35IW15_composition.csv`

---

### 2. `B_Meteorite.ipynb` — Chondrite Endmember Compositions

Establishes bulk compositions for four chondrite groups (**EH/EF, EL/EC, OC, CI**)
and compares them with **Earth** and **Mars** bulk compositions.

**Key steps:**
- Load meteorite and planetary mantle/core compositions from `ad.constants`
- Compute compound molar masses for oxides & metals using `MolarMassCalculator`
- Convert weight percent (wt%) to molar amounts via a reusable `wt_to_molar()` helper
- Normalize all molar amounts to MgO using `normalize_to_mgo()`
- Apply refractory element enrichment (Al₂O₃, CaO ×1.5) and volatile element
  depletion (SiO₂, Si ×0.6) corrections for the EF endmember
- Double-normalize to Mg and CI chondrites for cross-system comparison
  (mantle oxides and core metals mapped via `core_oxide_map`)
- Plot CI-chondrite-normalized spider diagram (mantle ■ vs. core ●)
- Export MgO-normalized molar compositions to CSV

**Output:**
- `Tables/A2_Meteorites_composition.csv`

