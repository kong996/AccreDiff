## 📓 Quickstart Notebooks

The notebooks are designed to be run **in order**. Each notebook builds on the previous one.

---

### 01 - Constants and Setup

**File:** `01_Constants_and_setup.ipynb`

The first notebook in the quickstart series. It loads and inspects all built-in constants and data structures provided by `accrediff`, and sets up the shared plotting style used throughout the series.

**What you will learn:**
- How to configure the plot style for publication-quality figures
- How to access the element atomic mass table (`Elements`)
- How to use `MolarMassCalculator` to compute molar masses of oxides and metals
- How to load meteorite bulk compositions (`EF`, `EC`, `OC`, `CI`) in wt%
- How to load planet mantle & core compositions (`Earth`, `Mars`) in wt%
- How to calculate bulk compositions of Earth & Mars from mantle/core data

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `Elements` | Atomic/molar mass table (updated with oxide compounds) |
| `Colors` | Color palette for Earth, Mars, EF, EC, OC, CI |
| `EF_bulk`, `EC_bulk`, `OC_bulk`, `CI_bulk` | Meteorite bulk compositions (wt%) |
| `Earth_mantle`, `Earth_core` | Earth mantle & core compositions (wt%) |
| `Mars_mantle`, `Mars_core` | Mars mantle & core compositions (wt%) |
| `Earth_bulk`, `Mars_bulk` | Calculated planet bulk compositions (wt%) |

> ➡️ **Next:** `02_Meteorite_bulk_composition.ipynb`

---

### 02 - Meteorite Bulk Composition

**File:** `02_Meteorite_bulk_composition.ipynb`

This notebook processes meteorite bulk compositions step by step, from raw weight percentages to CI-normalized spider diagrams, and exports the results to CSV.

**What you will learn:**
- How to convert wt% to molar amounts using element/oxide molar masses
- How to normalize compositions to MgO as a reference element
- How to apply geochemical enrichment/depletion corrections to EF composition
- How to normalize to both Mg and CI chondrites for cross-system comparison
- How to visualize mantle and core compositions as a CI-normalized spider diagram
- How to export processed compositions to CSV

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `EF/EC/OC/CI_bulk_molar` | Molar amounts converted from wt% |
| `EF/EC/OC/CI_bulk_mg` | MgO-normalized molar ratios |
| `mantle_CI` | CI-normalized mantle compositions for all systems |
| `core_CI` | CI-normalized core compositions for all systems |
| `./Tables/Meteorites_composition.csv` | Exported bulk composition table |

> ➡️ **Next:** `03_Embryo_bulk_composition.ipynb`

---

### 03 - Embryo Bulk Composition (~5 Myr)

**File:** `03_Embryo_bulk_composition.ipynb`

This notebook establishes the bulk composition of planetary embryos at the end of the gas disk phase (~5 Myr) by tracing collision histories from N-body simulations.

**What you will learn:**
- How to load N-body simulation snapshots (`aei` files) and convert mass to Earth units
- How to classify initial particle compositions by semi-major axis (EF / EC / OC / CI)
- How to load collision records and filter events within the gas disk lifetime (≤ 5 Myr)
- How to use `CollisionTracer` to reconstruct compositional mixing from full merger histories
- How to compute meteorite mass fractions and thermodynamic conditions (P, T)
- How to use `EarlyComUpdater` to assign geochemical bulk compositions to each particle
- How to filter embryos (m ≥ 0.05 M⊕) and visualize CI-normalized spider diagrams
- How to export the full particle dataset to CSV

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `A_dict` | N-body snapshots dict keyed by time string |
| `df_com` | Initial particle orbital & composition data at t = 0 Myr |
| `i_com_mapping` | Particle index → composition type mapping |
| `C_withgas` | Collision records within gas disk lifetime |
| `df_gas_all` | All particles at t = 5 Myr with composition, P, T |
| `df_embryos` | Filtered embryos (m ≥ 0.05 M⊕), sorted by mass |
| `./Tables/03_embryo_bulk_composition.csv` | Exported full particle dataset |

> ➡️ **Next:** `04_Accrection_history.ipynb`

---

### 04 - Accretion History (~5 Myr → 100 Myr)

**File:** `04_Accrection_history.ipynb`

This notebook reconstructs the post-gas-disk accretion history of planetary embryos by tracking mass evolution and tracing collision events across the full N-body simulation timeline.

**What you will learn:**
- How to load N-body snapshots and identify final planets (m ≥ 0.05 M⊕) at t = 100 Myr
- How to track each planet's mass and semi-major axis evolution across all time steps
- How to visualize full mass evolution (0 → 100 Myr) with bubble-size scatter plots
- How to load collision records and filter events after the gas disk phase (t > 5 Myr)
- How to use `CollisionTracer` to reconstruct complete merger trees for each planet
- How to visualize post-disk accretion events on a log-scale timeline
- How to export per-planet collision histories to Excel

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `A_dict` | N-body snapshots dict keyed by time string |
| `df_planet` | Final planets (m ≥ 0.05 M⊕) at t = 100 Myr |
| `planets_index` | Particle indices of final planets |
| `df_history` | Mass & semi-major axis evolution for each planet |
| `collision_df` | Post-disk collision records with `product_id` |
| `dict_history` | Per-planet full collision history DataFrames |
| `./Tables/04_accretion_history.xlsx` | Exported collision history (one sheet per planet) |

> ➡️ **Next:** `05_Impact_events.ipynb`

---

### 05 - Impact Events & Melting Conditions

**File:** `05_Impact_events.ipynb`

This notebook classifies each post-disk collision event by melting type and establishes the corresponding metal–silicate equilibration pressure (P) conditions throughout planetary accretion.

**What you will learn:**
- How to load per-planet collision histories from `04_accretion_history.xlsx`
- How to use `ImpactEventProcessor` to classify impacts as **small**, **global**, or **partial** melting
- How to propagate global melting thermal effects using a cooling duration parameter
- How to identify partial melting episodes using `ad.update_partial_melting_process()`
- How to visualize mass and equilibration pressure evolution colored by event type
- How to export fully classified impact event tables to Excel

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `data_dict` | Per-planet DataFrames with `events` and `P_equil` columns |
| `events` column | `small` / `global` / `partial` melting classification |
| `P_equil` column | Metal–silicate equilibration pressure at each event (GPa) |
| `./Tables/05_impact_events.xlsx` | Exported impact event table (one sheet per planet) |

> ➡️ **Next:** `06_Differentiation.ipynb`

---

### 06 - Multi-stage Core Formation & Differentiation

**File:** `06_Differentiation.ipynb`

This notebook simulates the multi-stage metal–silicate differentiation of growing planets by iterating through each impact event and solving for chemical equilibrium between metal and silicate phases.

**What you will learn:**
- How to load early embryo bulk compositions from `03_embryo_bulk_composition.csv`
- How to load classified impact event histories from `05_impact_events.xlsx`
- How to partition bulk material into re-equilibrating and non-re-equilibrating fractions based on event type (`small` / `global` / `partial`)
- How to compute partition coefficients (K_D) for Ni, Si, and O using `KDCalculator`
- How to solve for metal–silicate equilibrium using `KD_Params` + `ForwardKDOSolver`
- How to iteratively update target planet compositions after each collision event
- How to collect final planet compositions and normalize to Mg and CI chondrites
- How to visualize final mantle and core compositions as a CI-normalized spider diagram
- How to export the full differentiation history to Excel

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `df_dict` | Per-planet differentiation history DataFrames (one per sheet) |
| `df_planets` | Final bulk compositions of all planets at t = 100 Myr |
| `df_planets_Mg` | MgO-normalized compositions |
| `df_planets_CI` | Mg- & CI-normalized compositions for spider diagram |
| `./Tables/06_differentiation.xlsx` | Exported differentiation history (one sheet per planet) |

> ➡️ **Next:** `07_Comparison_bulk.ipynb`

---

### 07 - Bulk Composition Comparison

**File:** `07_Comparison_bulk.ipynb`

This notebook compares the modeled final planet compositions against the bulk compositions of Earth and Mars, providing a direct observational benchmark for the accretion-differentiation model.

**What you will learn:**
- How to load the per-planet differentiation history from `06_differentiation.xlsx`
- How to extract final bulk compositions at t = 100 Myr from the last row of each planet's history
- How to normalize major element abundances to MgO
- How to load Earth and Mars reference compositions from `ad.constants`
- How to normalize modeled planets and Mars to Earth bulk values for cross-comparison
- How to visualize a spider diagram with Earth = 1 reference line and ±10% gray band
- How to export Mg-normalized planet compositions to CSV

**Key outputs:**

| Variable | Description |
|----------|-------------|
| `df_planets` | Final bulk compositions of all planets at t = 100 Myr |
| `df_planets_Mg` | MgO-normalized compositions |
| `df_Earth_bulk_mg` | Earth Mg-normalized reference composition |
| `df_Mars_bulk_mg` | Mars Mg-normalized reference composition |
| `df_planets_Earth` | Planet + Mars compositions normalized to Earth bulk |
| `./Tables/07_bulk_comparison.csv` | Exported Mg-normalized planet compositions |