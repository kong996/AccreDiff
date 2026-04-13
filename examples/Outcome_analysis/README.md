# Outcome Analysis - AccreDiff Examples

## Overview

This directory contains example Jupyter notebooks for analyzing and comparing planetary compositions simulated by the AccreDiff core formation model.

## Contents

### A_bulk_comparison.ipynb

**Purpose:** Compare simulated planetary bulk compositions with Earth's reference composition.

**Key Features:**
- Load final planetary compositions from core formation simulations
- Normalize compositions by MgO content
- Compare to Earth and Mars reference compositions
- Visualize multi-element abundance comparisons
- Export results to CSV format for further analysis

**Main Steps:**
1. Load core formation model output (Excel format)
2. Extract final timestep compositions for each planet
3. Normalize all elements by MgO content
4. Apply Earth reference normalization
5. Generate scatter plot with log-scale visualization
6. Export processed data to CSV

**Output:**
- CSV file containing normalized planetary compositions
- Visualization showing planetary compositions relative to Earth

**Input Files Required:**
- `../Tables/{Model}_core_formation.xlsx` (Core formation model output)

**Output Files Generated:**
- `./Tables/{Model}.csv` (Normalized composition data)

---

### B_assemble.ipynb

**Purpose:** Aggregate and compare planetary parameters across multiple simulations to classify Earth-like and Mars-like analogs.

**Key Features:**
- Load multiple CSV simulation files from ./Tables folder
- Classify simulations by impactor mass fraction (m_e)
- Calculate Core Mass Fraction (CMF) and Mg# for each population
- Generate comparative strip plot visualization
- Statistical analysis (mean ± std) for both populations
- Export processed results to Excel format

**Main Steps:**
1. Load all CSV files from simulation output folder
2. Classify simulations:
   - **Earth analogs**: m_e between 0.7 and 1.3
   - **Mars analogs**: m_e between 0.05 and 0.3
   - **Others**: Remaining cases
3. Aggregate classified data into two populations
4. Calculate key parameters:
   - Core Mass Fraction (CMF) using `ad.compute_core_ratio()`
   - Mg# (Magnesium number) = Mg/(Mg+Fe) in mantle
5. Compute statistical summaries (mean & standard deviation)
6. Create three-metric strip plot visualization:
   - CMF (Core Mass Fraction)
   - Mg# (Magnesium number)
   - Mass (Relative to Earth)
7. Export aggregated results to Excel

**Visualization Details:**
- **Strip Plot Layout**: Three horizontal bands (CMF, Mg#, Mass)
- **Data Points**: Scatter plots with jitter (individual simulations)
- **Statistical Markers**: Diamond symbols with error bars (mean ± 1σ)
- **Reference Values**: Dashed lines showing Earth and Mars reference values
- **Uncertainty Bands**: Semitransparent colored regions (±10% uncertainty)
- **Nature Style**: High-resolution publication-ready formatting

**Input Files Required:**
- Multiple CSV files in `./Tables/` folder (from A_bulk_comparison.ipynb output)
- Each CSV must contain columns: `m_e`, `MgO`, `FeO`, `Al2O3`, `CaO`, `SiO2`, `NiO`, `Fe`, `Ni`, `Si`, `O`

**Output Files Generated:**
- `./Tables/assemble_results.xlsx` (Excel file with two sheets):
  - Sheet 1: `Earth_analogs` - Classified Earth-like simulations with CMF and Mg#
  - Sheet 2: `Mars_analogs` - Classified Mars-like simulations with CMF and Mg#
  - Both sheets include mean and std rows for statistical reference

**Parameters (Modifiable):**
```python
# File path
Path = './Tables/'

# Classification boundaries
Earth_like: m_e ∈ [0.7, 1.3]
Mars_like:  m_e ∈ [0.05, 0.3]

# Visualization parameters
H = 1.0    # Band height
G = 0.15   # Gap between bands

# Reference values (editable in code)
ref = {
    "Mass": {"Earth": 1.00, "Mars": 0.1074},
    "Mg#":  {"Earth": 0.89, "Mars": 0.775},    
    "CMF":  {"Earth": 0.32, "Mars": 0.20},
}

# Reference uncertainty bands (±10% range)
bands_ref_range = {
    "Mg#":  {"Earth": (0.801, 0.979), "Mars": (0.698, 0.853)},
    "CMF":  {"Earth": (0.288, 0.352), "Mars": (0.18, 0.22)},
}
```

---

### C_source.ipynb

**Purpose:** Analyze and compare source material composition (EF, EC, OC, CI fractions) across different collision models for Earth-like and Mars-like planetary bodies.

**Key Features:**
- Load multiple simulation data files (CSV format with `_source.csv` suffix)
- Classify simulations by planetary type (Earth_like / Mars_like)
- Calculate statistical summaries (mean ± standard deviation)
- Generate comparative visualization with scatter plots and error bars
- Export aggregated results to Excel format with multiple sheets

**Main Steps:**
1. Load all `_source.csv` files from simulation output directory
2. Extract and classify compositions by 'Class' column:
   - **Earth analogs**: Class = 'Earth_like'
   - **Mars analogs**: Class = 'Mars_like'
3. Aggregate classified data into two populations
4. Calculate key composition parameters:
   - EF (Enstatite Chondrite Fraction)
   - EC (Carbonaceous Chondrite Fraction)
   - OC (Ordinary Chondrite Fraction)
   - CI (CI Chondrite Fraction)
5. Compute statistical summaries (mean & standard deviation) for each population
6. Create comparative scatter plot visualization:
   - Individual data points for each composition parameter
   - Mean ± 1σ error bars for statistical reference
   - Side-by-side comparison of Earth and Mars analogs
7. Export aggregated results to Excel with multiple sheets

**Visualization Details:**
- **Plot Type**: Scatter plot with error bars (mean ± std)
- **Data Points**: Jittered scatter showing individual simulations
  - Earth analogs: Blue circles (α = 0.5, s = 20)
  - Mars analogs: Red squares (α = 0.5, s = 30)
- **Statistical Markers**: Mean ± standard deviation with error caps
  - Earth analogs: Blue circle markers (12pt)
  - Mars analogs: Red square markers (12pt)
- **Layout**: Four composition parameters (EF, EC, OC, CI) on x-axis
- **Nature Style**: High-resolution publication-ready formatting
  - Serif font (font.size = 12)
  - Axis line width = 1
  - No top/right spines
  - DPI = 300 for export quality

**Input Files Required:**
- Multiple `*_source.csv` files in `../Accrection/Tables/` directory (from source material tracking output)
- Each CSV must contain columns:
  - `Class`: Classification category ('Earth_like' or 'Mars_like')
  - `EF`, `EC`, `OC`, `CI`: Composition fractions
  - (Optional) `model`: Source model identifier

**Output Files Generated:**
- `./Tables/M03_source_composition.xlsx` (Excel file with four sheets):
  - Sheet 1: `Earth_analogs` - All Earth-like simulations with composition data
  - Sheet 2: `Mars_analogs` - All Mars-like simulations with composition data
  - Sheet 3: `Earth_stats` - Statistical summary (mean and std rows)
  - Sheet 4: `Mars_stats` - Statistical summary (mean and std rows)

**Parameters (Modifiable):**
```python
# Data directory path
data_dir = '../Accrection/Tables/'

# Input file filter
File pattern: '*_source.csv' (case-sensitive)

# Classification categories
Class column values:
  - 'Earth_like' → Earth analogs
  - 'Mars_like'  → Mars analogs

# Composition parameters to analyze
params = ['EF', 'EC', 'OC', 'CI']

# Visualization configuration
figsize = (6, 4)
dpi = 150 (display), 300 (export)
marker_size_earth = 12
marker_size_mars = 12
scatter_alpha = 0.5
jitter_range = [0.05, 0.15]

# Color scheme
Colors = {
    'Earth': "#1f77b4" (Blue),
    'Mars': "#d62728" (Red),
}
```

**Last Updated:** 2026-04-13  
**Version:** 1.0  
**Status:** Ready for use