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
---

### D_analysis.ipynb

**Purpose:** Analyze and compare planetary bulk compositions (normalized by reference standards) across Earth-like and Mars-like analog populations, with detailed meteorite reference data integration.

**Key Features:**
- Load and classify multiple simulation datasets by planetary mass (m_e)
- Normalize compositions by Earth and Mars reference standards
- Integrate meteorite composition reference data (EF, EC, OC, CI types)
- Generate multi-element composition comparison visualizations
- Calculate statistical summaries with log-scale error analysis
- Compute weight fractions and molar mass conversions
- Export statistical results for publication-quality analysis

**Main Steps:**
1. Load all simulation CSV files from `./Tables/` directory
2. Classify simulations by planetary mass ratio (m_e):
   - **Earth analogs**: m_e ∈ [0.7, 1.3]
   - **Mars analogs**: m_e ∈ [0.05, 0.3]
   - **Others**: Remaining cases
3. Load reference compositions:
   - Earth bulk composition (from `ad.constants.Earth_bulk_mg`)
   - Mars bulk composition (from `ad.constants.Mars_bulk_mg`)
   - Meteorite types (EF, EC, OC, CI) from external CSV file
4. Normalize all compositions by Mg content (Mg-normalized):
   - Compute common elemental columns
   - Divide by reference baseline values
   - Calculate mean ± standard deviation on log scale
5. Generate three visualization types:
   - **Combined Panel**: Earth-normalized (top) and Mars-normalized (bottom) side-by-side
   - **Earth Analogs**: Single panel Earth-normalized composition
   - **Mars Analogs**: Single panel Mars-normalized composition
6. Statistical analysis:
   - Calculate arithmetic mean for each element
   - Compute sample standard deviation (n-1)
   - Convert to log-scale statistics (log mean ± 1σ)
7. Weight fraction calculation:
   - Load molar masses for all elements and oxides
   - Convert atomic fractions to mass fractions
   - Normalize to total mass = 1.0
   - Compare simulated vs. reference values
8. Export statistical summary to console

**Visualization Details:**

**Panel Layout:**
- **Figure Size**: 12 × 6.75 inches (16:9 aspect ratio, publication-ready)
- **Subplots**: Two horizontal panels (Earth-normalized, Mars-normalized)
- **Title**: "Planetary Composition at 100 Myr"

**Data Representation:**
- **Simulation Points**: Gray scatter (n=total samples)
  - Jittered horizontally (offset range: ±0.2)
  - Size: 12pt, transparency: α = 0.85
  - Z-order: 3 (behind error bars)
  
- **Mean ± 1σ Error Bars**: Black diamond markers
  - Computed on log scale (geometric mean)
  - Error caps with width = 2pt
  - Size: 5pt, line width: 1.0pt
  - Z-order: 4 (in front of scatter)

- **Meteorite References**: Colored diamond markers
  - **EF** (Enstatite): Brown (#8c564b), offset: -0.25
  - **EC** (Carbonaceous): Orange (#ff7f0e), offset: -0.15
  - **OC** (Ordinary): Green (#2ca02c), offset: +0.15
  - **CI** (CI): Purple (#9467bd), offset: +0.25
  - Size: 50pt, edge: black 0.6pt, α = 0.8

- **Reference Bands**:
  - ±10% uncertainty shaded region (semi-transparent color band)
  - Reference line at normalized ratio = 1.0 (dashed, linestyle="--")
  - Color: Blue for Earth, Red for Mars

**Axis Configuration:**
- **Y-axis**: Logarithmic scale (log₁₀)
- **Y-limits**: 0.01 to 100 (Earth-normalized), 0.1 to 100 (Mars-normalized)
- **X-axis**: Element names (rotated 45°, right-aligned)
- **Grid**: Major and minor gridlines (dashed, α = 0.35)
- **Spines**: Top and right spines removed

**Legend:**
- **Location**: Right side, anchored outside plot area (bbox_to_anchor=(1.02, 0.5))
- **Format**: Vertical legend without frame
- **Items**:
  1. Simulation data point (gray circle, label: "simulation (n=X)")
  2. Mean ± 1σ (black diamond, label: "mean±1σ")
  3. Reference band (colored patch, label: "±10% [Earth/Mars]")
  4. Meteorite types (EF, EC, OC, CI as diamond markers)
- **Font Size**: 9pt
- **Spacing**: Tight (handletextpad=0.6, labelspacing=0.4)

**Style Parameters:**
- **Font Family**: Serif (Nature-style publication standard)
- **Font Size**: Main=12pt, axes labels=14pt, tick labels=12pt
- **Line Width**: Axes=1pt
- **DPI**: 300 (export quality), 150 (display)
- **Tick Direction**: Inward (xtick.direction="in", ytick.direction="in")

**Input Files Required:**
- Multiple CSV files in `./Tables/` directory (from A_bulk_comparison.ipynb output)
  - **File naming**: `{model_name}.csv` (case-sensitive)
  - **Required columns**: `m_e`, `MgO`, `Al2O3`, `CaO`, `FeO`, `NiO`, `SiO2`, `Fe`, `Ni`, `Si`, `O`
  - **Expected rows**: Individual planet simulations with composition data

- Meteorite composition reference file:
  - **Path**: `../Cosmochemistry/Tables/A2_Meteorites_composition.csv`
  - **Required rows**: EF, EC, OC, CI meteorite types
  - **Required columns**: Element or oxide names matching simulation data

**Output Files Generated:**
- Console printouts with statistical summaries
- Three visualization figures (displayed inline in Jupyter):
  1. Combined Earth + Mars panel (16:9 format)
  2. Earth analogs single panel (10 × 5 inches)
  3. Mars analogs single panel (10 × 5 inches)
- PNG export option (commented in code):
  ```python
  fig.savefig('composition_comparison_16x9.png', dpi=300, bbox_inches='tight')
  ```

**Parameters (Modifiable):**
```python
# ── Data Loading ────────────────────────────────────
data_folder = './Tables/'  # Path to CSV simulation files

# ── Classification Boundaries ───────────────────────
Earth_like: m_e ∈ [0.7, 1.3]
Mars_like:  m_e ∈ [0.05, 0.3]

# ── Reference Compositions ──────────────────────────
Earth_bulk_mg = ad.constants.Earth_bulk_mg
Mars_bulk_mg = ad.constants.Mars_bulk_mg
meteorites_file = '../Cosmochemistry/Tables/A2_Meteorites_composition.csv'

# ── Element Analysis ────────────────────────────────
cols = ['MgO', 'Al2O3', 'CaO', 'FeO', 'NiO', 'SiO2', 'Fe', 'Ni', 'Si', 'O']
mantle_items = ['MgO', 'Al2O3', 'CaO', 'FeO', 'NiO', 'SiO2']
core_items = ['Fe', 'Ni', 'Si', 'O']

# ── Meteorite Categories ────────────────────────────
meteorite_categories = ['EF', 'EC', 'OC', 'CI']
meteorite_offsets = {
    'EF': -0.25,  # Left offset for enstatite chondrites
    'EC': -0.15,  # Left offset for carbonaceous chondrites
    'OC':  0.15,  # Right offset for ordinary chondrites
    'CI':  0.25   # Far right offset for CI chondrites
}

# ── Normalization References ────────────────────────
ref_color_earth = 'blue'
ref_color_mars = 'red'
ref_uncertainty = 0.10  # ±10% band

# ── Visualization Configuration ─────────────────────
figsize_combined = (12, 6.75)      # 16:9 aspect ratio
figsize_single = (10, 5)           # Single panel
fontsize_title = 18                # Main title
fontsize_ylabel = 12               # Y-axis label
fontsize_legend = 9                # Legend text
scatter_size = 12                  # Simulation data points
meteorite_marker_size = 50         # Meteorite reference points
mean_marker_size = 5               # Mean ± 1σ markers
jitter_range = [-0.2, 0.2]         # Horizontal scatter range

# ── Scale Configuration ─────────────────────────────
y_scale = 'log'                    # Logarithmic y-axis
y_lim_earth = (0.01, 100)          # Y-limits for Earth-normalized
y_lim_mars = (0.1, 100)            # Y-limits for Mars-normalized

# ── Color Scheme ────────────────────────────────────
Colors = {
    'Earth': "#1f77b4",    # Blue
    'Mars': "#d62728",     # Red
    'EF': "#8c564b",       # Brown (Enstatite)
    'EC': "#ff7f0e",       # Orange (Carbonaceous)
    'OC': "#2ca02c",       # Green (Ordinary)
    'CI': "#9467bd",       # Purple (CI)
}

# ── Statistical Calculation ─────────────────────────
std_type = 'ddof=1'                # Sample standard deviation (n-1)
log_mean_method = 'geometric'      # Use geometric mean on log scale
```

**Workflow Integration:**
```
A_bulk_comparison.ipynb (PRIOR)
        ↓
    Produces: ./Tables/{model}.csv
        ↓
B_assemble.ipynb (OPTIONAL)
        ↓
    Classification reference
        ↓
D_analysis.ipynb (THIS NOTEBOOK) ← ANALYSIS & VISUALIZATION
        ↓
    Outputs: Statistics + Figures
```

**Python Functions Used:**
- `pd.read_csv()`: Load simulation data
- `pd.concat()`: Aggregate classified populations
- `pd.DataFrame.div()`: Normalize by reference composition
- `np.log10()`: Convert to log scale
- `ad.constants.Earth_bulk_mg`: Load Earth reference
- `ad.constants.Mars_bulk_mg`: Load Mars reference
- `ad.MolarMassCalculator()`: Calculate molar masses
- `plt.subplots()`: Create figure layout
- `ax.scatter()`: Plot individual data points
- `ax.errorbar()`: Plot mean ± 1σ with error bars
- `ax.axhspan()` / `ax.axhline()`: Reference band and line
- `Line2D()`: Custom legend handles
- `ax.legend()`: Configure legend display

**Dependencies:**
- `numpy`: Numerical operations (log scale calculations)
- `pandas`: Data manipulation and aggregation
- `matplotlib`: Visualization
- `accrediff`: AccreDiff constants and utilities

**Error Handling:**
- File not found checks for meteorite reference data
- Column validation (m_e must exist in simulation files)
- NaN handling in statistical calculations (dropna())
- Zero-value filtering for log scale (vals[vals > 0])

**Performance Notes:**
- Typical execution time: < 30 seconds for 1000+ simulations
- Memory usage: ~100 MB for large simulation sets
- Visualization rendering: ~5 seconds per figure

**Troubleshooting:**

| Issue | Solution |
|-------|----------|
| FileNotFoundError: Meteorite file | Check path: `../Cosmochemistry/Tables/A2_Meteorites_composition.csv` |
| KeyError: 'm_e' column missing | Ensure CSV files contain `m_e` column (from A_bulk_comparison.ipynb) |
| Empty Earth/Mars populations | Verify classification thresholds: Earth [0.7, 1.3], Mars [0.05, 0.3] |
| NaN in statistics | Check for zero or negative values in composition data |
| Legend overlaps plot | Adjust `bbox_to_anchor` parameter (currently =(1.02, 0.5)) |

**Publication Quality:**
- ✓ High-resolution DPI (300 for export)
- ✓ Nature-style formatting (serif font, minimal spines)
- ✓ Colorblind-friendly palette (Blue, Red, Brown, Orange, Green, Purple)
- ✓ Clear error representation (geometric mean ± 1σ on log scale)
- ✓ Reference standards explicitly marked
- ✓ Meteorite types clearly distinguished

**Last Updated:** 2026-04-13  
**Version:** 1.0  
**Status:** Ready for use