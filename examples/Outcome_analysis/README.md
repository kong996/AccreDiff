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

**Last Updated:** 2026-04-12  
**Version:** 1.0  
**Status:** Ready for use