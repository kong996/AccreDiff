# Differentiation Examples

## Overview

This directory contains Jupyter notebooks demonstrating the **multi-stage core formation** process during planetary accretion using the AccreDiff framework.

## Contents

### A_Major_Differentiation.ipynb

Simulates core-mantle differentiation through multiple collision events and produces final planetary compositions.

**Key Features:**
- Multi-collision event simulation
- KD-based core-mantle equilibration
- Support for three collision types: small, global, and partial impacts
- CI-normalized compositional analysis

**Outputs:**
- X/MgO ratio plot (composition normalized to MgO)
- CI-normalized spider diagram (mantle vs. core phases)
- Excel file with time-evolving collision event data

---

## Input Data Requirements

### 1. Early Composition Data
**File:** `{Model}_early_composition.csv`
- Format: CSV with embryo IDs as index
- Required columns: `i` (ID), `m_e` (mass), major element oxides and pure metals
- Elements: MgO, Al2O3, CaO, FeO, NiO, SiO2, Fe, Ni, Si, O

### 2. Collision Event Data
**File:** `{Model}_giant_impact.xlsx`
- Format: Multi-sheet Excel file
- Required columns per sheet:
  - `Time`: Collision timestamp
  - `target_id`: Target embryo identifier
  - `impactor_id`: Impactor embryo identifier
  - `events`: Collision type ('small', 'global', 'partial')
  - `P_equil`: Equilibration pressure (GPa)
  - `T_equil`: Equilibration temperature (K)
  - `k_mantle`: Mantle re-equilibration fraction (for global/partial impacts)

---

## Collision Event Types

| Type | Target Mantle | Impactor Core | Description |
|------|---------------|---------------|-------------|
| **Small** | No re-eq | Sinks (f_core=0.7) | Small impactor; minimal mantle mixing |
| **Global** | Partial (k_mantle×f_mantle) | Partial (f_core=0.7) | Massive impact; significant mantle re-equilibration |
| **Partial** | Partial (k_mantle) | Partial (f_core=0.3) | Intermediate impact; moderate re-equilibration |

---

## Output Files

All results are saved to `./Tables/`:

### Main Output
- **D1_C2C1B1A2_core_formation.xlsx**
  - Multi-sheet workbook (one sheet per collision scenario)
  - Each sheet contains collision-by-collision composition evolution
  - Final row represents the planet's composition at end of scenario

### Figures (displayed during execution)
1. **X/MgO Ratio Plot**: Element abundances relative to MgO
2. **CI-normalized Spider Diagram**: Comparison with chondritic reference

---