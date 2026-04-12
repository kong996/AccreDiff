# Impact Event Analysis Examples

This folder contains example notebooks for analyzing impact events during planetary
accretion by classifying collision-induced melting regimes and tracking magma pool
evolution, built on top of the `accrediff` package.

## Notebooks

### 1. `A_Small_Impact.ipynb` — Small Impact Event Detection

Establishes thermal and barometric conditions for impact events during planetary
accretion by performing binary classification of collision-induced melting regimes.

**Key steps:**
- Load N-body collision history with impact records and thermal state data
- Initialize `ImpactEventProcessor` with configuration parameters:
  - Small impact pressure coefficient: `p_ratio = 0.6` (relative to P_CMB)
  - Global melting recovery duration: `duration = 5 Myr`
  - Impact magnitude ratios: `ratio = [0.1, 0.5]`
- Classify impact events into two thermal regimes via collision energy analysis:
  - **Small impacts** (P < 0.6 × P_CMB): Localized melting, no system-wide effects
  - **Global melting** (P ≥ 0.6 × P_CMB): Magma ocean formation, planet-wide thermal anomaly
- Update global melting process evolution tracking thermal relaxation over ~5 Myr
- Extract equilibrium pressure (P_equil) and planetary mass at each impact moment
- Visualize mass and pressure evolution trajectories color-coded by event type
  - Black dots: Small impacts
  - Orange dots: Global melting events
  - Grey line: Continuous pressure/mass evolution
- Export binary-classified impact events to Excel

**Output:**
- `Tables/{Model}_small_impact.xlsx` (one sheet per planet)
  - Columns: `Time`, `Mass`, `P_equil`, `target_id`, `events` ('small'/'global'), ...

---

### 2. `B_Giant_Impact.ipynb` — Partial Melting Pool Refinement

Refines impact event classification from binary to ternary by identifying the
intermediate partial melting regime and tracking localized melt pool dynamics.

**Key steps:**
- Load pre-classified impact events from `A_Small_Impact.ipynb` output
- Apply `update_partial_melting_process()` algorithm to identify intermediate-energy
  impacts that generate transient partial melting pools (0 < f < 1)
- Distinguish three thermal regimes via pressure thresholds:
  - **Small impacts** (P < 0.6 × P_CMB): No melting
  - **Partial melting** (0.6 × P_CMB ≤ P < P_CMB): Localized melt pools with mixed
    solid-liquid state, cooling timescale ~5 Myr
  - **Global melting** (P ≥ P_CMB): Complete magma ocean, system-wide liquidation
- Update event classification column with ternary categories: 'small'/'partial'/'global'
- Extract partial melting pool characteristics (pressure range, temporal clustering,
  mass distribution)
- Visualize three-tier event hierarchy with logarithmic time scale for early-time
  resolution:
  - Black dots: Small impacts
  - Red stars (★): Partial melting events (distinctive marker for emphasis)
  - Orange dots: Global melting events
  - Grey line: Continuous pressure/mass evolution
- Export refined ternary-classified impact events to Excel with enhanced metadata

**Output:**
- `Tables/{Model}_giant_impact.xlsx` (one sheet per planet)
  - Columns: `Time`, `Mass`, `P_equil`, `target_id`, `events` ('small'/'partial'/'global'), ...
  - Partial melting pool statistics: time range, pressure characterization, mass inventory

---