# Accretion Examples

This folder contains example notebooks for reconstructing planetary accretion
histories and bulk compositions from N-body simulation outputs, built on top of
the `accrediff` package.

## Notebooks

### 1. `A_Embryos_composition.ipynb` — Early Embryo Compositions (Gas Disk Phase)

Reconstructs bulk compositions of planetary embryos at the end of the gas disk
phase (t ≤ 5 Myr) by tracing collision histories and mapping source materials.

**Key steps:**
- Load N-body snapshots and classify planetesimals by orbital zone (EF/EC/OC/CI)
  based on semi-major axis boundaries (1.0, 1.3, 3.0 AU)
- Partition particles into collision vs. original populations at t = 5 Myr
- Trace gas-disk collision histories using `CollisionTracer` to compute mixture
  fractions for merged particles
- Load meteorite endmember compositions and update elemental abundances via
  `EarlyComUpdater`
- Compute pressure–temperature conditions and identify embryos (m ≥ 0.1 M⊕)
- Plot linear X/MgO diagram and CI-normalized spider diagram (mantle ■ vs. core ●)
- Export all particle compositions to CSV

**Output:**
- `Tables/{model}_early_composition.csv`

---

### 2. `B_Collision_Track.ipynb` — Collision History Tracking

Reconstructs and visualizes the full collision/merger history for all final
planets (m ≥ 0.05 M⊕) across the entire simulation timespan.

**Key steps:**
- Load N-body snapshots and identify final planets at simulation end (t = 100 Myr)
- Track mass and semi-major axis evolution across all time steps
- Filter post-gas-disk collisions (t > 5 Myr) and resolve product particle IDs
  using `resolve_product_id`
- Trace complete merger trees for each planet via `CollisionTracer.trace_full_history()`
- Visualize mass evolution curves (0–100 Myr) and accretion event timelines (log-scale)
- Export collision histories to Excel (one sheet per planet)

**Output:**
- `Tables/{model}_collision_history.xlsx`

---

### 3. `C_Planetary_Source.ipynb` — Planetary Source Material Analysis

Determines fractional contributions from each chondritic reservoir (EF/EC/OC/CI)
to final Earth-like and Mars-like planets.

**Key steps:**
- Identify Earth-like (0.7–1.3 M⊕) and Mars-like (0.05–0.3 M⊕) planets
- Trace collision histories and map source materials back to initial orbital
  positions using `build_pl_source_dict`
- Build weighted ECDFs (`WeightedECDF`) and evaluate cumulative mass fractions
  at boundary orbits (1.0, 1.3, 3.0 AU)
- Convert cumulative fractions to absolute chondrite compositions:
  EF = F(1.0), EC = F(1.3)−F(1.0), OC = F(3.0)−F(1.3), CI = 1−F(3.0)
- Visualize cumulative source distributions and composition pie charts
- Export planetary source inventory to CSV

**Output:**
- `Tables/{model}_source.csv`

---

### 4. `D_Planetary_Evolution.ipynb` — Planetary Growth Evolution

Reconstructs temporal evolution of planetary accretion with real-time composition
tracking using mass-weighted conservative mixing.

**Key steps:**
- Build initial composition dictionary (particle ID → chondrite fractions) and
  update dynamically at each collision via mass-weighted averaging
- Trace full collision histories with composition tracking for each planet
- Fit exponential growth models to mass evolution curves using `GrowthFitter`
  and `GrowthModel`
- Compute key accretion timescales (40%, 60%, 80% of final mass) via linear
  interpolation
- Visualize mass growth curves with fitted models and stacked area plots of
  composition evolution over time
- Export growth timescale statistics to CSV

**Output:**
- `Tables/{model}_evolution_timescales.csv`