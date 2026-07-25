# EDACM Collision Outcome Example

This folder contains an example notebook for applying the optional EDACM-style
collision-outcome classifier to real N-body collision records, built on top of
the `accrediff` package.

## Notebook

### `EDACM_Collision.ipynb` — Real N-body EDACM Outcomes

Reads a GENGA-style collision file from the bundled `note_02` N-body dataset,
preprocesses the stored collision-state vectors, and appends EDACM
collision-outcome diagnostics using `ad.add_collision_outcomes()`.

**Key steps:**
- Load `Collisionsnote_02.dat` from `examples/dataset/N_body_dataset/note_02/`
  using the same column convention as the accretion notebooks
- Filter to terrestrial-body collisions with `indexi >= 2` and `indexj >= 2`
- Convert GENGA-like quantities to SI units:
  - mass: solar masses to kg
  - length and radius: AU to m
  - velocity: code units to m s^-1
- Compute relative impact velocity from the two velocity vectors
- Derive impact angle and impact parameter from relative position and velocity:
  `cos(theta) = |dr . dv| / (|dr| |dv|)` and `b = sin(theta)`
- Run the EDACM classifier to assign outcome labels such as `perfect_merger`,
  `hit_and_run`, `partial_accretion`, `erosion`, and `super_catastrophic`
- Summarize the outcome distribution and visualize outcomes in
  impact-velocity versus impact-angle space

**Input:**
- `../../dataset/N_body_dataset/note_02/Collisionsnote_02.dat`

**Main outputs in memory:**
- `df_collision_edacm`: collision table converted to EDACM-ready SI quantities
- `df_real_outcomes`: collision table enriched with EDACM diagnostics
- `outcome_summary`: count and fraction of each EDACM outcome class

**Diagnostic columns added by `ad.add_collision_outcomes()`:**
- `outcome_id`
- `outcome`
- `grazing`
- `catastrophic`
- `alpha_interact`
- `q_r`
- `q_rd_star`
- `v_threshold`
- `m_largest_remnant_fraction`
- `m_largest_remnant`

---

## Workflow Note

This notebook demonstrates the preprocessing and diagnostic-enrichment layer
only. It does not yet map non-perfect-merger outcomes to AccreDiff mantle/core
material budgets. That coupling should be implemented separately because it
affects mass conservation, reaccretion assumptions, and event-by-event
metal-silicate differentiation.

