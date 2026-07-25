# Gas Disk Migration Example

This folder contains an example notebook for inspecting and validating a
modified GENGA gas-disk prescription. The CUDA source file is provided as
`../../dataset/gas_config/gas.cu`, and the corresponding Python-side validation
utilities are implemented in `src/accrediff/gas_migration.py`.

## Notebook

### `gas_profile.ipynb` — Gas Disk Profile and Migration

Builds the Python-side version of the modified gas-disk prescription, inspects
the surface-density and torque profiles, and compares analytic migration tracks
with bundled GENGA-style simulation outputs generated using the same gas setup.

**Key steps:**
- Configure a gas disk and migrating particle with `ad.Gas_ModelConfig`
- Build the disk, torque, unit-system, and integrator objects with
  `ad.Gas_build_model`
- Document the connection between the Python validation module and the bundled
  GENGA replacement source file `gas.cu`
- Plot gas surface-density profiles at selected epochs
- Evaluate the static torque profile and identify torque zero-crossings
- Integrate migration tracks for particles starting at 0.5, 1.0, and 1.5 AU
- Compare analytic migration tracks with the bundled simulation files
  `One_01.csv`, `One_02.csv`, and `One_03.csv`

**Input:**
- `../../dataset/gas_config/gas.cu` (modified GENGA gas-disk source file)
- `../../dataset/gas_config/One_01.csv`
- `../../dataset/gas_config/One_02.csv`
- `../../dataset/gas_config/One_03.csv`

**Main outputs in memory:**
- `cfg`: gas-disk and particle configuration
- `model`: dictionary returned by `ad.Gas_build_model`
- `profile_t0`: radial torque and disk-property profile at `t=0`
- `theory_tracks`: analytic migration tracks from the gas module
- `simulation_tracks`: bundled GENGA-style comparison tracks

**Main AccreDiff APIs used:**
- `ad.Gas_ModelConfig`
- `ad.Gas_build_model`
- `Gas_DiskModel.sigma_cgs`
- `Gas_TorqueModel.profile`
- `Gas_MigrationIntegrator.integrate_full`

---

## Workflow Note

The bundled `gas.cu` file can be used to replace the corresponding gas-disk
configuration source file in GENGA. This notebook does not compile or rerun
GENGA; instead, it provides a Python-side profile, torque, and migration
validation workflow for the same prescription. In an extended workflow, the gas
module can be used to design gas-disk prescriptions, validate gas-force
configurations, or compare analytic migration tracks with GENGA outputs.
