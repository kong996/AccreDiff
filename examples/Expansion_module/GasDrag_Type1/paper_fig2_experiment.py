"""Named gas-disk experiments and saved-result matching for the Fig. 2 notebook."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from paper_fig2_disk import GasMigrationDiskAdapter, PaperFig2Disk

EXPERIMENTS = ("paper_reference", "gas_migration_matched", "gas_profile_baseline")


def make_disk(name: str):
    """Construct a fresh disk for one documented experiment."""
    if name == "paper_reference":
        return PaperFig2Disk()
    if name not in EXPERIMENTS:
        raise ValueError(f"unknown disk experiment {name!r}; choose from {EXPERIMENTS}")

    from accrediff.gas_migration import Gas_DiskModel, Gas_ModelConfig, Gas_UnitSystem

    common = {"Sigma0_cgs": 1700.0, "beta": 0.25, "tau_decay_myr": 1.0}
    if name == "gas_migration_matched":
        config = Gas_ModelConfig(a_alpha1=0.0, a_alpha2=1.5,
                                 def_h_1=0.025, **common)
    else:
        # The active gas_profile.ipynb baseline, not Gas_ModelConfig defaults.
        config = Gas_ModelConfig(a_alpha1=1.0, a_alpha2=0.5,
                                 def_h_1=0.0335743, alpha_visc=1e-4,
                                 prandtl=1.0, use_saturation=True, **common)
    units = Gas_UnitSystem()
    return GasMigrationDiskAdapter(Gas_DiskModel(config, units), units)


def output_path(directory: str | Path, name: str) -> Path:
    """Use distinct files so selecting a new disk never overwrites paper data."""
    if name not in EXPERIMENTS:
        raise ValueError(f"unknown disk experiment {name!r}")
    return Path(directory) / "outputs" / f"{name}_10myr.npz"


def matches_configuration(run, system, disk, gravity, drag, integration,
                          boundaries=None) -> bool:
    """Reject a saved run whose physical or numerical settings differ."""
    def normalized(value):
        return json.loads(json.dumps(value, sort_keys=True))

    if isinstance(disk, PaperFig2Disk):
        disk_description = {"type": type(disk).__name__, "config": asdict(disk)}
    elif isinstance(disk, GasMigrationDiskAdapter):
        disk_description = {"type": type(disk).__name__, "config": asdict(disk.disk.cfg)}
    else:
        return False
    expected = {
        "system": asdict(system),
        "integration": asdict(integration),
        "disk": disk_description,
        "gravity": asdict(gravity.config) if gravity else None,
        "boundaries": asdict(boundaries) if boundaries else None,
        "drag": ({"body": asdict(drag.planetesimal),
                  "flow": asdict(drag.flow), "law": drag.law,
                  "drag_coefficient": drag.cd} if drag else None),
    }
    return all(normalized(run.metadata.get(key)) == normalized(value)
               for key, value in expected.items())
