# accrediff/__init__.py

__version__ = "0.2.0"
__author__ = "Zhihui Kong"
__email__ = "kongzh0508@163.com"

# Keep constants as a submodule to avoid polluting the top-level namespace
from . import constant as constants

# --- Chemistry ---
from .chemistry import (
    MolarMassCalculator,
    KDCalculator,
    compute_core_ratio,
    compute_IW,
    compute_IW_mol,
)

# --- Accretion / impact processing ---
from .accretion import (
    CollisionTracer,
    GrowthModel,
    GrowthFitter,
    ImpactEventProcessor,
    resolve_product_id,
    build_pl_source_dict,
    update_partial_melting_process,
)

# --- Differentiation ---
from .differentiation import (
    KD_Params,
    KD_Result,
    Early_pressure,
    CMB_pressure,
    equil_pressure,
    P_to_T,
    EarlyComUpdater,
    ForwardKDOSolver,
    OLSolver,
)

# --- Melt scaling (Nakajima et al. 2021 style) ---
from .melt_model import MeltScalingModel

# --- Utilities / plotting ---
from .utils import (
    normalize_max, 
    power_law, 
    error_rate,
    compute_error_metrics, 
    WeightedECDF)
                    
#from .plotting import images_to_video
# --- Gas disk migration (Cresswell & Nelson 2008) ---
from .gas_migration import (
    Gas_ModelConfig,
    Gas_UnitSystem,
    Gas_DiskModel,
    Gas_TorqueModel,
    Gas_MigrationIntegrator,
    Gas_build_model,
)
__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "constants",
    # chemistry
    "MolarMassCalculator",
    "KDCalculator",
    "compute_core_ratio",
    "compute_IW",
    "compute_IW_mol",
    # accretion
    "CollisionTracer",
    "GrowthModel",
    "GrowthFitter",
    "ImpactEventProcessor",
    "resolve_product_id",
    "build_pl_source_dict",
    "update_partial_melting_process",
    # differentiation
    "KD_Params",
    "KD_Result",
    "Early_pressure",
    "CMB_pressure",
    "equil_pressure",
    "P_to_T",
    "EarlyComUpdater",
    "ForwardKDOSolver",
    "OLSolver",
    # melt scaling
    "MeltScalingModel",
    # utils / plotting
    "normalize_max",
    "power_law",
    "error_rate",
    "WeightedECDF",
    "compute_error_metrics",
    # "images_to_video",
]