# accrediff/__init__.py

__version__ = "0.1.0"
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
from .utils import normalize_max, power_law, error_rate, WeightedECDF
#from .plotting import images_to_video

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
    # "images_to_video",
]