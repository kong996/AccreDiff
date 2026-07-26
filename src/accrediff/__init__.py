# accrediff/__init__.py

__version__ = "1.0.0"
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

# --- Collision outcomes (EDACM-style classifier) ---
from .collision_outcomes import (
    CollisionOutcome,
    EDACMConfig,
    add_collision_outcomes,
    classify_edacm_outcomes,
    impact_parameter_from_angle,
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
    IWCompositionCalculator_v2,
)
from .differentiation_sulfur import (
    MA2001_FE_S_NI_EPSILON_1873K,
    MA2001_FE_S_NI_GAMMA0_1873K,
    MA2001_EPSILON_REFERENCE_T_K,
    MA2001_EPSILON_TEMPERATURE_SCALING,
    GAMMA_FEO_SILICATE,
    MODEL_ASSUMPTIONS,
    SulfurSolverError,
    PhysicalStateError,
    AmbiguousPhysicalRootError,
    NiPartitionNotReachableError,
    SulfurPartitionNotReachableError,
    OxygenKDNotReachableError,
    IWTargetNotReachableError,
    SulfurResidualReport,
    audit_physical_state,
    Ma2001TernaryActivity,
    SulfurKDCalculator,
    SulfurKDParams,
    SulfurKDResult,
    SUPPLEMENT_MAJOR_KD_PARAMS,
    supplement_major_kd,
    D_S_model_molar,
    D_S_model_rose_weston_molar,
    ForwardKDOSolverSulfur,
    SulfurOLSolver,
    make_sulfur_solver_factory,
    sulfur_model_for_event,
    sulfur_params_from_composition,
    SulfurDifferentiationEventDriver,
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
                    

# --- Gas disk migration (Cresswell & Nelson 2008) ---
from .gas_migration import (
    Gas_ModelConfig,
    Gas_UnitSystem,
    Gas_DiskModel,
    Gas_TorqueModel,
    #Gas_TorqueModel_V2,
    Gas_MigrationIntegrator,
    Gas_build_model,
    #Gas_build_model_v2,
)
from .plotting import (
    plot_IW_heatmap,
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
    # collision outcomes
    "CollisionOutcome",
    "EDACMConfig",
    "add_collision_outcomes",
    "classify_edacm_outcomes",
    "impact_parameter_from_angle",
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
    "IWCompositionCalculator_v2",
    "MA2001_FE_S_NI_EPSILON_1873K",
    "MA2001_FE_S_NI_GAMMA0_1873K",
    "MA2001_EPSILON_REFERENCE_T_K",
    "MA2001_EPSILON_TEMPERATURE_SCALING",
    "GAMMA_FEO_SILICATE",
    "MODEL_ASSUMPTIONS",
    "SulfurSolverError",
    "PhysicalStateError",
    "AmbiguousPhysicalRootError",
    "NiPartitionNotReachableError",
    "SulfurPartitionNotReachableError",
    "OxygenKDNotReachableError",
    "IWTargetNotReachableError",
    "SulfurResidualReport",
    "audit_physical_state",
    "Ma2001TernaryActivity",
    "SulfurKDCalculator",
    "SulfurKDParams",
    "SulfurKDResult",
    "SUPPLEMENT_MAJOR_KD_PARAMS",
    "supplement_major_kd",
    "D_S_model_molar",
    "D_S_model_rose_weston_molar",
    "ForwardKDOSolverSulfur",
    "SulfurOLSolver",
    "make_sulfur_solver_factory",
    "sulfur_model_for_event",
    "sulfur_params_from_composition",
    "SulfurDifferentiationEventDriver",
    # gas migration
    "Gas_ModelConfig",
    "Gas_UnitSystem",
    "Gas_DiskModel",
    "Gas_TorqueModel",
    #"Gas_TorqueModel_V2",
    "Gas_MigrationIntegrator",
    "Gas_build_model",
    #"Gas_build_model_v2",
    # melt scaling
    "MeltScalingModel",
    # utils / plotting
    "normalize_max",
    "power_law",
    "error_rate",
    "WeightedECDF",
    "compute_error_metrics",
    #"generate_migration_map_v2",
    #"plot_migration_map_v2",
    # "images_to_video",
    # plotting
    "plot_IW_heatmap",
]
