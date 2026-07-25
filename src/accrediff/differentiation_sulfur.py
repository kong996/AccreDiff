"""Independent Fe-Ni-S metal-silicate differentiation module.

This module does not inherit numerical or data classes from
``accrediff.differentiation``.  It intentionally retains the familiar field
layout (``Fe_t`` ... ``n`` on input and ``x_`` ... ``residual`` on output) so
existing accretion workflows can pass and consume compositions in the same
shape without coupling the two implementations.

The sulfur calculation contains two explicitly staged models:

* ``rose_weston`` is the primordial-differentiation (Brennan Eq. S8) model.
* ``boujibar`` is the composition-dependent model used in later events
  (Brennan Eq. S7; Boujibar et al. 2014 Eqs. 6 and 11).

All public compositions and returned sulfur partition coefficients are on a
molar basis.  Empirical mass-concentration coefficients are converted inside
the sulfur model functions.  Sulfur-sensitive exchange coefficients and
Fe-S-Ni activities are supplied by :class:`SulfurKDCalculator`, whose
Ma (2001) implementation is reproduced from the project reference notebook.

Current scientific scope
------------------------
* The Ma interaction parameters tabulated at 1873 K are used directly at all
  calculation temperatures.  No epsilon-temperature scaling is applied.
* ``gamma_FeO_silicate = 1.7`` is fixed following Brennan et al. supplementary
  Methods S1; ``gamma_Fe_metal`` is recalculated with Ma (2001).
* The sulfur correction to oxygen partitioning remains deliberately deferred.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from typing import Callable, Dict, Literal, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import toms748

from .chemistry import MolarMassCalculator
from .constant import Elements

SulfurModel = Literal["rose_weston", "boujibar"]

# Atomic masses come exclusively from ``constant.Elements``.  Component
# masses are derived once with the project's shared calculator rather than
# maintained as a second independent table.  ``AlO1.5`` is one half formula
# unit of Al2O3; the calculator intentionally accepts integer stoichiometry.
_MOLAR_MASS_CALCULATOR = MolarMassCalculator()
_COMPONENT_MOLAR_MASS = {
    "FeO": _MOLAR_MASS_CALCULATOR.molar_mass("FeO"),
    "NiO": _MOLAR_MASS_CALCULATOR.molar_mass("NiO"),
    "SiO2": _MOLAR_MASS_CALCULATOR.molar_mass("SiO2"),
    "MgO": _MOLAR_MASS_CALCULATOR.molar_mass("MgO"),
    "AlO1.5": 0.5 * _MOLAR_MASS_CALCULATOR.molar_mass("Al2O3"),
    "CaO": _MOLAR_MASS_CALCULATOR.molar_mass("CaO"),
}

# Brennan et al. supplementary Table S1: Ni is from Siebert et al. (2011),
# whereas Si and O are from Fischer et al. (2015).  This module deliberately
# uses this single Ni/Si/O parameter set at all pressures;
# the alternative P < 5 GPa branches in the general KDCalculator are outside
# the currently approved model scope.
SUPPLEMENT_MAJOR_KD_PARAMS = {
    "Ni": (0.35, 2934.0, -83.0),
    "Si": (1.3, -13500.0, 0.0),
    "O": (0.6, -3800.0, 22.0),
}

# Fe-S-Ni subset of ``epsilon_1873K.json`` used by
# ``GMD2026_v1/activity_coefficient/ternary_Ma2001.ipynb``.  The notebook
# identifies the dataset as Jennings (2021) at 1873 K.  The Ma equations use
# these dimensionless interaction parameters directly; no unverified
# temperature rescaling is applied here.
MA2001_FE_S_NI_EPSILON_1873K = {
    "S": {"S": -5.66, "Ni": 2.17},
    "Ni": {"S": 2.17, "Ni": 0.119},
}
MA2001_FE_S_NI_GAMMA0_1873K = {
    "S": 1.0,
    "Ni": 1.0,
}

# Model-scope records.  ``False`` is descriptive rather than a switch: no
# alternative scaling law is implemented until one is explicitly approved.
MA2001_EPSILON_REFERENCE_T_K = 1873.0
MA2001_EPSILON_TEMPERATURE_SCALING = False

# Brennan et al. (2020) supplementary Methods S1, p. 2:
# gamma_FeO^silicate was taken as 1.7 (Holzheid et al., 1996;
# O'Neill et al., 2002).  This is intentionally a module constant rather than
# a solver option in the current model.
GAMMA_FEO_SILICATE = 1.7

MODEL_ASSUMPTIONS = (
    "Ma epsilon values fixed at their 1873 K reference values",
    "gamma_FeO_silicate fixed at 1.7",
    "sulfur correction to oxygen partitioning omitted",
    "TiO2, Na2O, and K2O terms omitted from sulfide capacity",
    "single high-pressure Si/O parameter set used at all pressures",
    "Ma activities use a true Fe-Ni-S ternary normalization",
    "IW uses gamma_Fe times ternary-normalized X_Fe",
    "Si and O are excluded from the Fe-Ni-S activity subsystem",
    "published empirical relations may be extrapolated outside calibration ranges",
)


class SulfurSolverError(RuntimeError):
    """Base class for strict sulfur-solver failures."""


class PhysicalStateError(SulfurSolverError):
    """A trial composition violates conservation or phase admissibility."""


class AmbiguousPhysicalRootError(SulfurSolverError):
    """More than one physical algebraic branch exists without a branch hint."""


class NiPartitionNotReachableError(SulfurSolverError):
    """No physical Ni-partition root was found within the search budget."""


class SulfurPartitionNotReachableError(SulfurSolverError):
    """No physical sulfur-partition root was found within the search budget."""


class OxygenKDNotReachableError(SulfurSolverError):
    """No matching FeO state was found within the search budget."""


class IWTargetNotReachableError(SulfurSolverError):
    """No matching oxygen-loss state was found within the search budget."""


@dataclass(frozen=True)
class SulfurResidualReport:
    """Named conservation and equilibrium residuals for one phase state."""

    residuals: Dict[str, float]
    nonnegative: bool
    phase_totals_positive: bool
    all_finite: bool

    @property
    def physical(self) -> bool:
        return self.nonnegative and self.phase_totals_positive and self.all_finite

    @property
    def max_abs(self) -> float:
        finite = [
            abs(value)
            for value in self.residuals.values()
            if math.isfinite(value)
        ]
        return max(finite, default=0.0)


def _normalized_difference(
    left: float,
    right: float,
    floor: float = np.finfo(float).tiny,
) -> float:
    """Return a symmetric dimensionless residual for two non-negative terms."""

    return (left - right) / max(abs(left) + abs(right), floor)


def _mapping_from_state(state: object) -> Mapping[str, object]:
    if isinstance(state, Mapping):
        return state
    if hasattr(state, "__dict__"):
        return vars(state)
    raise TypeError("state must be a mapping or dataclass-like result.")


def audit_physical_state(
    state: object,
    params: "SulfurKDParams",
    *,
    KD_O_target: Optional[float] = None,
    IW_target: Optional[float] = None,
    IW_model: Optional[float] = None,
    physical_tol: float = 1e-12,
) -> SulfurResidualReport:
    """Independently audit conservation, partition equations, and positivity.

    The function deliberately reconstructs every equation from the returned
    amounts.  It does not trust any solver convergence flag.
    """

    values = _mapping_from_state(state)

    def amount(name: str) -> float:
        return float(values.get(name, 0.0))

    x_ = amount("x_")
    a_ = amount("a_")
    y_ = amount("y_")
    b_ = amount("b_")
    z_ = amount("z_")
    c_ = amount("c_")
    d_ = amount("d_")
    s_met = amount("S_met_")
    s_sil = amount("S_sil_")
    kd_ni = amount("KD_Ni_effective")
    d_s_mole = amount("D_S_mole")
    kd_o = amount("KD_O_")

    amounts = (x_, a_, y_, b_, z_, c_, d_, s_met, s_sil)
    inventory_scale = max(
        params.Fe_t
        + params.Ni_t
        + params.Si_t
        + params.u
        + params.m
        + params.n
        + params.S_t,
        np.finfo(float).tiny,
    )
    amount_tol = physical_tol * inventory_scale
    all_finite = all(math.isfinite(value) for value in amounts)
    nonnegative = all(value >= -amount_tol for value in amounts)

    metal_sum = a_ + b_ + c_ + d_ + s_met
    sil_sum = x_ + y_ + z_ + params.u + params.m + params.n + s_sil
    phase_totals_positive = (
        math.isfinite(metal_sum)
        and math.isfinite(sil_sum)
        and metal_sum > amount_tol
        and sil_sum > amount_tol
    )
    oxygen_total = params.Fe_t + params.Ni_t + 2.0 * params.Si_t - params.O_L

    residuals: Dict[str, float] = {
        "Fe_balance": (x_ + a_ - params.Fe_t)
        / max(params.Fe_t, amount_tol, np.finfo(float).tiny),
        "Ni_balance": (y_ + b_ - params.Ni_t)
        / max(params.Ni_t, amount_tol, np.finfo(float).tiny),
        "Si_balance": (z_ + c_ - params.Si_t)
        / max(params.Si_t, amount_tol, np.finfo(float).tiny),
        "O_balance": (
            x_ + y_ + 2.0 * z_ + d_ - oxygen_total
        )
        / max(abs(oxygen_total), amount_tol, np.finfo(float).tiny),
        "S_balance": (
            s_met + s_sil - params.S_t
        )
        / max(params.S_t, amount_tol, np.finfo(float).tiny),
    }

    if phase_totals_positive:
        residuals["Ni_equilibrium"] = _normalized_difference(
            b_ * x_,
            kd_ni * y_ * a_,
        )
        residuals["Si_equilibrium"] = _normalized_difference(
            c_ * x_**2 * metal_sum,
            params.KD_Si * z_ * a_**2 * sil_sum,
        )
        if params.S_t > 0.0:
            residuals["S_equilibrium"] = _normalized_difference(
                s_met * sil_sum,
                d_s_mole * s_sil * metal_sum,
            )
        else:
            residuals["S_equilibrium"] = 0.0
    else:
        residuals.update(
            Ni_equilibrium=math.inf,
            Si_equilibrium=math.inf,
            S_equilibrium=math.inf,
        )

    if KD_O_target is not None:
        if (
            math.isfinite(KD_O_target)
            and KD_O_target > 0.0
            and math.isfinite(kd_o)
            and kd_o > 0.0
        ):
            residuals["O_KD_log"] = math.log10(kd_o / KD_O_target)
        else:
            residuals["O_KD_log"] = math.inf
    if IW_target is not None:
        residuals["IW"] = (
            float(IW_model) - IW_target
            if IW_model is not None and math.isfinite(IW_model)
            else math.inf
        )

    all_finite = all_finite and all(
        math.isfinite(value) for value in residuals.values()
    )
    return SulfurResidualReport(
        residuals=residuals,
        nonnegative=nonnegative,
        phase_totals_positive=phase_totals_positive,
        all_finite=all_finite,
    )


def supplement_major_kd(element: Literal["Si", "O"], P: float, T: float) -> float:
    """Return the approved Brennan Table S1 Si or O exchange coefficient."""

    if not math.isfinite(P) or P < 0.0:
        raise ValueError("Pressure must be finite and non-negative.")
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("Temperature must be finite and positive.")
    a, b, c = SUPPLEMENT_MAJOR_KD_PARAMS[element]
    return 10.0 ** (a + b / T + c * P / T)


class Ma2001TernaryActivity:
    """Ma (2001) ternary activity-coefficient model (Eqs. 25 and 26).

    Components 1, 2 and 3 are respectively ``solvent``, ``solute2`` and
    ``solute3``.  The two supplied solute mole fractions must obey
    ``x2 >= 0``, ``x3 >= 0`` and ``x2 + x3 < 1``; the solvent fraction is
    calculated as ``1 - x2 - x3``.

    This is the tested, dependency-free form of the implementation in
    ``GMD2026_v1/activity_coefficient/ternary_Ma2001.ipynb``.  All logarithms
    are natural logarithms.  A missing cross-interaction parameter is treated
    as zero and a missing infinite-dilution activity coefficient as one, just
    as in that notebook.
    """

    def __init__(
        self,
        *,
        solvent: str,
        solute2: str,
        solute3: str,
        x2: float,
        x3: float,
        epsilon_db: Mapping[str, Mapping[str, float]],
        gamma0_db: Optional[Mapping[str, float]] = None,
    ) -> None:
        self.solvent = str(solvent)
        self.solute2 = str(solute2)
        self.solute3 = str(solute3)
        if len({self.solvent, self.solute2, self.solute3}) != 3:
            raise ValueError("The solvent and two solutes must be distinct.")
        self.x2 = float(x2)
        self.x3 = float(x3)
        if (
            not math.isfinite(self.x2)
            or not math.isfinite(self.x3)
            or self.x2 < 0.0
            or self.x3 < 0.0
            or self.x2 + self.x3 >= 1.0
        ):
            raise ValueError(
                "Ternary mole fractions require x2 >= 0, x3 >= 0, "
                "and x2 + x3 < 1."
            )
        self.x1 = 1.0 - self.x2 - self.x3
        self.epsilon_db = {
            str(i): {str(j): float(value) for j, value in row.items()}
            for i, row in epsilon_db.items()
        }
        self.gamma0_db = {
            str(component): float(value)
            for component, value in (gamma0_db or {}).items()
        }
        if any(
            not math.isfinite(value)
            for row in self.epsilon_db.values()
            for value in row.values()
        ):
            raise ValueError("Interaction parameters must be finite.")
        if any(
            not math.isfinite(value) or value <= 0.0
            for value in self.gamma0_db.values()
        ):
            raise ValueError("Infinite-dilution activity coefficients must be positive.")

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(solvent={self.solvent!r}, "
            f"solute2={self.solute2!r}, solute3={self.solute3!r}, "
            f"x1={self.x1!r}, x2={self.x2!r}, x3={self.x3!r})"
        )

    def _epsilon(self, first: str, second: str) -> float:
        """Return a symmetric interaction parameter, defaulting to zero."""

        first_row = self.epsilon_db.get(first, {})
        if second in first_row:
            return first_row[second]
        return self.epsilon_db.get(second, {}).get(first, 0.0)

    def _gamma0(self, component: str) -> float:
        return self.gamma0_db.get(component, 1.0)

    @staticmethod
    def _ln_one_minus_x_over_x(x: float) -> float:
        """Stable evaluation of ``ln(1-x)/x``, including its limit at zero."""

        if x == 0.0:
            return -1.0
        return math.log1p(-x) / x

    def ln_gamma_solvent(self) -> float:
        """Return ``ln(gamma_1)`` from Ma (2001) Eq. 25."""

        x2, x3 = self.x2, self.x3
        eps22 = self._epsilon(self.solute2, self.solute2)
        eps33 = self._epsilon(self.solute3, self.solute3)
        eps23 = self._epsilon(self.solute2, self.solute3)
        return (
            eps22 * (x2 + math.log1p(-x2))
            + eps33 * (x3 + math.log1p(-x3))
            + eps23
            * x2
            * x3
            * (1.0 - 1.0 / (1.0 - x2) - 1.0 / (1.0 - x3))
            - 0.5
            * eps23
            * x2**2
            * x3**2
            * (
                3.0 / (1.0 - x2)
                + 3.0 / (1.0 - x3)
                + x2 / (1.0 - x2) ** 2
                + x3 / (1.0 - x3) ** 2
                - 3.0
            )
        )

    def _ln_gamma_solute(
        self,
        solute: str,
        x_solute: float,
        other_solute: str,
        x_other: float,
    ) -> float:
        """Return ``ln(gamma_i)`` from Ma (2001) Eq. 26."""

        eps_ii = self._epsilon(solute, solute)
        eps_ki = self._epsilon(other_solute, solute)
        return (
            self.ln_gamma_solvent()
            + math.log(self._gamma0(solute))
            - eps_ii * math.log1p(-x_solute)
            - eps_ki
            * x_other
            * (
                1.0
                + self._ln_one_minus_x_over_x(x_other)
                - 1.0 / (1.0 - x_solute)
            )
            + eps_ki
            * x_other**2
            * x_solute
            * (
                1.0 / (1.0 - x_solute)
                + 1.0 / (1.0 - x_other)
                + x_solute / (2.0 * (1.0 - x_solute) ** 2)
                - 1.0
            )
        )

    def ln_gamma(self, component: str) -> float:
        """Return the natural-log activity coefficient of one component."""

        if component == self.solvent:
            return self.ln_gamma_solvent()
        if component == self.solute2:
            return self._ln_gamma_solute(
                self.solute2,
                self.x2,
                self.solute3,
                self.x3,
            )
        if component == self.solute3:
            return self._ln_gamma_solute(
                self.solute3,
                self.x3,
                self.solute2,
                self.x2,
            )
        raise ValueError(f"Component {component!r} is not in this ternary system.")

    def gamma(self, component: str) -> float:
        """Return the activity coefficient of one component."""

        return math.exp(self.ln_gamma(component))

    def activity(self, component: str) -> float:
        """Return ``a_i = gamma_i * X_i`` for one component."""

        composition = {
            self.solvent: self.x1,
            self.solute2: self.x2,
            self.solute3: self.x3,
        }
        if component not in composition:
            raise ValueError(f"Component {component!r} is not in this ternary system.")
        return self.gamma(component) * composition[component]

    def results(self) -> Dict[str, Dict[str, float]]:
        """Return the same composition/activity mapping as the source notebook."""

        components = (self.solvent, self.solute2, self.solute3)
        composition = {
            self.solvent: self.x1,
            self.solute2: self.x2,
            self.solute3: self.x3,
        }
        return {
            "composition": composition,
            "ln_gamma": {name: self.ln_gamma(name) for name in components},
            "gamma": {name: self.gamma(name) for name in components},
            "activity": {name: self.activity(name) for name in components},
        }


class SulfurKDCalculator:
    """Calculate base exchange coefficients and Fe-Ni-S activities.

    The class mirrors :class:`accrediff.chemistry.KDCalculator`: empirical
    ``(a, b, c)`` tuples are stored in ``params`` and evaluated as
    ``log10(KD_base) = a + b/T + c*P/T``.  :meth:`get_KD` deliberately
    returns only that composition-independent base value for Ni, Si, or O.

    The Ma (2001) activity model is exposed separately.  The Ni partition
    equation applies ``gamma_Fe / gamma_Ni`` to the externally stored base
    coefficient.  ``X_S`` and ``X_Ni`` passed to the activity methods must be
    normalized within the Fe-Ni-S subsystem; Fe is the remaining fraction.

    Temperature scope
    -----------------
    The epsilon values are the fixed 1873 K values used by
    ``ternary_Ma2001.ipynb``.  Per the current modeling decision, they are
    applied directly at every ``T``: epsilon-temperature scaling is not
    implemented. ``T`` affects only the empirical KD pressure-temperature
    term.  The public validation remains ``P >= 0`` and ``T > 0`` so callers
    may deliberately extrapolate beyond published calibration ranges.
    """

    def __init__(
        self,
        *,
        epsilon_db: Optional[Mapping[str, Mapping[str, float]]] = None,
        gamma0_db: Optional[Mapping[str, float]] = None,
    ) -> None:
        self.params: Dict[str, Tuple[float, float, float]] = dict(
            SUPPLEMENT_MAJOR_KD_PARAMS
        )
        epsilon_source = epsilon_db or MA2001_FE_S_NI_EPSILON_1873K
        gamma0_source = gamma0_db or MA2001_FE_S_NI_GAMMA0_1873K
        self.epsilon_db = {
            str(i): {str(j): float(value) for j, value in row.items()}
            for i, row in epsilon_source.items()
        }
        self.gamma0_db = {
            str(component): float(value)
            for component, value in gamma0_source.items()
        }

    def activity_model(
        self,
        *,
        X_Ni: float,
        X_S: float,
    ) -> Ma2001TernaryActivity:
        """Build the Fe-S-Ni activity model for a metal composition."""

        return Ma2001TernaryActivity(
            solvent="Fe",
            solute2="S",
            solute3="Ni",
            x2=X_S,
            x3=X_Ni,
            epsilon_db=self.epsilon_db,
            gamma0_db=self.gamma0_db,
        )

    def get_activity_coefficients(
        self,
        *,
        X_Ni: float,
        X_S: float,
    ) -> Dict[str, float]:
        """Return Ma-model ``gamma_Fe``, ``gamma_S`` and ``gamma_Ni``."""

        model = self.activity_model(X_Ni=X_Ni, X_S=X_S)
        return {
            component: model.gamma(component)
            for component in ("Fe", "S", "Ni")
        }

    def get_ni_activity_factor(
        self,
        *,
        X_Ni: float,
        X_S: float,
    ) -> float:
        """Return the Fe-Ni-S correction ``gamma_Fe / gamma_Ni``."""

        model = self.activity_model(X_Ni=X_Ni, X_S=X_S)
        return math.exp(model.ln_gamma("Fe") - model.ln_gamma("Ni"))

    def get_KD(
        self,
        element: str,
        P: float,
        T: float,
    ) -> float:
        """Return the composition-independent base KD for Ni, Si, or O."""

        if element not in self.params:
            available = ", ".join(sorted(self.params))
            raise ValueError(f"Unknown element {element!r}; available: {available}.")
        if not math.isfinite(P) or P < 0.0:
            raise ValueError("Pressure must be finite and non-negative.")
        if not math.isfinite(T) or T <= 0.0:
            raise ValueError("Temperature must be finite and positive.")

        a, b, c = self.params[element]
        return 10.0 ** (a + b / T + c * P / T)


@dataclass
class SulfurKDParams:
    """Standalone solver input with the same field layout as ``KD_Params``.

    All composition fields are molar amounts. ``KD_Ni`` and ``KD_Si`` are
    externally supplied base exchange coefficients. The solver applies the
    Fe-Ni-S activity factor to ``KD_Ni`` inside the Ni partition equation.
    """

    Fe_t: float
    Ni_t: float
    Si_t: float
    O_L: float
    KD_Ni: float
    KD_Si: float
    u: float = 0.0
    m: float = 0.0
    n: float = 0.0
    S_t: float = 0.0

    def __post_init__(self) -> None:
        composition_fields = (
            "Fe_t",
            "Ni_t",
            "Si_t",
            "u",
            "m",
            "n",
            "S_t",
        )
        for name in composition_fields:
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if not math.isfinite(self.O_L):
            raise ValueError("O_L must be finite.")
        oxygen_total = self.Fe_t + self.Ni_t + 2.0 * self.Si_t - self.O_L
        if oxygen_total < -1e-12 * max(
            self.Fe_t + self.Ni_t + 2.0 * self.Si_t,
            np.finfo(float).tiny,
        ):
            raise ValueError(
                "O_L implies a negative distributable oxygen inventory."
            )
        for name in ("KD_Ni", "KD_Si"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")


@dataclass
class SulfurKDResult:
    """Standalone result preserving the original ``KD_Result`` field order.

    ``sulfur_residual`` is the logarithmic equilibrium residual
    ``ln(X_S^metal / X_S^silicate) - ln(D_S_mole)``.  This definition remains
    well conditioned when either sulfur phase fraction is extremely small.
    """

    x_: float
    a_: float
    y_: float
    b_: float
    z_: float
    c_: float
    d_: float
    KD_O_: float
    residual: float
    S_met_: float = 0.0
    S_sil_: float = 0.0
    D_S_weight: float = 0.0
    D_S_mole: float = 0.0
    sulfur_residual: float = 0.0
    sulfur_iterations: int = 0
    sulfur_converged: bool = True
    sulfur_model: str = "none"
    coupled_residual: float = 0.0
    coupled_iterations: int = 0
    coupled_converged: bool = True
    KD_Ni_effective: float = 0.0
    KD_Ni_base: float = 0.0
    KD_Ni_activity_factor: float = 1.0
    gamma_Fe_metal: float = 1.0
    gamma_Ni_metal: float = 1.0
    gamma_S_metal: float = 1.0
    X_Fe_ternary: float = 1.0
    X_Ni_ternary: float = 0.0
    X_S_ternary: float = 0.0
    major_converged: bool = False
    oxygen_converged: bool = False
    iw_converged: bool = False
    physical_converged: bool = False
    overall_converged: bool = False
    failure_reason: Optional[str] = None
    residuals: Dict[str, float] = field(default_factory=dict)
    model_assumptions: Tuple[str, ...] = MODEL_ASSUMPTIONS


def _validate_phase_fractions(values: np.ndarray, phase: str) -> None:
    if np.any(~np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
        raise ValueError(f"{phase} mole fractions must be finite values in [0, 1].")
    if not np.isclose(float(values.sum()), 1.0, rtol=0.0, atol=1e-8):
        raise ValueError(f"{phase} mole fractions must sum to one.")


def _log_D_S_model_rose_weston_molar(
    P: float,
    T: float,
    *,
    X_Fe_met: float,
    X_Ni_met: float,
    X_O_met: float,
    X_Si_met: float,
    X_S_met: float,
    X_FeO_sil: float,
    X_NiO_sil: float,
    X_SiO2_sil: float,
    X_MgO_sil: float,
    X_AlO1_5_sil: float,
    X_CaO_sil: float,
    X_S_sil: float,
    X_C_met: float = 0.0,
) -> float:
    """Return ``ln(D_S_mole)`` for the Rose-Weston sulfur relation.

    Rose-Weston et al. (2009) Eq. 16, used as Eq. S8 in the Brennan
    supplementary calculation, gives a mass-concentration coefficient:

    ``D_S_weight = w_S_metal / w_S_silicate``.

    This function accepts whole-phase mole fractions, calculates both
    phase-average molar masses, and returns

    ``D_S_mole = D_S_weight * MW_metal / MW_silicate``.

    The restricted regression dataset covers P >= 2.8 GPa and approximately
    2233--2693 K, with a_FeO = 0.10 +/- 0.02 and most runs at fO2 >= IW-2.
    These limits are documentation only: any finite ``P >= 0`` and ``T > 0``
    remains accepted for deliberate extrapolation.
    """

    if not math.isfinite(P) or P < 0.0:
        raise ValueError("Pressure must be finite and non-negative.")
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("Temperature must be finite and positive.")
    metal_X = np.asarray(
        [X_Fe_met, X_Ni_met, X_O_met, X_Si_met, X_S_met, X_C_met],
        dtype=float,
    )
    silicate_X = np.asarray(
        [
            X_FeO_sil,
            X_NiO_sil,
            X_SiO2_sil,
            X_MgO_sil,
            X_AlO1_5_sil,
            X_CaO_sil,
            X_S_sil,
        ],
        dtype=float,
    )
    _validate_phase_fractions(metal_X, "Metal")
    _validate_phase_fractions(silicate_X, "Silicate")

    mw_metal = (
        X_Fe_met * Elements["Fe"]
        + X_Ni_met * Elements["Ni"]
        + X_O_met * Elements["O"]
        + X_Si_met * Elements["Si"]
        + X_S_met * Elements["S"]
        + X_C_met * Elements["C"]
    )
    mw_silicate = (
        X_FeO_sil * _COMPONENT_MOLAR_MASS["FeO"]
        + X_NiO_sil * _COMPONENT_MOLAR_MASS["NiO"]
        + X_SiO2_sil * _COMPONENT_MOLAR_MASS["SiO2"]
        + X_MgO_sil * _COMPONENT_MOLAR_MASS["MgO"]
        + X_AlO1_5_sil * _COMPONENT_MOLAR_MASS["AlO1.5"]
        + X_CaO_sil * _COMPONENT_MOLAR_MASS["CaO"]
        + X_S_sil * Elements["S"]
    )
    if (
        not math.isfinite(mw_metal)
        or not math.isfinite(mw_silicate)
        or mw_metal <= 0.0
        or mw_silicate <= 0.0
    ):
        raise ValueError("Phase-average molar masses must be finite and positive.")
    log_d_weight = math.log(10.0) * (
        -4.37 + 13686.0 / T + 217.49 * P / T
    )
    log_d_mole = log_d_weight + math.log(mw_metal) - math.log(mw_silicate)
    if not math.isfinite(log_d_mole):
        raise ValueError("Rose-Weston log(D_S) is not finite.")
    return log_d_mole


def D_S_model_rose_weston_molar(
    P: float,
    T: float,
    *,
    X_Fe_met: float,
    X_Ni_met: float,
    X_O_met: float,
    X_Si_met: float,
    X_S_met: float,
    X_FeO_sil: float,
    X_NiO_sil: float,
    X_SiO2_sil: float,
    X_MgO_sil: float,
    X_AlO1_5_sil: float,
    X_CaO_sil: float,
    X_S_sil: float,
    X_C_met: float = 0.0,
) -> float:
    """Return Rose-Weston ``D_S`` on a molar basis.

    Inputs are whole-phase mole fractions.  The empirical mass-concentration
    coefficient is converted internally to a molar coefficient.  See
    :func:`_log_D_S_model_rose_weston_molar` for the documented calibration
    range; extrapolation remains allowed.
    """

    log_d_mole = _log_D_S_model_rose_weston_molar(
        P,
        T,
        X_Fe_met=X_Fe_met,
        X_Ni_met=X_Ni_met,
        X_O_met=X_O_met,
        X_Si_met=X_Si_met,
        X_S_met=X_S_met,
        X_FeO_sil=X_FeO_sil,
        X_NiO_sil=X_NiO_sil,
        X_SiO2_sil=X_SiO2_sil,
        X_MgO_sil=X_MgO_sil,
        X_AlO1_5_sil=X_AlO1_5_sil,
        X_CaO_sil=X_CaO_sil,
        X_S_sil=X_S_sil,
        X_C_met=X_C_met,
    )
    try:
        d_mole = math.exp(log_d_mole)
    except OverflowError as exc:
        raise ValueError(
            "Rose-Weston D_S exceeds the representable floating-point range."
        ) from exc
    if not math.isfinite(d_mole) or d_mole <= 0.0:
        raise ValueError(
            "Rose-Weston D_S lies outside the positive floating-point range."
        )
    return d_mole


def _log_D_S_model_boujibar_molar(
    P: float,
    T: float,
    *,
    X_FeO_oxide: float,
    X_CaO_oxide: float,
    X_MgO_oxide: float,
    X_FeO_sil_phase: float,
    X_Fe_met: float,
    X_Si_met: float,
    X_Ni_met: float,
    X_O_met: float,
    MW_metal: float,
    MW_silicate: float,
    X_C_met: float = 0.0,
) -> float:
    """Return ``ln(D_S_mole)`` for Boujibar's sulfur relation.

    Inputs are mole fractions.  Internally:

    * Eq. 6 uses oxide-normalized mole fractions, with sulfur excluded;
    * Eq. 11 uses mass fractions, converted here from whole-phase mole
      fractions and phase-average molar masses;
    * Eq. 11's mass-concentration result is converted back to ``D_S_mole``.

    Model scope
    -----------
    The ``TiO2``, ``Na2O`` and ``K2O`` terms in the complete sulfide-capacity
    expression are deliberately omitted until those components are tracked.
    Carbon defaults to zero because the major-element solver does not track C.
    Boujibar et al. report calibration over approximately 1--25 GPa,
    1623--2733 K, and IW-6.6 to IW+0.4. These limits are not enforced; the
    function accepts finite ``P >= 0`` and ``T > 0`` for controlled
    extrapolation.
    """

    if not math.isfinite(P) or P < 0.0:
        raise ValueError("Pressure must be finite and non-negative.")
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("Temperature must be finite and positive.")
    if (
        not math.isfinite(MW_metal)
        or not math.isfinite(MW_silicate)
        or MW_metal <= 0.0
        or MW_silicate <= 0.0
    ):
        raise ValueError("Phase-average molar masses must be finite and positive.")

    oxide_values = np.asarray(
        [X_FeO_oxide, X_CaO_oxide, X_MgO_oxide],
        dtype=float,
    )
    metal_values = np.asarray(
        [X_Fe_met, X_Si_met, X_Ni_met, X_O_met, X_C_met],
        dtype=float,
    )
    if np.any(~np.isfinite(oxide_values)) or np.any(
        (oxide_values < 0.0) | (oxide_values > 1.0)
    ):
        raise ValueError("Oxide mole fractions must be finite values in [0, 1].")
    if float(oxide_values.sum()) > 1.0 + 1e-10:
        raise ValueError("Tracked oxide mole fractions cannot sum to more than one.")
    if np.any(~np.isfinite(metal_values)) or np.any(
        (metal_values < 0.0) | (metal_values > 1.0)
    ):
        raise ValueError("Metal mole fractions must be finite values in [0, 1].")
    if float(metal_values.sum()) > 1.0 + 1e-10:
        raise ValueError("Tracked metal mole fractions cannot sum to more than one.")
    if (
        not math.isfinite(X_FeO_sil_phase)
        or X_FeO_sil_phase <= 0.0
        or X_FeO_sil_phase > 1.0
    ):
        raise ValueError(
            "Whole-silicate FeO mole fraction must be finite and in (0, 1]."
        )

    x_feo_oxide = float(X_FeO_oxide)
    x_cao_oxide = float(X_CaO_oxide)
    x_mgo_oxide = float(X_MgO_oxide)

    w_feo = (
        X_FeO_sil_phase * _COMPONENT_MOLAR_MASS["FeO"] / MW_silicate
    )
    w_fe = X_Fe_met * Elements["Fe"] / MW_metal
    w_si = X_Si_met * Elements["Si"] / MW_metal
    w_ni = X_Ni_met * Elements["Ni"] / MW_metal
    w_o = X_O_met * Elements["O"] / MW_metal
    w_c = X_C_met * Elements["C"] / MW_metal
    mass_fractions = {
        "FeO_silicate": w_feo,
        "Fe_metal": w_fe,
        "Si_metal": w_si,
        "Ni_metal": w_ni,
        "O_metal": w_o,
        "C_metal": w_c,
    }
    for name, value in mass_fractions.items():
        if not math.isfinite(value) or value < 0.0 or value >= 1.0:
            raise ValueError(f"{name} mass fraction must be finite and in [0, 1).")
    if w_feo <= 0.0:
        raise ValueError("FeO mass fraction must be positive in Boujibar Eq. 11.")
    if w_fe + w_si + w_ni + w_o + w_c > 1.0 + 1e-10:
        raise ValueError("Tracked metal mass fractions cannot sum to more than one.")

    # Boujibar et al. (2014), Eq. 6.
    # TODO(science): add +0.77 X_TiO2 +0.75 (X_Na2O + X_K2O) when tracked.
    log_cs = -5.704 + 3.15 * x_feo_oxide + 2.65 * x_cao_oxide + 0.12 * x_mgo_oxide

    # Boujibar et al. (2014), Eq. 11.  All composition terms here are mass
    # fractions, including FeO in silicate.
    log_one_minus_si = math.log10(1.0 - w_si)
    log_d_weight = (
        math.log10(w_feo)
        - log_cs
        + 405.0 / T
        + 136.0 * P / T
        + 32.0 * log_one_minus_si
        + 181.0 * log_one_minus_si**2
        + 305.0 * log_one_minus_si**3
        + 30.2 * math.log10(1.0 - w_c)
        + 1.13 * math.log10(1.0 - w_fe)
        + 10.7 * math.log10(1.0 - w_ni)
        + 31.4 * math.log10(1.0 - w_o)
        - 3.72
    )
    log_d_mole = (
        math.log(10.0) * log_d_weight
        + math.log(MW_metal)
        - math.log(MW_silicate)
    )
    if not math.isfinite(log_d_mole):
        raise ValueError("Boujibar log(D_S) is not finite.")
    return log_d_mole


def D_S_model_molar(
    P: float,
    T: float,
    *,
    X_FeO_oxide: float,
    X_CaO_oxide: float,
    X_MgO_oxide: float,
    X_FeO_sil_phase: float,
    X_Fe_met: float,
    X_Si_met: float,
    X_Ni_met: float,
    X_O_met: float,
    MW_metal: float,
    MW_silicate: float,
    X_C_met: float = 0.0,
) -> float:
    """Return Boujibar ``D_S`` through the existing molar interface.

    The empirical calculation is performed in log space and converted to a
    molar coefficient only at this public API boundary.
    """

    log_d_mole = _log_D_S_model_boujibar_molar(
        P,
        T,
        X_FeO_oxide=X_FeO_oxide,
        X_CaO_oxide=X_CaO_oxide,
        X_MgO_oxide=X_MgO_oxide,
        X_FeO_sil_phase=X_FeO_sil_phase,
        X_Fe_met=X_Fe_met,
        X_Si_met=X_Si_met,
        X_Ni_met=X_Ni_met,
        X_O_met=X_O_met,
        MW_metal=MW_metal,
        MW_silicate=MW_silicate,
        X_C_met=X_C_met,
    )
    try:
        d_mole = math.exp(log_d_mole)
    except OverflowError as exc:
        raise ValueError(
            "Boujibar D_S exceeds the representable floating-point range."
        ) from exc
    if not math.isfinite(d_mole) or d_mole <= 0.0:
        raise ValueError(
            "Boujibar D_S lies outside the positive floating-point range."
        )
    return d_mole


class ForwardKDOSolverSulfur:
    """Strict nested Fe-Ni-S metal-silicate partition solver.

    Conserved amounts are never clipped or transformed.  For fixed FeO and
    metal sulfur, NiO is obtained from a bounded root of the Ni exchange
    equation.  Every Ni trial solves the Si exchange quadratic and rejects
    roots that violate Si or O conservation.  Sulfur and oxygen exchange are
    then solved on their own physical intervals.

    ``nonneg``, ``enforce_z_box`` and ``enforce_d_nonneg`` remain accepted
    only so existing callers do not break.  They no longer alter conserved
    variables in the strict implementation.
    """

    def __init__(
        self,
        params: SulfurKDParams,
        *,
        P: float,
        T: float,
        sulfur_model: SulfurModel,
        sulfur_kd_calculator: Optional[SulfurKDCalculator] = None,
        sulfur_tol: float = 1e-12,
        sulfur_max_iter: int = 200,
        sulfur_grid_N: int = 20,
        coupled_max_iter: int = 100,
        coupled_tol: float = 1e-12,
        nonneg: Literal["clip", "softplus", "none"] = "clip",
        enforce_z_box: bool = False,
        enforce_d_nonneg: bool = False,
        physical_mode: Literal["strict", "diagnostic"] = "strict",
        allow_nearest: bool = False,
        ni_grid_N: int = 24,
        physical_tol: float = 1e-12,
        root_search: Literal["fast", "thorough"] = "fast",
    ) -> None:
        if not math.isfinite(P) or P < 0.0:
            raise ValueError("Pressure must be finite and non-negative.")
        if not math.isfinite(T) or T <= 0.0:
            raise ValueError("Temperature must be finite and positive.")
        if sulfur_model not in {"rose_weston", "boujibar"}:
            raise ValueError("sulfur_model must be 'rose_weston' or 'boujibar'.")
        if nonneg not in {"clip", "softplus", "none"}:
            raise ValueError("nonneg must be 'clip', 'softplus', or 'none'.")
        if physical_mode not in {"strict", "diagnostic"}:
            raise ValueError("physical_mode must be 'strict' or 'diagnostic'.")
        if root_search not in {"fast", "thorough"}:
            raise ValueError("root_search must be 'fast' or 'thorough'.")
        if sulfur_tol <= 0.0 or coupled_tol <= 0.0 or physical_tol <= 0.0:
            raise ValueError("Convergence and physical tolerances must be positive.")
        if (
            sulfur_max_iter < 1
            or coupled_max_iter < 1
            or sulfur_grid_N < 2
            or ni_grid_N < 2
        ):
            raise ValueError("Iteration limits and grid sizes must be positive.")
        if params.Fe_t <= 0.0:
            raise ValueError("The Fe-buffered differentiation model requires Fe_t > 0.")

        self.p = params
        self.P = float(P)
        self.T = float(T)
        self.sulfur_model = sulfur_model
        self.sulfur_kd_calculator = sulfur_kd_calculator or SulfurKDCalculator()
        self.sulfur_tol = float(sulfur_tol)
        self.sulfur_max_iter = int(sulfur_max_iter)
        self.sulfur_grid_N = int(sulfur_grid_N)
        self.coupled_max_iter = int(coupled_max_iter)
        self.coupled_tol = float(coupled_tol)
        self.ni_grid_N = int(ni_grid_N)
        self.physical_tol = float(physical_tol)
        self.physical_mode = physical_mode
        self.allow_nearest = bool(allow_nearest or physical_mode == "diagnostic")
        self.root_search = root_search
        self._inventory_scale = max(
            params.Fe_t
            + params.Ni_t
            + params.Si_t
            + params.u
            + params.m
            + params.n
            + params.S_t,
            np.finfo(float).tiny,
        )
        self._amount_tol = self.physical_tol * self._inventory_scale
        # Candidate filtering and final acceptance must use exactly the same
        # dimensionless Si-equilibrium tolerance.  A separate polynomial
        # residual previously made physical acceptance depend on warm-start
        # evaluation order for ill-conditioned, very small KD_Si states.
        self._major_equilibrium_tol = max(self.coupled_tol, 1e-10)
        self._last_y_hint: Optional[float] = None
        self._last_s_log_ratio_hint: Optional[float] = None
        self._si_root_cache: Dict[
            Tuple[float, ...],
            Tuple[float, float],
        ] = {}

        # Compatibility-only attributes.  No projection is performed.
        self.nonneg = nonneg
        self.enforce_z_box = enforce_z_box
        self.enforce_d_nonneg = enforce_d_nonneg

    @property
    def sulfur_params(self) -> SulfurKDParams:
        return self.p

    def _effective_search_subdivisions(self, requested: int) -> int:
        """Return the actual connected-search subdivision budget."""

        requested = max(int(requested), 2)
        return min(requested, 64) if self.root_search == "fast" else requested

    def _total_O(self) -> float:
        oxygen_total = (
            self.p.Fe_t + self.p.Ni_t + 2.0 * self.p.Si_t - self.p.O_L
        )
        if oxygen_total < -self._amount_tol:
            raise PhysicalStateError("Distributable oxygen inventory is negative.")
        return 0.0 if oxygen_total < 0.0 else oxygen_total

    @staticmethod
    def _stable_quadratic_roots(
        A: float,
        B: float,
        C: float,
        *,
        relative_tol: float = 1e-14,
        absolute_tol: float = 0.0,
    ) -> Tuple[float, ...]:
        """Return real roots using a cancellation-resistant formulation."""

        if not all(math.isfinite(value) for value in (A, B, C)):
            raise PhysicalStateError("Si quadratic coefficients are not finite.")
        # A, B, and C carry different powers of the composition scale, so
        # comparing their raw magnitudes would break invariance when all
        # inventories are rescaled. Exact zero is the only dimensionally
        # meaningful degeneracy test at this level.
        if A == 0.0:
            if B == 0.0:
                if C == 0.0:
                    raise AmbiguousPhysicalRootError(
                        "Si exchange equation is algebraically indeterminate."
                    )
                return ()
            return (-C / B,)

        discriminant = B * B - 4.0 * A * C
        discriminant_scale = max(
            B * B,
            abs(4.0 * A * C),
            np.finfo(float).tiny,
        )
        if discriminant < -relative_tol * discriminant_scale:
            return ()
        if discriminant < 0.0:
            discriminant = 0.0
        sqrt_discriminant = math.sqrt(discriminant)
        q = -0.5 * (B + math.copysign(sqrt_discriminant, B))
        if q == 0.0:
            return (-B / (2.0 * A),)
        roots = (q / A, C / q)
        if math.isclose(
            roots[0],
            roots[1],
            rel_tol=relative_tol,
            abs_tol=absolute_tol,
        ):
            return (roots[0],)
        return roots

    def _snap_numerical_boundary(
        self,
        value: float,
        lower: float,
        upper: float,
        *,
        name: str,
    ) -> float:
        scale = max(abs(lower), abs(upper), self._inventory_scale)
        tolerance = self.physical_tol * scale
        if value < lower - tolerance or value > upper + tolerance:
            raise PhysicalStateError(
                f"{name}={value!r} lies outside [{lower!r}, {upper!r}]."
            )
        if value < lower:
            return lower
        if value > upper:
            return upper
        return value

    def _solve_z(
        self,
        *,
        x_: float,
        y_: float,
        a_: float,
        b_: float,
        s_met: float,
        s_sil: float,
        z_hint: Optional[float] = None,
    ) -> Tuple[float, float]:
        """Solve and physically filter ``(SiO2_silicate, Si_metal)``.

        The quadratic is written in metal Si, ``c = Si_t - z``, rather than
        directly in silicate SiO2.  This retains full precision when the
        physical solution has ``z`` extremely close to ``Si_t``.
        """

        p = self.p
        cache_key = (
            x_,
            y_,
            a_,
            b_,
            s_met,
            s_sil,
            math.nan if z_hint is None else z_hint,
        )
        # NaN cannot be looked up reliably, so use a stable sentinel for the
        # default branch hint.
        if z_hint is None:
            cache_key = (*cache_key[:-1], -math.inf)
        cached = self._si_root_cache.get(cache_key)
        if cached is not None:
            return cached

        def remember(root: Tuple[float, float]) -> Tuple[float, float]:
            if len(self._si_root_cache) >= 8192:
                self._si_root_cache.clear()
            self._si_root_cache[cache_key] = root
            return root

        alpha = p.Si_t
        if alpha == 0.0:
            d_ = self._total_O() - x_ - y_
            if d_ < -self._amount_tol:
                raise PhysicalStateError("No oxygen remains after FeO and NiO.")
            return remember((0.0, 0.0))

        x2 = x_ * x_
        a2 = a_ * a_
        # With z = alpha-c:
        # d = O_total-x-y-2*alpha+2c,
        # N_metal = metal0+3c, and N_silicate = sil0-c.
        d0 = self._total_O() - x_ - y_ - 2.0 * alpha
        metal0 = a_ + b_ + s_met + d0
        sil0 = x_ + y_ + alpha + p.u + p.m + p.n + s_sil
        A = 3.0 * x2 - p.KD_Si * a2
        B = x2 * metal0 + p.KD_Si * a2 * (alpha + sil0)
        C = -p.KD_Si * a2 * alpha * sil0
        candidates = self._stable_quadratic_roots(
            A,
            B,
            C,
            absolute_tol=self._amount_tol,
        )
        physical: list[Tuple[float, float]] = []
        for candidate in candidates:
            if not math.isfinite(candidate):
                continue
            # Polish the analytic quadratic root before applying the physical
            # audit.  This is especially important when KD_Si is very small and
            # the desired root lies extremely close to the Si_t boundary.
            for _ in range(2):
                polynomial = A * candidate * candidate + B * candidate + C
                derivative = 2.0 * A * candidate + B
                if derivative == 0.0 or not math.isfinite(derivative):
                    break
                corrected = candidate - polynomial / derivative
                if not math.isfinite(corrected):
                    break
                candidate = corrected
            try:
                c_ = self._snap_numerical_boundary(
                    candidate,
                    0.0,
                    alpha,
                    name="Si_metal",
                )
            except PhysicalStateError:
                continue
            z_ = alpha - c_
            d_ = self._total_O() - x_ - y_ - 2.0 * z_
            if d_ < -self.physical_tol * max(
                self._total_O(),
                self._inventory_scale,
            ):
                continue
            if d_ < 0.0:
                d_ = 0.0
            metal_sum = a_ + b_ + c_ + d_ + s_met
            sil_sum = z_ + x_ + y_ + p.u + p.m + p.n + s_sil
            if metal_sum <= self._amount_tol or sil_sum <= self._amount_tol:
                continue
            si_equilibrium_residual = _normalized_difference(
                c_ * x2 * metal_sum,
                p.KD_Si * z_ * a2 * sil_sum,
            )
            if abs(si_equilibrium_residual) > self._major_equilibrium_tol:
                continue
            physical.append((z_, c_))

        if not physical:
            raise PhysicalStateError(
                "The Si exchange equation has no root satisfying Si and O conservation."
            )
        if len(physical) == 1:
            return remember(physical[0])

        hint = alpha if z_hint is None else z_hint
        distances = sorted(
            (abs(root[0] - hint), root)
            for root in physical
        )
        if len(distances) > 1 and math.isclose(
            distances[0][0],
            distances[1][0],
            rel_tol=0.0,
            abs_tol=self.physical_tol,
        ):
            raise AmbiguousPhysicalRootError(
                "Two Si roots are equally compatible with the branch hint."
            )
        return remember(distances[0][1])

    # Compatibility wrapper used by earlier diagnostic code.
    def _f_z_with_sulfur(
        self,
        x_: float,
        y_: float,
        a_: float,
        b_: float,
        s_met: float,
        s_sil: float,
    ) -> float:
        z_, _ = self._solve_z(
            x_=x_,
            y_=y_,
            a_=a_,
            b_=b_,
            s_met=s_met,
            s_sil=s_sil,
            z_hint=self.p.Si_t,
        )
        return z_

    def _state_for_y(
        self,
        *,
        x_: float,
        s_met: float,
        y_: float,
        s_sil: Optional[float] = None,
        z_hint: Optional[float] = None,
    ) -> Dict[str, float]:
        """Build one exact-conservation state for a NiO trial."""

        p = self.p
        x_ = self._snap_numerical_boundary(x_, 0.0, p.Fe_t, name="FeO_silicate")
        s_met = self._snap_numerical_boundary(
            s_met,
            0.0,
            p.S_t,
            name="S_metal",
        )
        y_ = self._snap_numerical_boundary(y_, 0.0, p.Ni_t, name="NiO_silicate")
        a_ = p.Fe_t - x_
        b_ = p.Ni_t - y_
        if s_sil is None:
            s_sil = p.S_t - s_met
        else:
            s_sil = self._snap_numerical_boundary(
                s_sil,
                0.0,
                p.S_t,
                name="S_silicate",
            )
            if abs((s_met + s_sil) - p.S_t) > self._amount_tol:
                raise PhysicalStateError("Explicit sulfur phase amounts do not conserve S.")
        z_, c_ = self._solve_z(
            x_=x_,
            y_=y_,
            a_=a_,
            b_=b_,
            s_met=s_met,
            s_sil=s_sil,
            z_hint=z_hint,
        )
        d_ = self._total_O() - x_ - y_ - 2.0 * z_
        d_scale = max(self._total_O(), self._inventory_scale)
        if d_ < -self.physical_tol * d_scale:
            raise PhysicalStateError("Metal oxygen is negative.")
        if d_ < 0.0:
            d_ = 0.0  # round-off only; larger violations were rejected above

        metal_sum = a_ + b_ + c_ + d_ + s_met
        sil_sum = x_ + y_ + z_ + p.u + p.m + p.n + s_sil
        if metal_sum <= self._amount_tol or sil_sum <= self._amount_tol:
            raise PhysicalStateError("Both phases must have positive molar totals.")
        ternary_sum = a_ + b_ + s_met
        if ternary_sum <= 0.0:
            raise PhysicalStateError("The Fe-Ni-S activity subsystem is empty.")
        x_fe_ternary = a_ / ternary_sum
        x_ni_ternary = b_ / ternary_sum
        x_s_ternary = s_met / ternary_sum
        if (
            x_fe_ternary <= 0.0
            or x_ni_ternary < 0.0
            or x_s_ternary < 0.0
        ):
            raise PhysicalStateError("Fe-Ni-S ternary composition is not admissible.")

        # Construct and evaluate the Ma model once per trial.  Si and O remain
        # part of the whole metal phase, but are intentionally excluded from
        # this Fe-Ni-S activity subsystem.
        activity_model = self.sulfur_kd_calculator.activity_model(
            X_Ni=x_ni_ternary,
            X_S=x_s_ternary,
        )
        ln_gamma_fe = activity_model.ln_gamma("Fe")
        ln_gamma_ni = activity_model.ln_gamma("Ni")
        ln_gamma_s = activity_model.ln_gamma("S")
        activities = {
            "Fe": math.exp(ln_gamma_fe),
            "Ni": math.exp(ln_gamma_ni),
            "S": math.exp(ln_gamma_s),
        }
        activity_factor = math.exp(ln_gamma_fe - ln_gamma_ni)
        kd_ni = p.KD_Ni * activity_factor
        ni_residual = _normalized_difference(
            b_ * x_,
            kd_ni * y_ * a_,
        )
        return {
            "x_": x_,
            "a_": a_,
            "y_": y_,
            "b_": b_,
            "z_": z_,
            "c_": c_,
            "d_": d_,
            "S_met_": s_met,
            "S_sil_": s_sil,
            "KD_Ni_effective": kd_ni,
            "KD_Ni_base": p.KD_Ni,
            "KD_Ni_activity_factor": activity_factor,
            "coupled_residual": ni_residual,
            "gamma_Fe_metal": activities["Fe"],
            "gamma_Ni_metal": activities["Ni"],
            "gamma_S_metal": activities["S"],
            "X_Fe_ternary": x_fe_ternary,
            "X_Ni_ternary": x_ni_ternary,
            "X_S_ternary": x_s_ternary,
        }

    @staticmethod
    def _find_brackets(
        samples: Sequence[Tuple[float, float, Dict[str, float]]],
        *,
        tolerance: float,
    ) -> Tuple[
        list[Tuple[float, float, Dict[str, float]]],
        list[Tuple[Tuple[float, float, Dict[str, float]], Tuple[float, float, Dict[str, float]]]],
    ]:
        exact: list[Tuple[float, float, Dict[str, float]]] = []
        brackets: list[
            Tuple[
                Tuple[float, float, Dict[str, float]],
                Tuple[float, float, Dict[str, float]],
            ]
        ] = []
        previous: Optional[Tuple[float, float, Dict[str, float]]] = None
        for sample in samples:
            _, residual, _ = sample
            if abs(residual) <= tolerance:
                exact.append(sample)
            if previous is not None and previous[1] * residual < 0.0:
                brackets.append((previous, sample))
            previous = sample
        return exact, brackets

    def _discover_connected_roots(
        self,
        sample: Callable[
            [float],
            Optional[Tuple[float, float, Dict[str, float]]],
        ],
        lower: float,
        upper: float,
        *,
        tolerance: float,
        subdivisions: int,
        anchors: Sequence[float] = (),
    ) -> Tuple[
        list[Tuple[float, float, Dict[str, float]]],
        list[
            Tuple[
                Tuple[float, float, Dict[str, float]],
                Tuple[float, float, Dict[str, float]],
            ]
        ],
        list[Tuple[float, float, Dict[str, float]]],
    ]:
        """Sample a bounded interval without joining across invalid states.

        Fast mode evaluates both ends, the midpoint, and supplied warm-start
        anchors. A bounded fallback scan is used only when those trials do not
        locate a root. Thorough mode performs that scan immediately. Brackets
        are formed only between adjacent attempted coordinates that are both
        physical, so an observed invalid gap can never be crossed.
        """

        attempted: Dict[
            float,
            Optional[Tuple[float, float, Dict[str, float]]],
        ] = {}

        def attempt(coordinate: float) -> None:
            coordinate = min(max(float(coordinate), lower), upper)
            if coordinate not in attempted:
                attempted[coordinate] = sample(coordinate)

        def classify() -> Tuple[
            list[Tuple[float, float, Dict[str, float]]],
            list[
                Tuple[
                    Tuple[float, float, Dict[str, float]],
                    Tuple[float, float, Dict[str, float]],
                ]
            ],
            list[Tuple[float, float, Dict[str, float]]],
        ]:
            exact: list[Tuple[float, float, Dict[str, float]]] = []
            brackets: list[
                Tuple[
                    Tuple[float, float, Dict[str, float]],
                    Tuple[float, float, Dict[str, float]],
                ]
            ] = []
            valid = [
                value
                for _, value in sorted(attempted.items())
                if value is not None
            ]
            ordered = sorted(attempted)
            for coordinate in ordered:
                current = attempted[coordinate]
                if current is not None and abs(current[1]) <= tolerance:
                    exact.append(current)
            for left_x, right_x in zip(ordered, ordered[1:]):
                left = attempted[left_x]
                right = attempted[right_x]
                if left is None or right is None:
                    continue
                if left[1] * right[1] < 0.0:
                    brackets.append((left, right))
            return exact, brackets, valid

        attempt(lower)
        attempt(upper)
        attempt(0.5 * (lower + upper))
        for anchor in anchors:
            attempt(anchor)
            if self.root_search == "fast":
                exact, brackets, valid = classify()
                if exact or brackets:
                    return exact, brackets, valid
        exact, brackets, valid = classify()
        if self.root_search == "thorough" or (not exact and not brackets):
            bounded_n = self._effective_search_subdivisions(subdivisions)
            for index in range(1, bounded_n):
                attempt(lower + (upper - lower) * index / bounded_n)
            exact, brackets, valid = classify()
        return exact, brackets, valid

    @staticmethod
    def _solve_bracketed_state_root(
        evaluate: Callable[[float], Tuple[float, Dict[str, float]]],
        left: Tuple[float, float, Dict[str, float]],
        right: Tuple[float, float, Dict[str, float]],
        *,
        residual_tol: float,
        coordinate_tol: float,
        max_iter: int,
    ) -> Tuple[float, float, Dict[str, float], int, bool]:
        lo, flo, slo = left
        hi, fhi, shi = right
        if flo * fhi > 0.0:
            raise ValueError("Bracketed root solve requires a sign change.")

        cache: Dict[float, Tuple[float, Dict[str, float]]] = {
            lo: (flo, slo),
            hi: (fhi, shi),
        }

        def objective(coordinate: float) -> float:
            if coordinate not in cache:
                cache[coordinate] = evaluate(coordinate)
            return cache[coordinate][0]

        if abs(flo) <= residual_tol:
            return lo, flo, slo, 0, True
        if abs(fhi) <= residual_tol:
            return hi, fhi, shi, 0, True

        iterations = 0
        try:
            root, info = toms748(
                objective,
                lo,
                hi,
                xtol=max(coordinate_tol, np.finfo(float).tiny),
                rtol=4.0 * np.finfo(float).eps,
                maxiter=max_iter,
                full_output=True,
                disp=False,
            )
            iterations = int(info.iterations)
            objective(root)
        except (ArithmeticError, RuntimeError, ValueError):
            # A supposedly connected physical bracket can still contain an
            # undiscovered invalid point. Fall back to safeguarded bisection
            # so the caller receives the same strict failure semantics.
            for iterations in range(1, max_iter + 1):
                mid = 0.5 * (lo + hi)
                if mid == lo or mid == hi:
                    break
                fm, sm = evaluate(mid)
                cache[mid] = (fm, sm)
                if abs(fm) <= residual_tol:
                    return mid, fm, sm, iterations, True
                if flo * fm <= 0.0:
                    hi, fhi, shi = mid, fm, sm
                else:
                    lo, flo, slo = mid, fm, sm

        best_x, (best_f, best_state) = min(
            cache.items(),
            key=lambda item: abs(item[1][0]),
        )
        return (
            best_x,
            best_f,
            best_state,
            iterations,
            abs(best_f) <= residual_tol,
        )

    def _solve_ni_state(
        self,
        x_: float,
        s_met: float,
        s_sil: Optional[float] = None,
    ) -> Dict[str, float]:
        """Solve the composition-dependent Ni exchange equation on its box."""

        p = self.p
        if p.Ni_t == 0.0:
            state = self._state_for_y(
                x_=x_,
                s_met=s_met,
                s_sil=s_sil,
                y_=0.0,
            )
            state.update(
                coupled_iterations=0,
                coupled_converged=abs(state["coupled_residual"]) <= self.coupled_tol,
            )
            return state

        def sample(y_value: float) -> Optional[Tuple[float, float, Dict[str, float]]]:
            try:
                state = self._state_for_y(
                    x_=x_,
                    s_met=s_met,
                    s_sil=s_sil,
                    y_=y_value,
                    z_hint=p.Si_t,
                )
            except SulfurSolverError:
                return None
            return y_value, state["coupled_residual"], state

        denominator = (p.Fe_t - x_) * p.KD_Ni + x_
        y_hint = (
            x_ * p.Ni_t / denominator
            if denominator > 0.0
            else 0.5 * p.Ni_t
        )
        anchors = [y_hint]
        if self._last_y_hint is not None:
            anchors.insert(0, self._last_y_hint)
        search_subdivisions = self._effective_search_subdivisions(self.ni_grid_N)
        exact, brackets, valid_samples = self._discover_connected_roots(
            sample,
            0.0,
            p.Ni_t,
            tolerance=self.coupled_tol,
            subdivisions=search_subdivisions,
            anchors=anchors,
        )

        if not valid_samples:
            raise NiPartitionNotReachableError(
                f"No physical NiO trial state was found within the "
                f"{self.root_search} search budget "
                f"({search_subdivisions} subdivisions)."
            )
        if exact:
            y_value, _, state = min(exact, key=lambda item: item[0])
            state.update(coupled_iterations=0, coupled_converged=True)
            self._last_y_hint = y_value
            return state

        if brackets:
            left, right = min(
                brackets,
                key=lambda pair: abs(0.5 * (pair[0][0] + pair[1][0]) - y_hint),
            )

            def evaluate(y_value: float) -> Tuple[float, Dict[str, float]]:
                state = self._state_for_y(
                    x_=x_,
                    s_met=s_met,
                    s_sil=s_sil,
                    y_=y_value,
                    z_hint=p.Si_t,
                )
                return state["coupled_residual"], state

            _, residual, state, iterations, converged = self._solve_bracketed_state_root(
                evaluate,
                left,
                right,
                residual_tol=self.coupled_tol,
                coordinate_tol=self.physical_tol * max(
                    p.Ni_t,
                    self._inventory_scale,
                ),
                max_iter=self.coupled_max_iter,
            )
            state.update(
                coupled_residual=residual,
                coupled_iterations=iterations,
                coupled_converged=converged,
            )
            if converged or self.allow_nearest:
                self._last_y_hint = state["y_"]
                return state
            raise NiPartitionNotReachableError(
                f"Ni root did not converge; best residual={residual:.3e}."
            )

        best = min(valid_samples, key=lambda item: abs(item[1]))
        state = best[2]
        state.update(coupled_iterations=0, coupled_converged=False)
        if self.allow_nearest:
            self._last_y_hint = state["y_"]
            return state
        raise NiPartitionNotReachableError(
            f"No connected Ni root was found within the {self.root_search} "
            f"search budget ({search_subdivisions} subdivisions); "
            f"best sampled residual={best[1]:.3e}."
        )

    def _major_state_for_sulfur(self, x_init: float, s_met: float) -> Dict[str, float]:
        """Compatibility name for the strict fixed-S major-element solve."""

        return self._solve_ni_state(x_init, s_met)

    def _phase_properties(
        self,
        sol: Dict[str, float],
        p: SulfurKDParams,
    ) -> Dict[str, float]:
        amount_names = (
            "x_",
            "a_",
            "y_",
            "b_",
            "z_",
            "c_",
            "d_",
            "S_met_",
            "S_sil_",
        )
        if any(
            not math.isfinite(sol[name]) or sol[name] < -self._amount_tol
            for name in amount_names
        ):
            raise PhysicalStateError("Phase amounts must be finite and non-negative.")
        metal_sum = sol["a_"] + sol["b_"] + sol["c_"] + sol["d_"] + sol["S_met_"]
        sil_sum = (
            sol["x_"]
            + sol["y_"]
            + sol["z_"]
            + p.u
            + p.m
            + p.n
            + sol["S_sil_"]
        )
        oxide_sum = sol["x_"] + sol["y_"] + sol["z_"] + p.u + p.m + p.n
        if (
            not math.isfinite(metal_sum)
            or not math.isfinite(sil_sum)
            or metal_sum <= 0.0
            or sil_sum <= 0.0
        ):
            raise PhysicalStateError(
                "Metal and silicate phase totals must be finite and positive."
            )
        metal_mass = (
            sol["a_"] * Elements["Fe"]
            + sol["b_"] * Elements["Ni"]
            + sol["c_"] * Elements["Si"]
            + sol["d_"] * Elements["O"]
            + sol["S_met_"] * Elements["S"]
        )
        sil_mass = (
            sol["x_"] * _COMPONENT_MOLAR_MASS["FeO"]
            + sol["y_"] * _COMPONENT_MOLAR_MASS["NiO"]
            + sol["z_"] * _COMPONENT_MOLAR_MASS["SiO2"]
            + p.u * _COMPONENT_MOLAR_MASS["MgO"]
            + p.m * _COMPONENT_MOLAR_MASS["AlO1.5"]
            + p.n * _COMPONENT_MOLAR_MASS["CaO"]
            + sol["S_sil_"] * Elements["S"]
        )
        return {
            "metal_sum": metal_sum,
            "sil_sum": sil_sum,
            "oxide_sum": oxide_sum,
            "MW_metal": metal_mass / metal_sum,
            "MW_silicate": sil_mass / sil_sum,
        }

    def _sulfur_coefficients(
        self,
        sol: Dict[str, float],
    ) -> Tuple[float, float, float]:
        """Return ``(D_weight, D_mole, ln(D_mole))`` for one state."""

        p = self.p
        phase = self._phase_properties(sol, p)
        metal_sum = phase["metal_sum"]
        sil_sum = phase["sil_sum"]
        oxide_sum = phase["oxide_sum"]
        if self.sulfur_model == "rose_weston":
            log_d_mole = _log_D_S_model_rose_weston_molar(
                self.P,
                self.T,
                X_Fe_met=sol["a_"] / metal_sum,
                X_Ni_met=sol["b_"] / metal_sum,
                X_O_met=sol["d_"] / metal_sum,
                X_Si_met=sol["c_"] / metal_sum,
                X_S_met=sol["S_met_"] / metal_sum,
                X_FeO_sil=sol["x_"] / sil_sum,
                X_NiO_sil=sol["y_"] / sil_sum,
                X_SiO2_sil=sol["z_"] / sil_sum,
                X_MgO_sil=p.u / sil_sum,
                X_AlO1_5_sil=p.m / sil_sum,
                X_CaO_sil=p.n / sil_sum,
                X_S_sil=sol["S_sil_"] / sil_sum,
            )
        else:
            if not math.isfinite(oxide_sum) or oxide_sum <= 0.0:
                raise PhysicalStateError(
                    "Boujibar sulfide capacity requires a positive oxide total."
                )
            log_d_mole = _log_D_S_model_boujibar_molar(
                self.P,
                self.T,
                X_FeO_oxide=sol["x_"] / oxide_sum,
                X_CaO_oxide=p.n / oxide_sum,
                X_MgO_oxide=p.u / oxide_sum,
                X_FeO_sil_phase=sol["x_"] / sil_sum,
                X_Fe_met=sol["a_"] / metal_sum,
                X_Si_met=sol["c_"] / metal_sum,
                X_Ni_met=sol["b_"] / metal_sum,
                X_O_met=sol["d_"] / metal_sum,
                MW_metal=phase["MW_metal"],
                MW_silicate=phase["MW_silicate"],
            )
        try:
            d_mole = math.exp(log_d_mole)
            d_weight = math.exp(
                log_d_mole
                + math.log(phase["MW_silicate"])
                - math.log(phase["MW_metal"])
            )
        except OverflowError as exc:
            raise PhysicalStateError(
                "D_S exceeds the representable output range at this state."
            ) from exc
        if (
            not math.isfinite(d_mole)
            or not math.isfinite(d_weight)
            or d_mole <= 0.0
            or d_weight <= 0.0
        ):
            raise PhysicalStateError(
                "D_S lies outside the positive representable output range."
            )
        return d_weight, d_mole, log_d_mole

    @staticmethod
    def _sulfur_amounts_from_log_ratio(
        sulfur_total: float,
        log_ratio: float,
    ) -> Tuple[float, float]:
        """Return conserved S amounts without subtractive cancellation.

        ``log_ratio = ln(S_metal / S_silicate)``.  The branch-specific
        formulas retain a representable trace phase even when the ratio is
        hundreds of natural-log units from unity.
        """

        if sulfur_total <= 0.0 or not math.isfinite(sulfur_total):
            raise PhysicalStateError(
                "Positive finite total sulfur is required in log-ratio space."
            )
        if not math.isfinite(log_ratio):
            raise PhysicalStateError("Sulfur log ratio must be finite.")

        absolute_ratio = abs(log_ratio)
        inverse_ratio = math.exp(-absolute_ratio)
        log_denominator = math.log1p(inverse_ratio)
        major = sulfur_total / (1.0 + inverse_ratio)
        # Form the trace amount from logarithms. This avoids prematurely
        # underflowing exp(-q) before multiplication by a large inventory.
        trace = math.exp(
            math.log(sulfur_total) - absolute_ratio - log_denominator
        )
        return (major, trace) if log_ratio >= 0.0 else (trace, major)

    def _evaluate_sulfur(
        self,
        x_: float,
        s_met: float,
        s_sil: Optional[float] = None,
    ) -> Dict[str, float]:
        state = self._solve_ni_state(x_, s_met, s_sil=s_sil)
        if not state["coupled_converged"] and not self.allow_nearest:
            raise NiPartitionNotReachableError("Ni did not converge at sulfur trial.")
        d_weight, d_mole, log_d_mole = self._sulfur_coefficients(state)
        phase = self._phase_properties(state, self.p)
        if state["S_met_"] <= 0.0:
            residual = -math.inf
        elif state["S_sil_"] <= 0.0:
            residual = math.inf
        else:
            residual = (
                math.log(state["S_met_"] / state["S_sil_"])
                + math.log(phase["sil_sum"] / phase["metal_sum"])
                - log_d_mole
            )
        state.update(
            D_S_weight=d_weight,
            D_S_mole=d_mole,
            _log_D_S_mole=log_d_mole,
            sulfur_residual=residual,
        )
        return state

    def _evaluate_sulfur_log_ratio(
        self,
        x_: float,
        log_ratio: float,
    ) -> Dict[str, float]:
        s_met, s_sil = self._sulfur_amounts_from_log_ratio(
            self.p.S_t,
            log_ratio,
        )
        state = self._evaluate_sulfur(
            x_,
            s_met,
            s_sil=s_sil,
        )
        # Use the coordinate itself rather than recomputing S_met/S_sil;
        # this remains accurate at ratios too extreme for direct division.
        phase = self._phase_properties(state, self.p)
        state["sulfur_residual"] = (
            log_ratio
            + math.log(phase["sil_sum"] / phase["metal_sum"])
            - state["_log_D_S_mole"]
        )
        return state

    def _finalize_forward_state(
        self,
        state: Dict[str, float],
        *,
        sulfur_iterations: int,
        sulfur_converged: bool,
        failure_reason: Optional[str] = None,
    ) -> Dict[str, float]:
        state.update(
            sulfur_iterations=sulfur_iterations,
            sulfur_converged=sulfur_converged,
            sulfur_model=self.sulfur_model,
            failure_reason=failure_reason,
        )
        report = audit_physical_state(
            state,
            self.p,
            physical_tol=self.physical_tol,
        )
        major_names = (
            "Fe_balance",
            "Ni_balance",
            "Si_balance",
            "O_balance",
            "Ni_equilibrium",
            "Si_equilibrium",
        )
        major_converged = report.physical and all(
            abs(report.residuals[name]) <= self._major_equilibrium_tol
            for name in major_names
        )
        physical_converged = report.physical
        overall = major_converged and sulfur_converged and physical_converged
        if failure_reason is None and not overall:
            if not physical_converged:
                failure_reason = "Phase amounts failed physical admissibility."
            elif not major_converged:
                failure_reason = "Major-element partition residual did not converge."
            elif not sulfur_converged:
                failure_reason = "Sulfur partition residual did not converge."
        state.update(
            residuals=dict(report.residuals),
            major_converged=major_converged,
            oxygen_converged=False,
            iw_converged=False,
            physical_converged=physical_converged,
            overall_converged=overall,
            failure_reason=failure_reason,
            model_assumptions=MODEL_ASSUMPTIONS,
        )
        if not overall and not self.allow_nearest:
            raise PhysicalStateError(
                failure_reason
                or f"Forward state failed strict residual audit: {report.residuals!r}"
            )
        return state

    def forward_solve(self, x_init: float) -> Dict[str, float]:
        """Solve Ni, Si, O, and sulfur at one physically admissible FeO amount."""

        p = self.p
        if not math.isfinite(x_init):
            raise ValueError("FeO trial must be finite.")
        x_upper = min(p.Fe_t, self._total_O())
        x_ = self._snap_numerical_boundary(
            float(x_init),
            0.0,
            x_upper,
            name="FeO_silicate",
        )

        if p.S_t == 0.0:
            state = self._solve_ni_state(x_, 0.0)
            state.update(
                D_S_weight=0.0,
                D_S_mole=0.0,
                sulfur_residual=0.0,
            )
            return self._finalize_forward_state(
                state,
                sulfur_iterations=0,
                sulfur_converged=bool(state["coupled_converged"]),
            )

        # q = ln(S_metal/S_silicate) avoids deleting legitimate solutions
        # arbitrarily close to either sulfur boundary. Its admissible interval
        # is determined by the smallest positive amount that can be represented
        # for this sulfur inventory, rather than by an arbitrary fixed number.
        minimum_positive = float(np.nextafter(0.0, 1.0))
        q_capacity = math.log(p.S_t) - math.log(minimum_positive)
        if not math.isfinite(q_capacity) or q_capacity <= 0.0:
            raise SulfurPartitionNotReachableError(
                "Total sulfur is too small to represent two positive phase amounts."
            )
        q_limit = math.nextafter(q_capacity, 0.0)

        def sample(q_value: float) -> Optional[Tuple[float, float, Dict[str, float]]]:
            try:
                state = self._evaluate_sulfur_log_ratio(x_, q_value)
            except (SulfurSolverError, ValueError, ArithmeticError):
                return None
            return q_value, state["sulfur_residual"], state

        anchors: list[float] = []
        if self._last_s_log_ratio_hint is not None:
            anchors.append(self._last_s_log_ratio_hint)
        # Most experimental D_S values correspond to |q| of order 1--20.
        # Geometric anchors expand to the machine-derived boundary, providing
        # bounded extreme-ratio coverage without a dense global grid.
        magnitude = 1.0
        while magnitude < q_limit:
            anchors.extend((-magnitude, magnitude))
            magnitude *= 2.0
        search_subdivisions = self._effective_search_subdivisions(
            self.sulfur_grid_N
        )
        exact, brackets, valid_samples = self._discover_connected_roots(
            sample,
            -q_limit,
            q_limit,
            tolerance=self.sulfur_tol,
            subdivisions=search_subdivisions,
            anchors=anchors,
        )

        if not valid_samples:
            raise SulfurPartitionNotReachableError(
                f"No physical sulfur trial state was found within the "
                f"{self.root_search} search budget "
                f"({search_subdivisions} subdivisions) and representable "
                f"q range [{-q_limit:.6g}, {q_limit:.6g}]."
            )
        if exact:
            q_value, _, state = min(
                exact,
                key=lambda item: abs(
                    item[0]
                    - (
                        self._last_s_log_ratio_hint
                        if self._last_s_log_ratio_hint is not None
                        else 0.0
                    )
                ),
            )
            self._last_s_log_ratio_hint = q_value
            return self._finalize_forward_state(
                state,
                sulfur_iterations=0,
                sulfur_converged=True,
            )

        if brackets:
            q_hint = (
                self._last_s_log_ratio_hint
                if self._last_s_log_ratio_hint is not None
                else 0.0
            )
            left, right = min(
                brackets,
                key=lambda pair: abs(
                    0.5 * (pair[0][0] + pair[1][0]) - q_hint
                ),
            )

            def evaluate(q_value: float) -> Tuple[float, Dict[str, float]]:
                state = self._evaluate_sulfur_log_ratio(x_, q_value)
                return state["sulfur_residual"], state

            q_value, residual, state, iterations, converged = (
                self._solve_bracketed_state_root(
                    evaluate,
                    left,
                    right,
                    residual_tol=self.sulfur_tol,
                    coordinate_tol=max(self.physical_tol, 1e-14),
                    max_iter=self.sulfur_max_iter,
                )
            )
            state["sulfur_residual"] = residual
            if converged or self.allow_nearest:
                self._last_s_log_ratio_hint = q_value
                return self._finalize_forward_state(
                    state,
                    sulfur_iterations=iterations,
                    sulfur_converged=converged and bool(state["coupled_converged"]),
                    failure_reason=None if converged else "Sulfur root did not converge.",
                )
            raise SulfurPartitionNotReachableError(
                f"Sulfur root did not converge; best residual={residual:.3e}."
            )

        best = min(valid_samples, key=lambda item: abs(item[1]))
        if self.allow_nearest:
            return self._finalize_forward_state(
                best[2],
                sulfur_iterations=0,
                sulfur_converged=False,
                failure_reason="No sulfur sign-changing bracket.",
            )
        raise SulfurPartitionNotReachableError(
            f"No connected sulfur root was found within the {self.root_search} "
            f"search budget ({search_subdivisions} subdivisions) and "
            f"representable q range [{-q_limit:.6g}, {q_limit:.6g}]; "
            f"best sampled residual={best[1]:.3e}."
        )

    def compute_KD_O_from_x(self, x_: float) -> Tuple[float, Dict[str, float]]:
        """Compute Brennan Eq. S4 only from a strictly physical nested state."""

        state = self.forward_solve(x_)
        if not state["overall_converged"]:
            raise PhysicalStateError("Nested state is not fully converged.")
        phase = self._phase_properties(state, self.p)
        x_feo = state["x_"] / phase["sil_sum"]
        if x_feo <= 0.0:
            raise PhysicalStateError("Oxygen KD is undefined at zero silicate FeO.")

        # TODO(science): Fischer's explicit sulfur correction to oxygen
        # partitioning remains intentionally deferred.
        kd_o = (
            state["a_"]
            * state["d_"]
            / (x_feo * phase["metal_sum"] ** 2)
        )
        if not math.isfinite(kd_o) or kd_o < 0.0:
            raise PhysicalStateError("Calculated oxygen KD is not finite/non-negative.")
        return kd_o, state

    def _result(
        self,
        kd_o: float,
        target: float,
        state: Dict[str, float],
        *,
        oxygen_tol: float,
        oxygen_converged: bool,
        failure_reason: Optional[str] = None,
    ) -> SulfurKDResult:
        audit_state = dict(state)
        audit_state["KD_O_"] = kd_o
        report = audit_physical_state(
            audit_state,
            self.p,
            KD_O_target=target,
            physical_tol=self.physical_tol,
        )
        log_residual = report.residuals.get("O_KD_log", math.inf)
        oxygen_converged = (
            oxygen_converged
            and math.isfinite(log_residual)
            and abs(log_residual) <= oxygen_tol
        )
        overall = (
            bool(state["overall_converged"])
            and oxygen_converged
            and report.physical
        )
        if not overall and not self.allow_nearest:
            raise OxygenKDNotReachableError(
                failure_reason
                or f"Oxygen KD residual failed strict audit: {log_residual:.3e}."
            )
        return SulfurKDResult(
            x_=state["x_"],
            a_=state["a_"],
            y_=state["y_"],
            b_=state["b_"],
            z_=state["z_"],
            c_=state["c_"],
            d_=state["d_"],
            KD_O_=kd_o,
            residual=kd_o - target,
            S_met_=state["S_met_"],
            S_sil_=state["S_sil_"],
            D_S_weight=state["D_S_weight"],
            D_S_mole=state["D_S_mole"],
            sulfur_residual=state["sulfur_residual"],
            sulfur_iterations=int(state["sulfur_iterations"]),
            sulfur_converged=bool(state["sulfur_converged"]),
            sulfur_model=str(state["sulfur_model"]),
            coupled_residual=state["coupled_residual"],
            coupled_iterations=int(state["coupled_iterations"]),
            coupled_converged=bool(state["coupled_converged"]),
            KD_Ni_effective=state["KD_Ni_effective"],
            KD_Ni_base=state["KD_Ni_base"],
            KD_Ni_activity_factor=state["KD_Ni_activity_factor"],
            gamma_Fe_metal=state["gamma_Fe_metal"],
            gamma_Ni_metal=state["gamma_Ni_metal"],
            gamma_S_metal=state["gamma_S_metal"],
            X_Fe_ternary=state["X_Fe_ternary"],
            X_Ni_ternary=state["X_Ni_ternary"],
            X_S_ternary=state["X_S_ternary"],
            major_converged=bool(state["major_converged"]),
            oxygen_converged=oxygen_converged,
            physical_converged=report.physical,
            overall_converged=overall,
            failure_reason=failure_reason,
            residuals=dict(report.residuals),
            model_assumptions=MODEL_ASSUMPTIONS,
        )

    def solve_x_for_KD_O(
        self,
        KD_O_target: float,
        *,
        tol: float = 1e-12,
        max_iter: int = 200,
        grid_N: int = 200,
        x_hint: Optional[float] = None,
    ) -> SulfurKDResult:
        """Strictly solve ``log10(KD_O / KD_O_target) = 0``."""

        if not math.isfinite(KD_O_target) or KD_O_target <= 0.0:
            raise ValueError("KD_O_target must be finite and positive.")
        if not math.isfinite(tol) or tol <= 0.0:
            raise ValueError("tol must be finite and positive.")
        if max_iter < 1 or grid_N < 2:
            raise ValueError("max_iter must be positive and grid_N >= 2.")

        x_upper = min(self.p.Fe_t, self._total_O())
        if x_upper <= 0.0:
            raise OxygenKDNotReachableError("No positive FeO interval exists.")
        x_edge = max(np.nextafter(0.0, 1.0), x_upper * 1e-14)
        valid_samples: list[Tuple[float, float, Dict[str, float]]] = []
        sample_by_x: Dict[float, Tuple[float, float, Dict[str, float]]] = {}
        kd_by_x: Dict[float, float] = {}
        attempted: Dict[float, bool] = {}

        def sample(x_value: float) -> None:
            if x_value in attempted:
                return
            attempted[x_value] = False
            try:
                kd_o, state = self.compute_KD_O_from_x(x_value)
            except (SulfurSolverError, ValueError, ArithmeticError, FloatingPointError):
                return
            if kd_o <= 0.0:
                return
            residual = math.log10(kd_o / KD_O_target)
            kd_by_x[x_value] = kd_o
            sampled = (x_value, residual, state)
            valid_samples.append(sampled)
            sample_by_x[x_value] = sampled
            attempted[x_value] = True

        hint = self.p.Fe_t * 0.5 if x_hint is None else float(x_hint)
        # Evaluate the caller's continuation state first. A set followed by
        # sorting previously destroyed this priority and allowed unrelated
        # boundary probes to change warm-start rounding before an exact hint.
        primary_x = list(
            dict.fromkeys(
                (
                    min(max(hint, x_edge), x_upper * (1.0 - 1e-10)),
                    x_edge,
                    x_upper * (1.0 - 1e-10),
                    x_upper * 0.25,
                    x_upper * 0.50,
                    x_upper * 0.75,
                )
            )
        )
        for x_value in primary_x:
            sample(x_value)

        # A physical interval commonly terminates where metal oxygen reaches
        # zero.  KD_O can cross its target in a very narrow region immediately
        # before that boundary, so refine valid/invalid transitions instead
        # of globally increasing the grid density.
        def refine_valid_boundaries() -> None:
            ordered_attempts = sorted(attempted)
            transitions = [
                (left, right)
                for left, right in zip(ordered_attempts, ordered_attempts[1:])
                if attempted[left] != attempted[right]
            ]
            transitions.sort(
                key=lambda pair: abs(
                    sample_by_x[
                        pair[0] if attempted[pair[0]] else pair[1]
                    ][1]
                )
            )
            for left_x, right_x in transitions:
                if attempted[left_x]:
                    valid = sample_by_x[left_x]
                    invalid_x = right_x
                else:
                    invalid_x = left_x
                    valid = sample_by_x[right_x]
                for _ in range(min(max_iter, 32)):
                    midpoint = 0.5 * (valid[0] + invalid_x)
                    if midpoint == valid[0] or midpoint == invalid_x:
                        break
                    sample(midpoint)
                    if attempted[midpoint]:
                        midpoint_sample = sample_by_x[midpoint]
                        # The purpose of this refinement is to locate an
                        # oxygen-KD bracket, not to resolve the physical
                        # boundary itself.  Stop as soon as the target is
                        # bracketed; continuing would repeat dozens of full
                        # nested sulfur/Ni solves without improving the root.
                        if (
                            abs(midpoint_sample[1]) <= tol
                            or valid[1] * midpoint_sample[1] < 0.0
                        ):
                            break
                        valid = midpoint_sample
                    else:
                        invalid_x = midpoint
                    if abs(invalid_x - valid[0]) <= (
                        self.physical_tol
                        * max(x_upper, self._inventory_scale)
                    ):
                        break

        def classify_connected() -> Tuple[
            list[Tuple[float, float, Dict[str, float]]],
            list[
                Tuple[
                    Tuple[float, float, Dict[str, float]],
                    Tuple[float, float, Dict[str, float]],
                ]
            ],
        ]:
            ordered_attempts = sorted(attempted)
            exact_samples = [
                sample_by_x[x_value]
                for x_value in ordered_attempts
                if attempted[x_value] and abs(sample_by_x[x_value][1]) <= tol
            ]
            connected_brackets = []
            for left_x, right_x in zip(
                ordered_attempts,
                ordered_attempts[1:],
            ):
                if not attempted[left_x] or not attempted[right_x]:
                    continue
                left = sample_by_x[left_x]
                right = sample_by_x[right_x]
                if left[1] * right[1] < 0.0:
                    connected_brackets.append((left, right))
            return exact_samples, connected_brackets

        def finish(
            sampled: Tuple[float, float, Dict[str, float]],
        ) -> SulfurKDResult:
            return self._result(
                kd_by_x[sampled[0]],
                KD_O_target,
                sampled[2],
                oxygen_tol=tol,
                oxygen_converged=abs(sampled[1]) <= tol,
                failure_reason=(
                    None
                    if abs(sampled[1]) <= tol
                    else "Oxygen KD root did not converge."
                ),
            )

        attempted_brackets: set[Tuple[float, float]] = set()

        def try_current_roots() -> Optional[SulfurKDResult]:
            exact_samples, connected_brackets = classify_connected()
            if exact_samples:
                chosen = min(
                    exact_samples,
                    key=(
                        (lambda item: abs(item[0] - x_hint))
                        if x_hint is not None
                        else (lambda item: item[0])
                    ),
                )
                return finish(chosen)

            for left, right in sorted(
                connected_brackets,
                key=lambda pair: abs(
                    0.5 * (pair[0][0] + pair[1][0]) - hint
                ),
            ):
                marker = (left[0], right[0])
                if marker in attempted_brackets:
                    continue
                attempted_brackets.add(marker)

                def evaluate_x(
                    x_value: float,
                ) -> Tuple[float, Dict[str, float]]:
                    sample(x_value)
                    if not attempted[x_value]:
                        raise PhysicalStateError(
                            "Oxygen bracket crossed an invalid physical state."
                        )
                    sampled = sample_by_x[x_value]
                    return sampled[1], sampled[2]

                try:
                    root_x, _, _, _, converged = self._solve_bracketed_state_root(
                        evaluate_x,
                        left,
                        right,
                        residual_tol=tol,
                        coordinate_tol=self.physical_tol
                        * max(x_upper, self._inventory_scale),
                        max_iter=max_iter,
                    )
                except (SulfurSolverError, ValueError, ArithmeticError):
                    continue
                if converged:
                    return finish(sample_by_x[root_x])
            return None

        refine_valid_boundaries()
        solved = try_current_roots()
        if solved is not None:
            return solved

        # Fallback only: a small bounded scan identifies genuinely
        # non-monotonic or disconnected intervals after boundary refinement.
        fallback_n = (
            min(max(int(grid_N), 8), 24)
            if self.root_search == "fast"
            else max(int(grid_N), 8)
        )
        for index in range(1, fallback_n):
            sample(x_edge + (x_upper - x_edge) * index / fallback_n)
        for _ in range(3):
            refine_valid_boundaries()
            solved = try_current_roots()
            if solved is not None:
                return solved

        valid_samples.sort(key=lambda item: item[0])
        if not valid_samples:
            raise OxygenKDNotReachableError(
                f"No physically valid positive-KD oxygen state was found "
                f"within the {self.root_search} search budget "
                f"({fallback_n} subdivisions)."
            )

        best = min(valid_samples, key=lambda item: abs(item[1]))
        if self.allow_nearest:
            x_value, _, state = best
            return self._result(
                kd_by_x[x_value],
                KD_O_target,
                state,
                oxygen_tol=tol,
                oxygen_converged=False,
                failure_reason="No oxygen-KD sign-changing bracket.",
            )
        achieved = [kd_by_x[item[0]] for item in valid_samples]
        raise OxygenKDNotReachableError(
            f"No oxygen-KD root was found within the {self.root_search} "
            f"search budget ({fallback_n} subdivisions). "
            f"Achievable sampled KD range is [{min(achieved):.6g}, "
            f"{max(achieved):.6g}], target={KD_O_target:.6g}."
        )


class SulfurOLSolver:
    """Strict outer root solver for oxygen loss at a target Delta IW."""

    def __init__(
        self,
        base_params: SulfurKDParams,
        KD_O_target: float,
        IW_target: float,
        solver_factory: Callable[[SulfurKDParams], ForwardKDOSolverSulfur],
        *,
        sulfur_kd_calculator: Optional[SulfurKDCalculator] = None,
        precision: Literal["fast", "normal", "high"] = "normal",
        outer_tol: Optional[float] = None,
        max_outer_iter: int = 300,
        outer_grid_N: int = 16,
        allow_nearest: bool = False,
        root_search: Literal["fast", "thorough"] = "fast",
    ) -> None:
        precision_presets = {
            "fast": {
                "inner_tol": 1e-6,
                "inner_max_iter": 100,
                "inner_grid_N": 50,
                "outer_tol": 1e-4,
            },
            "normal": {
                "inner_tol": 1e-10,
                "inner_max_iter": 200,
                "inner_grid_N": 200,
                "outer_tol": 1e-6,
            },
            "high": {
                "inner_tol": 1e-12,
                "inner_max_iter": 500,
                "inner_grid_N": 400,
                "outer_tol": 1e-8,
            },
        }
        if precision not in precision_presets:
            raise ValueError("precision must be 'fast', 'normal', or 'high'.")
        if not math.isfinite(KD_O_target) or KD_O_target <= 0.0:
            raise ValueError("KD_O_target must be finite and positive.")
        if not math.isfinite(IW_target):
            raise ValueError("IW_target must be finite.")
        if max_outer_iter < 1 or outer_grid_N < 2:
            raise ValueError("Outer iteration limit and grid size must be positive.")
        if root_search not in {"fast", "thorough"}:
            raise ValueError("root_search must be 'fast' or 'thorough'.")
        preset = precision_presets[precision]
        self.base_params = base_params
        self.KD_O_target = float(KD_O_target)
        self.IW_target = float(IW_target)
        self.solver_factory = solver_factory
        self.sulfur_kd_calculator = sulfur_kd_calculator or SulfurKDCalculator()
        self.gamma_FeO_sil = GAMMA_FEO_SILICATE
        self.inner_tol = float(preset["inner_tol"])
        self.inner_max_iter = int(preset["inner_max_iter"])
        self.inner_grid_N = int(preset["inner_grid_N"])
        self.outer_tol = (
            float(outer_tol)
            if outer_tol is not None
            else float(preset["outer_tol"])
        )
        if not math.isfinite(self.outer_tol) or self.outer_tol <= 0.0:
            raise ValueError("outer_tol must be finite and positive.")
        self.max_outer_iter = int(max_outer_iter)
        self.outer_grid_N = int(outer_grid_N)
        self.allow_nearest = bool(allow_nearest)
        self.root_search = root_search
        self._cached_result: Optional[
            Tuple[float, Dict[str, float], float]
        ] = None
        self._last_x_hint: Optional[float] = None
        self._last_inner_state: Optional[Dict[str, float]] = None

    def _compute_IW(self, res_dict: Dict[str, float]) -> float:
        p = self.base_params
        sil_sum = (
            res_dict.get("x_", 0.0)
            + res_dict.get("y_", 0.0)
            + res_dict.get("z_", 0.0)
            + p.u
            + p.m
            + p.n
            + res_dict.get("S_sil_", 0.0)
        )
        metal_sum = (
            res_dict.get("a_", 0.0)
            + res_dict.get("b_", 0.0)
            + res_dict.get("c_", 0.0)
            + res_dict.get("d_", 0.0)
            + res_dict.get("S_met_", 0.0)
        )
        if sil_sum <= 0.0 or metal_sum <= 0.0:
            raise PhysicalStateError("IW requires two positive phase totals.")
        a_feo = GAMMA_FEO_SILICATE * res_dict.get("x_", 0.0) / sil_sum
        gamma_fe = float(res_dict.get("gamma_Fe_metal", math.nan))
        if not math.isfinite(gamma_fe) or gamma_fe <= 0.0:
            ternary_sum = (
                res_dict.get("a_", 0.0)
                + res_dict.get("b_", 0.0)
                + res_dict.get("S_met_", 0.0)
            )
            if ternary_sum <= 0.0:
                raise PhysicalStateError("IW requires a positive Fe-Ni-S subsystem.")
            gamma_fe = self.sulfur_kd_calculator.get_activity_coefficients(
                X_Ni=res_dict.get("b_", 0.0) / ternary_sum,
                X_S=res_dict.get("S_met_", 0.0) / ternary_sum,
            )["Fe"]
        ternary_sum = (
            res_dict.get("a_", 0.0)
            + res_dict.get("b_", 0.0)
            + res_dict.get("S_met_", 0.0)
        )
        if ternary_sum <= 0.0:
            raise PhysicalStateError("IW requires a positive Fe-Ni-S subsystem.")
        # True Fe-Ni-S approximation: gamma_Fe and X_Fe use the same ternary
        # composition basis. Si and O remain in mass balance and partition
        # equations but do not dilute the Fe activity of this subsystem.
        a_fe = gamma_fe * res_dict.get("a_", 0.0) / ternary_sum
        if a_feo <= 0.0 or a_fe <= 0.0:
            raise PhysicalStateError("IW activities must be positive.")
        iw = 2.0 * math.log10(a_feo / a_fe)
        if not math.isfinite(iw):
            raise PhysicalStateError("Calculated Delta IW is not finite.")
        return iw

    def _evaluate_signed(
        self,
        O_L_guess: float,
    ) -> Tuple[float, Dict[str, float], float]:
        new_params = replace(self.base_params, O_L=float(O_L_guess))
        solver = self.solver_factory(new_params)
        # The outer search mode is authoritative for the complete nested solve.
        # In particular, thorough mode must not leave the inner O/Ni/S searches
        # at their default fast caps.
        solver.root_search = self.root_search
        if self._last_inner_state is not None:
            solver._last_y_hint = self._last_inner_state.get("y_")
            s_met = self._last_inner_state.get("S_met_", 0.0)
            s_sil = self._last_inner_state.get("S_sil_", 0.0)
            if s_met > 0.0 and s_sil > 0.0:
                solver._last_s_log_ratio_hint = math.log(s_met / s_sil)
        result = solver.solve_x_for_KD_O(
            self.KD_O_target,
            tol=self.inner_tol,
            max_iter=self.inner_max_iter,
            grid_N=self.inner_grid_N,
            x_hint=self._last_x_hint,
        )
        self._last_x_hint = result.x_
        res_dict = dict(vars(result))
        self._last_inner_state = res_dict
        if not result.overall_converged or not result.oxygen_converged:
            raise OxygenKDNotReachableError(
                "Outer IW evaluation received a non-converged oxygen state."
            )
        iw_model = self._compute_IW(res_dict)
        signed = iw_model - self.IW_target
        report = audit_physical_state(
            res_dict,
            new_params,
            KD_O_target=self.KD_O_target,
            IW_target=self.IW_target,
            IW_model=iw_model,
        )
        if not report.physical:
            raise PhysicalStateError("Outer IW state failed physical audit.")
        res_dict["residuals"] = dict(report.residuals)
        return signed, res_dict, iw_model

    def _residual(self, O_L_guess: float) -> float:
        """Compatibility objective: infinity for every invalid nested state."""

        try:
            signed, _, _ = self._evaluate_signed(O_L_guess)
            return abs(signed)
        except (ArithmeticError, RuntimeError, TypeError, ValueError):
            return math.inf

    def _cache(
        self,
        O_L: float,
        res_dict: Dict[str, float],
        iw_model: float,
        *,
        converged: bool,
        failure_reason: Optional[str] = None,
    ) -> None:
        res_dict = dict(res_dict)
        residuals = dict(res_dict.get("residuals", {}))
        residuals["IW"] = iw_model - self.IW_target
        res_dict.update(
            iw_converged=converged,
            overall_converged=bool(res_dict.get("overall_converged", False))
            and converged,
            failure_reason=failure_reason,
            residuals=residuals,
        )
        self._cached_result = (O_L, res_dict, iw_model)

    def solve(
        self,
        O_bounds: Optional[Tuple[float, float]] = None,
    ) -> Tuple[float, float]:
        """Solve the signed IW residual on physically feasible O_L intervals."""

        stoichiometric_oxygen = (
            self.base_params.Fe_t
            + self.base_params.Ni_t
            + 2.0 * self.base_params.Si_t
        )
        if O_bounds is None:
            lo = 0.0
            hi = stoichiometric_oxygen * (1.0 - 1e-12)
        else:
            lo, hi = map(float, O_bounds)
            if not math.isfinite(lo) or not math.isfinite(hi) or lo >= hi:
                raise ValueError("O_bounds must be two finite increasing values.")
            if hi > stoichiometric_oxygen:
                hi = stoichiometric_oxygen * (1.0 - 1e-12)
            if lo >= hi:
                raise ValueError("O_bounds contain no non-negative-oxygen states.")

        evaluation_cache: Dict[
            float,
            Optional[Tuple[float, Dict[str, float], float]],
        ] = {}

        def evaluate_cached(
            oxygen_loss: float,
        ) -> Optional[Tuple[float, Dict[str, float], float]]:
            if oxygen_loss in evaluation_cache:
                return evaluation_cache[oxygen_loss]
            try:
                evaluated = self._evaluate_signed(oxygen_loss)
            except (ArithmeticError, RuntimeError, TypeError, ValueError):
                evaluated = None
            evaluation_cache[oxygen_loss] = evaluated
            return evaluated

        def collect_valid() -> list[Tuple[float, float, Dict[str, float], float]]:
            return sorted(
                (
                    oxygen_loss,
                    evaluated[0],
                    evaluated[1],
                    evaluated[2],
                )
                for oxygen_loss, evaluated in evaluation_cache.items()
                if evaluated is not None
            )

        def classify() -> Tuple[
            list[Tuple[float, float, Dict[str, float], float]],
            list[Tuple[float, float, Dict[str, float], float]],
            list[
                Tuple[
                    Tuple[float, float, Dict[str, float], float],
                    Tuple[float, float, Dict[str, float], float],
                ]
            ],
        ]:
            ordered = sorted(evaluation_cache.items())
            samples = [
                (
                    oxygen_loss,
                    evaluated[0],
                    evaluated[1],
                    evaluated[2],
                )
                for oxygen_loss, evaluated in ordered
                if evaluated is not None
            ]
            exact_samples = [
                sample for sample in samples if abs(sample[1]) <= self.outer_tol
            ]
            bracket_samples = []
            for (left_o, left_evaluated), (right_o, right_evaluated) in zip(
                ordered,
                ordered[1:],
            ):
                if left_evaluated is None or right_evaluated is None:
                    continue
                left = (
                    left_o,
                    left_evaluated[0],
                    left_evaluated[1],
                    left_evaluated[2],
                )
                right = (
                    right_o,
                    right_evaluated[0],
                    right_evaluated[1],
                    right_evaluated[2],
                )
                if left[1] * right[1] < 0.0:
                    bracket_samples.append((left, right))
            return samples, exact_samples, bracket_samples

        def attempt_connected_bracket(
            left: Tuple[float, float, Dict[str, float], float],
            right: Tuple[float, float, Dict[str, float], float],
        ) -> Optional[Tuple[float, float, Dict[str, float], float]]:
            """Bisect a bracket, rejecting it if physical feasibility breaks."""

            lo_o, lo_f, lo_state, lo_iw = left
            hi_o, hi_f, hi_state, hi_iw = right
            if lo_f * hi_f >= 0.0:
                return None

            def objective(oxygen_loss: float) -> float:
                evaluated = evaluate_cached(oxygen_loss)
                if evaluated is None:
                    raise PhysicalStateError(
                        "IW bracket crossed an invalid physical state."
                    )
                return evaluated[0]

            try:
                root, _ = toms748(
                    objective,
                    lo_o,
                    hi_o,
                    xtol=max(
                        1e-12
                        * max(
                            abs(lo_o),
                            abs(hi_o),
                            np.finfo(float).tiny,
                        ),
                        np.finfo(float).tiny,
                    ),
                    rtol=4.0 * np.finfo(float).eps,
                    maxiter=self.max_outer_iter,
                    full_output=True,
                    disp=False,
                )
            except (ArithmeticError, RuntimeError, ValueError):
                return None
            evaluated = evaluate_cached(root)
            if evaluated is None or abs(evaluated[0]) > self.outer_tol:
                return None
            return root, evaluated[0], evaluated[1], evaluated[2]

        # Probe the supplied state first.  This is both the most likely root
        # during event-to-event continuation and the best warm start for the
        # nested oxygen solve.  Candidate brackets are tested immediately,
        # so the search stops as soon as a connected physical root is found.
        base_o = min(max(self.base_params.O_L, lo), hi)
        primary_o = [
            base_o,
            hi,
            lo,
            lo + 0.50 * (hi - lo),
            lo + 0.25 * (hi - lo),
            lo + 0.75 * (hi - lo),
        ]
        outer_subdivisions = (
            min(self.outer_grid_N, 24)
            if self.root_search == "fast"
            else self.outer_grid_N
        )
        fallback_o = sorted(
            (
                lo + (hi - lo) * index / outer_subdivisions
                for index in range(1, outer_subdivisions)
            ),
            key=lambda value: abs(value - base_o),
        )
        probe_o = list(dict.fromkeys([*primary_o, *fallback_o]))
        attempted_brackets: set[Tuple[float, float]] = set()
        probe_index = 0
        secant_probes = 0
        while probe_index < len(probe_o):
            oxygen_loss = probe_o[probe_index]
            probe_index += 1
            evaluate_cached(oxygen_loss)
            valid_samples, exact, brackets = classify()
            if exact:
                chosen = min(
                    exact,
                    key=lambda item: abs(item[0] - self.base_params.O_L),
                )
                self._cache(chosen[0], chosen[2], chosen[3], converged=True)
                return chosen[0], abs(chosen[1])

            # Two nearby valid IW evaluations provide a cheap local secant
            # predictor.  Event-to-event continuation normally makes this
            # much faster than walking a global O_L grid.  The predicted point
            # is still passed through the full strict nested solve.
            secant_inserted = False
            if len(valid_samples) >= 2 and secant_probes < 8:
                nearest = sorted(
                    valid_samples,
                    key=lambda item: abs(item[1]),
                )[:2]
                first, second = nearest
                denominator = second[1] - first[1]
                if abs(denominator) > 1e-15:
                    predicted = second[0] - second[1] * (
                        second[0] - first[0]
                    ) / denominator
                    if (
                        lo < predicted < hi
                        and predicted not in evaluation_cache
                        and predicted not in probe_o[probe_index:]
                    ):
                        probe_o.insert(probe_index, predicted)
                        secant_probes += 1
                        secant_inserted = True

            # Give the inexpensive secant prediction one strict evaluation
            # before spending nested solves on bisection.  Bisection remains
            # the guaranteed fallback when the secant stalls.
            if secant_inserted:
                continue

            for left, right in sorted(
                brackets,
                key=lambda pair: abs(
                    0.5 * (pair[0][0] + pair[1][0]) - self.base_params.O_L
                ),
            ):
                marker = (left[0], right[0])
                if marker in attempted_brackets:
                    continue
                attempted_brackets.add(marker)
                root = attempt_connected_bracket(left, right)
                if root is not None:
                    self._cache(root[0], root[2], root[3], converged=True)
                    return root[0], abs(root[1])

        # A connected feasible O_L branch can be much narrower than the
        # fallback spacing.  Refine only observed valid/invalid transitions;
        # this locates such branches without multiplying the full nested
        # solve by a dense global grid.
        ordered_cache = sorted(evaluation_cache.items())
        transitions = [
            (left, right)
            for left, right in zip(ordered_cache, ordered_cache[1:])
            if (left[1] is None) != (right[1] is None)
        ]
        def transition_priority(
            pair: Tuple[
                Tuple[
                    float,
                    Optional[Tuple[float, Dict[str, float], float]],
                ],
                Tuple[
                    float,
                    Optional[Tuple[float, Dict[str, float], float]],
                ],
            ],
        ) -> Tuple[float, float]:
            valid_evaluated = (
                pair[0][1] if pair[0][1] is not None else pair[1][1]
            )
            assert valid_evaluated is not None
            return (
                abs(valid_evaluated[0]),
                min(
                    abs(pair[0][0] - base_o),
                    abs(pair[1][0] - base_o),
                ),
            )

        transitions.sort(key=transition_priority)
        for left_item, right_item in transitions:
            invalid_o: float
            valid: Tuple[float, float, Dict[str, float], float]
            if left_item[1] is None:
                invalid_o = left_item[0]
                right_evaluated = right_item[1]
                assert right_evaluated is not None
                valid = (
                    right_item[0],
                    right_evaluated[0],
                    right_evaluated[1],
                    right_evaluated[2],
                )
            else:
                invalid_o = right_item[0]
                left_evaluated = left_item[1]
                assert left_evaluated is not None
                valid = (
                    left_item[0],
                    left_evaluated[0],
                    left_evaluated[1],
                    left_evaluated[2],
                )

            for _ in range(min(self.max_outer_iter, 32)):
                midpoint = 0.5 * (valid[0] + invalid_o)
                if midpoint == valid[0] or midpoint == invalid_o:
                    break
                evaluated = evaluate_cached(midpoint)
                if evaluated is None:
                    invalid_o = midpoint
                else:
                    midpoint_sample = (
                        midpoint,
                        evaluated[0],
                        evaluated[1],
                        evaluated[2],
                    )
                    if abs(midpoint_sample[1]) <= self.outer_tol:
                        self._cache(
                            midpoint_sample[0],
                            midpoint_sample[2],
                            midpoint_sample[3],
                            converged=True,
                        )
                        return midpoint_sample[0], abs(midpoint_sample[1])
                    if valid[1] * midpoint_sample[1] < 0.0:
                        bracket = (
                            (valid, midpoint_sample)
                            if valid[0] < midpoint_sample[0]
                            else (midpoint_sample, valid)
                        )
                        root = attempt_connected_bracket(*bracket)
                        if root is not None:
                            self._cache(
                                root[0],
                                root[2],
                                root[3],
                                converged=True,
                            )
                            return root[0], abs(root[1])
                    valid = midpoint_sample
                if abs(valid[0] - invalid_o) <= (
                    1e-10
                    * max(stoichiometric_oxygen, np.finfo(float).tiny)
                ):
                    break

        valid_samples = collect_valid()
        if not valid_samples:
            raise IWTargetNotReachableError(
                f"No physical O_L state completed the nested partition solve "
                f"within the {self.root_search} search budget "
                f"({outer_subdivisions} subdivisions)."
            )
        best = min(valid_samples, key=lambda item: abs(item[1]))
        if self.allow_nearest:
            self._cache(
                best[0],
                best[2],
                best[3],
                converged=False,
                failure_reason="No IW sign-changing bracket.",
            )
            return best[0], abs(best[1])
        iw_values = [sample[3] for sample in valid_samples]
        raise IWTargetNotReachableError(
            f"No IW root was found within the {self.root_search} search "
            f"budget ({outer_subdivisions} subdivisions). "
            f"Achievable sampled IW range is [{min(iw_values):.6g}, "
            f"{max(iw_values):.6g}], target={self.IW_target:.6g}."
        )

    def get_final_result(self) -> Tuple[Dict[str, float], float]:
        """Return the phase result exactly corresponding to the solved O_L."""

        if self._cached_result is None:
            raise RuntimeError("No cached result. Call solve() successfully first.")
        _, result, iw_model = self._cached_result
        return result, iw_model


def make_sulfur_solver_factory(
    *,
    P: float,
    T: float,
    sulfur_model: SulfurModel,
    **solver_kwargs,
) -> Callable[[SulfurKDParams], ForwardKDOSolverSulfur]:
    """Build a parameter-to-solver factory for one differentiation event."""

    def factory(params: SulfurKDParams) -> ForwardKDOSolverSulfur:
        return ForwardKDOSolverSulfur(
            params,
            P=P,
            T=T,
            sulfur_model=sulfur_model,
            **solver_kwargs,
        )

    return factory


def sulfur_model_for_event(event_index: int) -> SulfurModel:
    """Return the paper-prescribed sulfur model for an accretion event.

    Event zero is primordial differentiation and uses Rose-Weston Eq. S8.
    Every subsequent event uses the composition-dependent Boujibar Eq. S7.
    """

    if isinstance(event_index, bool) or not isinstance(event_index, int):
        raise TypeError("event_index must be an integer.")
    if event_index < 0:
        raise ValueError("event_index must be non-negative.")
    return "rose_weston" if event_index == 0 else "boujibar"


def sulfur_params_from_composition(
    composition: Mapping[str, float],
    *,
    P: float,
    T: float,
    O_L: Optional[float] = None,
    major_kd_calculator: Optional[SulfurKDCalculator] = None,
) -> SulfurKDParams:
    """Build event input parameters from a two-phase molar composition.

    This is the explicit state-transfer boundary between accretion events.
    The previous event's ``FeO``/``Fe``, ``NiO``/``Ni``, ``SiO2``/``Si`` and
    ``S_sil``/``S_met`` values are collapsed back to bulk element totals.
    A collision driver may first sum or mass-weight several predecessor
    compositions, then call this function for the new event.

    ``major_kd_calculator`` provides the base Ni and Si coefficients at the
    state-transfer boundary. Passing it explicitly lets a workflow share one
    :class:`SulfurKDCalculator`; composition-dependent activity corrections
    are not folded into these stored inputs.
    """

    def value(name: str) -> float:
        result = float(composition.get(name, 0.0))
        if not math.isfinite(result) or result < 0.0:
            raise ValueError(f"Composition value {name!r} must be finite and non-negative.")
        return result

    fe, feo = value("Fe"), value("FeO")
    ni, nio = value("Ni"), value("NiO")
    si, sio2 = value("Si"), value("SiO2")
    oxygen = value("O")
    sulfur_total = (
        value("S")
        if "S" in composition
        else value("S_sil") + value("S_met")
    )
    oxygen_loss = (
        fe + ni + 2.0 * si - oxygen
        if O_L is None
        else float(O_L)
    )
    if not math.isfinite(oxygen_loss):
        raise ValueError("O_L must be finite.")

    kd_calc = major_kd_calculator or SulfurKDCalculator()
    return SulfurKDParams(
        Fe_t=feo + fe,
        Ni_t=nio + ni,
        Si_t=sio2 + si,
        O_L=oxygen_loss,
        KD_Ni=kd_calc.get_KD("Ni", P, T),
        KD_Si=kd_calc.get_KD("Si", P, T),
        u=value("MgO"),
        m=2.0 * value("Al2O3") if "Al2O3" in composition else value("AlO1.5"),
        n=value("CaO"),
        S_t=sulfur_total,
    )


class SulfurDifferentiationEventDriver:
    """Create differentiation solvers with explicit S8-then-S7 staging.

    The class deliberately does not retain or mutate a hidden chemical state.
    Callers pass the bulk composition for each event as
    :class:`SulfurKDParams`, preserving the input structure of
    the established accretion workflow without inheriting its implementation. Use
    :func:`sulfur_params_from_composition` to convert the preceding event's
    phase result (or a collision mixture) back into those bulk parameters.
    A single externally constructed ``SulfurKDCalculator`` is shared by the
    nested Ni-S calculation and the IW activity calculation.
    """

    def __init__(
        self,
        *,
        sulfur_kd_calculator: Optional[SulfurKDCalculator] = None,
        **solver_kwargs,
    ) -> None:
        self.sulfur_kd_calculator = sulfur_kd_calculator or SulfurKDCalculator()
        self.solver_kwargs = dict(solver_kwargs)

    def make_solver_factory(
        self,
        *,
        event_index: int,
        P: float,
        T: float,
    ) -> Callable[[SulfurKDParams], ForwardKDOSolverSulfur]:
        """Return a parameter-to-solver factory for one fixed event."""

        return make_sulfur_solver_factory(
            P=P,
            T=T,
            sulfur_model=sulfur_model_for_event(event_index),
            sulfur_kd_calculator=self.sulfur_kd_calculator,
            **self.solver_kwargs,
        )

    def make_ol_solver(
        self,
        params: SulfurKDParams,
        *,
        event_index: int,
        P: float,
        T: float,
        KD_O_target: float,
        IW_target: float,
        precision: Literal["fast", "normal", "high"] = "normal",
        outer_tol: Optional[float] = None,
        max_outer_iter: int = 300,
        outer_grid_N: int = 16,
        allow_nearest: bool = False,
        root_search: Literal["fast", "thorough"] = "fast",
    ) -> SulfurOLSolver:
        """Construct the standalone outer solver for one event."""

        nested_solver_kwargs = dict(self.solver_kwargs)
        nested_solver_kwargs["root_search"] = root_search
        factory = make_sulfur_solver_factory(
            P=P,
            T=T,
            sulfur_model=sulfur_model_for_event(event_index),
            sulfur_kd_calculator=self.sulfur_kd_calculator,
            **nested_solver_kwargs,
        )
        return SulfurOLSolver(
            params,
            KD_O_target,
            IW_target,
            factory,
            sulfur_kd_calculator=self.sulfur_kd_calculator,
            precision=precision,
            outer_tol=outer_tol,
            max_outer_iter=max_outer_iter,
            outer_grid_N=outer_grid_N,
            allow_nearest=allow_nearest,
            root_search=root_search,
        )


__all__ = [
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
]
