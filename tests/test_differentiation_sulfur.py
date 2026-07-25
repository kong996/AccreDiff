import math
from dataclasses import replace

import pytest

from accrediff.differentiation import ForwardKDOSolver, KD_Params
from accrediff.chemistry import MolarMassCalculator
from accrediff.constant import Elements
from accrediff.differentiation_sulfur import (
    D_S_model_molar,
    D_S_model_rose_weston_molar,
    ForwardKDOSolverSulfur,
    GAMMA_FEO_SILICATE,
    MA2001_EPSILON_REFERENCE_T_K,
    MA2001_EPSILON_TEMPERATURE_SCALING,
    MA2001_FE_S_NI_EPSILON_1873K,
    MA2001_FE_S_NI_GAMMA0_1873K,
    Ma2001TernaryActivity,
    OxygenKDNotReachableError,
    SulfurPartitionNotReachableError,
    SulfurDifferentiationEventDriver,
    SulfurKDCalculator,
    SulfurKDParams,
    SulfurOLSolver,
    SUPPLEMENT_MAJOR_KD_PARAMS,
    supplement_major_kd,
    make_sulfur_solver_factory,
    sulfur_model_for_event,
    sulfur_params_from_composition,
    audit_physical_state,
)

MOLAR_MASS = MolarMassCalculator()
COMPONENT_MOLAR_MASS = {
    "FeO": MOLAR_MASS.molar_mass("FeO"),
    "NiO": MOLAR_MASS.molar_mass("NiO"),
    "SiO2": MOLAR_MASS.molar_mass("SiO2"),
    "MgO": MOLAR_MASS.molar_mass("MgO"),
    "AlO1.5": 0.5 * MOLAR_MASS.molar_mass("Al2O3"),
    "CaO": MOLAR_MASS.molar_mass("CaO"),
}


def base_params():
    return KD_Params(
        Fe_t=0.30,
        Ni_t=0.02,
        Si_t=0.20,
        O_L=0.10,
        KD_Ni=20.0,
        KD_Si=0.01,
        u=0.40,
        m=0.05,
        n=0.03,
    )


def sulfur_params(S_t=0.01):
    values = vars(base_params()).copy()
    values["KD_Ni"] = SulfurKDCalculator().get_KD("Ni", 20.0, 2500.0)
    return SulfurKDParams(**values, S_t=S_t)


def solver(params, model="rose_weston", **kwargs):
    return ForwardKDOSolverSulfur(
        params,
        P=20.0,
        T=2500.0,
        sulfur_model=model,
        nonneg="clip",
        enforce_z_box=True,
        enforce_d_nonneg=True,
        sulfur_grid_N=8,
        **kwargs,
    )


def test_standalone_types_keep_input_output_shape_without_inheritance():
    assert not issubclass(SulfurKDParams, KD_Params)
    assert not issubclass(ForwardKDOSolverSulfur, ForwardKDOSolver)
    assert list(vars(sulfur_params()).keys())[: len(vars(base_params()))] == list(
        vars(base_params()).keys()
    )

    actual = solver(sulfur_params(0.0)).solve_x_for_KD_O(0.01, grid_N=30)
    original_result_fields = [
        "x_",
        "a_",
        "y_",
        "b_",
        "z_",
        "c_",
        "d_",
        "KD_O_",
        "residual",
    ]
    assert list(vars(actual).keys())[: len(original_result_fields)] == original_result_fields
    assert actual.S_met_ == 0.0
    assert actual.S_sil_ == 0.0
    assert actual.sulfur_converged
    assert actual.KD_Ni_effective > 0.0
    assert actual.gamma_Fe_metal > 0.0


@pytest.mark.parametrize("model", ["rose_weston", "boujibar"])
def test_nested_sulfur_solution_conserves_s_and_matches_reported_d(model):
    params = sulfur_params()
    result = solver(params, model).solve_x_for_KD_O(0.01, grid_N=25)
    metal_sum = result.a_ + result.b_ + result.c_ + result.d_ + result.S_met_
    sil_sum = (
        result.x_
        + result.y_
        + result.z_
        + params.u
        + params.m
        + params.n
        + result.S_sil_
    )
    observed_d = (result.S_met_ / metal_sum) / (result.S_sil_ / sil_sum)

    assert result.S_met_ + result.S_sil_ == pytest.approx(params.S_t, abs=1e-15)
    assert observed_d == pytest.approx(result.D_S_mole, rel=2e-10)
    assert abs(result.sulfur_residual) <= 1e-12
    assert result.sulfur_converged
    assert result.sulfur_model == model


def test_rose_weston_function_converts_weight_d_to_molar_d():
    metal = dict(Fe=0.80, Ni=0.05, O=0.03, Si=0.10, S=0.02)
    sil = dict(FeO=0.15, NiO=0.01, SiO2=0.35, MgO=0.35, Al=0.08, CaO=0.05, S=0.01)
    d_mole = D_S_model_rose_weston_molar(
        20.0,
        2500.0,
        X_Fe_met=metal["Fe"],
        X_Ni_met=metal["Ni"],
        X_O_met=metal["O"],
        X_Si_met=metal["Si"],
        X_S_met=metal["S"],
        X_FeO_sil=sil["FeO"],
        X_NiO_sil=sil["NiO"],
        X_SiO2_sil=sil["SiO2"],
        X_MgO_sil=sil["MgO"],
        X_AlO1_5_sil=sil["Al"],
        X_CaO_sil=sil["CaO"],
        X_S_sil=sil["S"],
    )
    mw_metal = (
        metal["Fe"] * Elements["Fe"]
        + metal["Ni"] * Elements["Ni"]
        + metal["O"] * Elements["O"]
        + metal["Si"] * Elements["Si"]
        + metal["S"] * Elements["S"]
    )
    mw_sil = (
        sil["FeO"] * COMPONENT_MOLAR_MASS["FeO"]
        + sil["NiO"] * COMPONENT_MOLAR_MASS["NiO"]
        + sil["SiO2"] * COMPONENT_MOLAR_MASS["SiO2"]
        + sil["MgO"] * COMPONENT_MOLAR_MASS["MgO"]
        + sil["Al"] * COMPONENT_MOLAR_MASS["AlO1.5"]
        + sil["CaO"] * COMPONENT_MOLAR_MASS["CaO"]
        + sil["S"] * Elements["S"]
    )
    expected_weight_d = 10 ** (-4.37 + 13686 / 2500.0 + 217.49 * 20.0 / 2500.0)
    assert d_mole == pytest.approx(expected_weight_d * mw_metal / mw_sil)


def test_boujibar_function_returns_molar_coefficient():
    # The expected expression is independently evaluated from mass fractions.
    mw_metal = (
        0.80 * Elements["Fe"]
        + 0.10 * Elements["Si"]
        + 0.05 * Elements["Ni"]
        + 0.03 * Elements["O"]
        + 0.02 * Elements["S"]
    )
    kwargs = dict(
        X_FeO_oxide=0.20,
        X_CaO_oxide=0.05,
        X_MgO_oxide=0.40,
        X_FeO_sil_phase=0.18,
        X_Fe_met=0.80,
        X_Si_met=0.10,
        X_Ni_met=0.05,
        X_O_met=0.03,
        MW_metal=mw_metal,
        MW_silicate=60.0,
    )
    actual = D_S_model_molar(20.0, 2500.0, **kwargs)
    w_feo = (
        kwargs["X_FeO_sil_phase"]
        * COMPONENT_MOLAR_MASS["FeO"]
        / kwargs["MW_silicate"]
    )
    w_fe = kwargs["X_Fe_met"] * Elements["Fe"] / kwargs["MW_metal"]
    w_si = kwargs["X_Si_met"] * Elements["Si"] / kwargs["MW_metal"]
    w_ni = kwargs["X_Ni_met"] * Elements["Ni"] / kwargs["MW_metal"]
    w_o = kwargs["X_O_met"] * Elements["O"] / kwargs["MW_metal"]
    log_cs = -5.704 + 3.15 * 0.20 + 2.65 * 0.05 + 0.12 * 0.40
    lsi = math.log10(1.0 - w_si)
    log_d_weight = (
        math.log10(w_feo)
        - log_cs
        + 405 / 2500.0
        + 136 * 20.0 / 2500.0
        + 32 * lsi
        + 181 * lsi**2
        + 305 * lsi**3
        + 1.13 * math.log10(1 - w_fe)
        + 10.7 * math.log10(1 - w_ni)
        + 31.4 * math.log10(1 - w_o)
        - 3.72
    )
    expected = 10**log_d_weight * kwargs["MW_metal"] / kwargs["MW_silicate"]
    assert actual == pytest.approx(expected)


def test_boujibar_rejects_inconsistent_phase_average_molar_mass():
    with pytest.raises(ValueError, match="mass fractions cannot sum"):
        D_S_model_molar(
            20.0,
            2500.0,
            X_FeO_oxide=0.20,
            X_CaO_oxide=0.05,
            X_MgO_oxide=0.40,
            X_FeO_sil_phase=0.18,
            X_Fe_met=0.80,
            X_Si_met=0.10,
            X_Ni_met=0.05,
            X_O_met=0.03,
            MW_metal=50.0,
            MW_silicate=60.0,
        )


@pytest.mark.parametrize(
    ("x_s", "x_ni", "expected_ln_gamma", "expected_gamma"),
    [
        (
            0.10,
            0.05,
            {
                "Fe": 0.017460457307125154,
                "S": -0.45485323887595386,
                "Ni": 0.2649104471686812,
            },
            {
                "Fe": 1.017613782165915,
                "S": 0.63454108707812,
                "Ni": 1.3033142550704722,
            },
        ),
        (
            0.20,
            0.01,
            {
                "Fe": 0.12549986504167485,
                "S": -1.110196936960816,
                "Ni": 0.6163993955402624,
            },
            {
                "Fe": 1.1337150159562366,
                "S": 0.3294940651254005,
                "Ni": 1.8522468126984213,
            },
        ),
    ],
)
def test_ma2001_activity_matches_reference_notebook(
    x_s,
    x_ni,
    expected_ln_gamma,
    expected_gamma,
):
    model = Ma2001TernaryActivity(
        solvent="Fe",
        solute2="S",
        solute3="Ni",
        x2=x_s,
        x3=x_ni,
        epsilon_db=MA2001_FE_S_NI_EPSILON_1873K,
        gamma0_db=MA2001_FE_S_NI_GAMMA0_1873K,
    )
    result = model.results()

    for component in ("Fe", "S", "Ni"):
        assert result["ln_gamma"][component] == pytest.approx(
            expected_ln_gamma[component],
            rel=1e-14,
        )
        assert result["gamma"][component] == pytest.approx(
            expected_gamma[component],
            rel=1e-14,
        )


def test_sulfur_kd_calculator_separates_base_kd_and_ni_correction():
    calculator = SulfurKDCalculator()
    base = calculator.get_KD("Ni", 20.0, 2500.0)
    factor = calculator.get_ni_activity_factor(
        X_Ni=0.05,
        X_S=0.10,
    )
    assert base == pytest.approx(
        10 ** (0.35 + 2934.0 / 2500.0 - 83.0 * 20.0 / 2500.0),
        rel=1e-14,
    )
    assert base * factor == pytest.approx(5.651110960895377, rel=1e-14)
    assert calculator.get_KD("Si", 20.0, 2500.0) == pytest.approx(
        supplement_major_kd("Si", 20.0, 2500.0)
    )
    assert calculator.get_KD("O", 20.0, 2500.0) == pytest.approx(
        supplement_major_kd("O", 20.0, 2500.0)
    )
    assert MA2001_EPSILON_REFERENCE_T_K == 1873.0
    assert MA2001_EPSILON_TEMPERATURE_SCALING is False


def test_si_and_oxygen_names_remain_compatible():
    result = solver(sulfur_params()).solve_x_for_KD_O(0.01, grid_N=20)
    assert result.z_ + result.c_ == pytest.approx(sulfur_params().Si_t)
    expected_o = (
        sulfur_params().Fe_t
        + sulfur_params().Ni_t
        + 2 * sulfur_params().Si_t
        - sulfur_params().O_L
    )
    assert result.x_ + result.y_ + 2 * result.z_ + result.d_ == pytest.approx(expected_o)


def test_ni_s_feedback_is_always_active_and_reported():
    params = sulfur_params()
    result = solver(params).solve_x_for_KD_O(0.01, grid_N=20)
    ternary_sum = result.a_ + result.b_ + result.S_met_
    expected_gamma = SulfurKDCalculator().get_activity_coefficients(
        X_Ni=result.b_ / ternary_sum,
        X_S=result.S_met_ / ternary_sum,
    )
    assert result.KD_Ni_base == pytest.approx(params.KD_Ni)
    assert result.KD_Ni_effective == pytest.approx(
        result.KD_Ni_base * result.KD_Ni_activity_factor
    )
    assert result.X_Fe_ternary + result.X_Ni_ternary + result.X_S_ternary == (
        pytest.approx(1.0)
    )
    assert result.X_Fe_ternary == pytest.approx(result.a_ / ternary_sum)
    assert result.gamma_Fe_metal == pytest.approx(expected_gamma["Fe"])
    assert result.gamma_Ni_metal == pytest.approx(expected_gamma["Ni"])
    assert result.gamma_S_metal == pytest.approx(expected_gamma["S"])


def test_externally_supplied_base_ni_kd_changes_partition_solution():
    params = sulfur_params()
    low = solver(replace(params, KD_Ni=0.5 * params.KD_Ni)).forward_solve(0.20)
    high = solver(replace(params, KD_Ni=2.0 * params.KD_Ni)).forward_solve(0.20)

    assert low["KD_Ni_base"] == pytest.approx(0.5 * params.KD_Ni)
    assert high["KD_Ni_base"] == pytest.approx(2.0 * params.KD_Ni)
    assert low["y_"] > high["y_"]


def test_activity_model_uses_true_fe_ni_s_normalization_once_per_state():
    class CountingCalculator(SulfurKDCalculator):
        def __init__(self):
            super().__init__()
            self.calls = 0

        def activity_model(self, *, X_Ni, X_S):
            self.calls += 1
            return super().activity_model(X_Ni=X_Ni, X_S=X_S)

    params = sulfur_params()
    calculator = CountingCalculator()
    sulfur_solver = solver(params, sulfur_kd_calculator=calculator)
    state = sulfur_solver._state_for_y(x_=0.20, s_met=0.008, y_=0.004)
    ternary_sum = state["a_"] + state["b_"] + state["S_met_"]

    assert calculator.calls == 1
    assert state["X_Fe_ternary"] == pytest.approx(state["a_"] / ternary_sum)
    assert state["X_Ni_ternary"] == pytest.approx(state["b_"] / ternary_sum)
    assert state["X_S_ternary"] == pytest.approx(state["S_met_"] / ternary_sum)
    assert sum(
        state[name]
        for name in ("X_Fe_ternary", "X_Ni_ternary", "X_S_ternary")
    ) == pytest.approx(1.0)


def test_solver_is_invariant_to_bulk_inventory_scale():
    reference = sulfur_params()
    reference_result = solver(reference).solve_x_for_KD_O(
        0.01,
        tol=1e-10,
        grid_N=20,
    )
    amount_fields = {
        "Fe_t",
        "Ni_t",
        "Si_t",
        "O_L",
        "u",
        "m",
        "n",
        "S_t",
    }

    for scale in (1e-12, 1e-9, 1e-6, 1e3):
        scaled_values = {
            name: value * scale if name in amount_fields else value
            for name, value in vars(reference).items()
        }
        scaled_result = solver(SulfurKDParams(**scaled_values)).solve_x_for_KD_O(
            0.01,
            tol=1e-10,
            grid_N=20,
        )
        for field in ("x_", "y_", "z_", "S_met_", "S_sil_"):
            assert getattr(scaled_result, field) / scale == pytest.approx(
                getattr(reference_result, field),
                rel=2e-10,
                abs=1e-14,
            )


def test_calibration_ranges_are_documentation_not_validation():
    calculator = SulfurKDCalculator()
    assert calculator.get_KD("Ni", 0.0, 1000.0) > 0.0
    assert D_S_model_rose_weston_molar(
        0.0,
        1000.0,
        X_Fe_met=0.80,
        X_Ni_met=0.05,
        X_O_met=0.03,
        X_Si_met=0.10,
        X_S_met=0.02,
        X_FeO_sil=0.15,
        X_NiO_sil=0.01,
        X_SiO2_sil=0.35,
        X_MgO_sil=0.35,
        X_AlO1_5_sil=0.08,
        X_CaO_sil=0.05,
        X_S_sil=0.01,
    ) > 0.0


def test_default_and_explicit_sulfur_calculators_match():
    params = sulfur_params()
    default = solver(params).solve_x_for_KD_O(
        0.01,
        grid_N=20,
    )
    injected = solver(
        params,
        sulfur_kd_calculator=SulfurKDCalculator(),
    ).solve_x_for_KD_O(0.01, grid_N=20)

    for field in ("x_", "y_", "z_", "S_met_", "D_S_mole"):
        assert getattr(injected, field) == pytest.approx(
            getattr(default, field),
            rel=1e-13,
            abs=1e-15,
        )


def test_iw_can_obtain_gamma_fe_from_injected_ma_calculator():
    params = sulfur_params()
    calculator = SulfurKDCalculator()
    factory = make_sulfur_solver_factory(
        P=20.0,
        T=2500.0,
        sulfur_model="rose_weston",
        sulfur_kd_calculator=calculator,
        nonneg="clip",
        enforce_z_box=True,
        enforce_d_nonneg=True,
        sulfur_grid_N=8,
    )
    result = vars(factory(params).solve_x_for_KD_O(0.01, grid_N=20))
    outer = SulfurOLSolver(
        params,
        0.01,
        -2.0,
        factory,
        sulfur_kd_calculator=calculator,
    )
    sil_sum = (
        result["x_"]
        + result["y_"]
        + result["z_"]
        + params.u
        + params.m
        + params.n
        + result["S_sil_"]
    )
    metal_sum = (
        result["a_"]
        + result["b_"]
        + result["c_"]
        + result["d_"]
        + result["S_met_"]
    )
    ternary_sum = result["a_"] + result["b_"] + result["S_met_"]
    gamma_fe = calculator.get_activity_coefficients(
        X_Ni=result["b_"] / ternary_sum,
        X_S=result["S_met_"] / ternary_sum,
    )["Fe"]
    expected = 2.0 * math.log10(
        (GAMMA_FEO_SILICATE * result["x_"] / sil_sum)
        / (gamma_fe * result["a_"] / ternary_sum)
    )
    assert GAMMA_FEO_SILICATE == 1.7
    assert outer.gamma_FeO_sil == 1.7
    assert outer._compute_IW(result) == pytest.approx(expected, rel=1e-14)


def test_outer_solver_recovers_known_connected_iw_root():
    params = sulfur_params()
    factory = make_sulfur_solver_factory(
        P=20.0,
        T=2500.0,
        sulfur_model="rose_weston",
        sulfur_grid_N=8,
    )
    # Independently generated from the strict inner solution at O_L = 0.08.
    iw_target = -0.4721086233783975
    outer = SulfurOLSolver(
        params,
        0.01,
        iw_target,
        factory,
        precision="fast",
        outer_grid_N=8,
    )
    oxygen_loss, iw_error = outer.solve((0.06, 0.11))
    result, iw_model = outer.get_final_result()
    report = audit_physical_state(
        result,
        replace(params, O_L=oxygen_loss),
        KD_O_target=0.01,
        IW_target=iw_target,
        IW_model=iw_model,
    )

    assert oxygen_loss == pytest.approx(0.08, abs=5e-5)
    assert iw_error < 1e-4
    assert report.max_abs < 1e-4
    assert result["overall_converged"]


def test_event_driver_uses_s8_then_s7_without_hidden_switching():
    driver = SulfurDifferentiationEventDriver(
        nonneg="clip",
        enforce_z_box=True,
        enforce_d_nonneg=True,
        sulfur_grid_N=8,
    )
    primordial = driver.make_solver_factory(event_index=0, P=20.0, T=2500.0)(
        sulfur_params()
    )
    later = driver.make_solver_factory(event_index=1, P=20.0, T=2500.0)(
        sulfur_params()
    )

    assert sulfur_model_for_event(0) == "rose_weston"
    assert sulfur_model_for_event(1) == "boujibar"
    assert sulfur_model_for_event(20) == "boujibar"
    assert primordial.sulfur_model == "rose_weston"
    assert later.sulfur_model == "boujibar"
    assert primordial.forward_solve(0.20)["sulfur_model"] == "rose_weston"
    assert later.forward_solve(0.20)["sulfur_model"] == "boujibar"


def test_previous_phase_result_can_be_rebuilt_as_next_event_bulk_params():
    previous = solver(sulfur_params()).solve_x_for_KD_O(0.01, grid_N=20)
    composition = {
        "FeO": previous.x_,
        "NiO": previous.y_,
        "SiO2": previous.z_,
        "Fe": previous.a_,
        "Ni": previous.b_,
        "Si": previous.c_,
        "O": previous.d_,
        "S_sil": previous.S_sil_,
        "S_met": previous.S_met_,
        "MgO": sulfur_params().u,
        "AlO1.5": sulfur_params().m,
        "CaO": sulfur_params().n,
    }
    rebuilt = sulfur_params_from_composition(
        composition,
        P=20.0,
        T=2500.0,
        O_L=sulfur_params().O_L,
    )

    assert rebuilt.Fe_t == pytest.approx(sulfur_params().Fe_t)
    assert rebuilt.Ni_t == pytest.approx(sulfur_params().Ni_t)
    assert rebuilt.Si_t == pytest.approx(sulfur_params().Si_t)
    assert rebuilt.S_t == pytest.approx(sulfur_params().S_t)
    assert rebuilt.u == pytest.approx(sulfur_params().u)
    assert rebuilt.m == pytest.approx(sulfur_params().m)
    assert rebuilt.n == pytest.approx(sulfur_params().n)


def test_low_pressure_state_transfer_uses_supplement_si_parameters():
    composition = {
        "FeO": 0.20,
        "NiO": 0.01,
        "SiO2": 0.30,
        "Fe": 0.10,
        "Ni": 0.01,
        "Si": 0.02,
        "O": 0.001,
        "S": 0.005,
        "MgO": 0.30,
        "Al2O3": 0.04,
        "CaO": 0.02,
    }
    params = sulfur_params_from_composition(composition, P=0.1, T=2000.0)
    a, b, c = SUPPLEMENT_MAJOR_KD_PARAMS["Si"]
    assert params.KD_Si == pytest.approx(10 ** (a + b / 2000.0 + c * 0.1 / 2000.0))
    assert supplement_major_kd("O", 0.1, 2000.0) == pytest.approx(
        10 ** (0.6 - 3800.0 / 2000.0 + 22.0 * 0.1 / 2000.0)
    )


def test_state_transfer_accepts_an_external_major_kd_calculator():
    class StubKDCalculator:
        def get_KD(self, element, P, T):
            assert (P, T) == (0.1, 2000.0)
            return {"Ni": 123.0, "Si": 0.456}[element]

    composition = {
        "FeO": 0.20,
        "NiO": 0.01,
        "SiO2": 0.30,
        "Fe": 0.10,
        "Ni": 0.01,
        "Si": 0.02,
        "O": 0.001,
        "S": 0.005,
        "MgO": 0.30,
        "Al2O3": 0.04,
        "CaO": 0.02,
    }
    params = sulfur_params_from_composition(
        composition,
        P=0.1,
        T=2000.0,
        major_kd_calculator=StubKDCalculator(),
    )
    assert params.KD_Ni == 123.0
    assert params.KD_Si == 0.456


def test_oxygen_exchange_uses_fixed_fischer_definition():
    params = sulfur_params()
    sulfur_solver = solver(params)
    kd_o, result = sulfur_solver.compute_KD_O_from_x(0.2)
    sil_sum = (
        result["x_"]
        + result["y_"]
        + result["z_"]
        + params.u
        + params.m
        + params.n
        + result["S_sil_"]
    )
    metal_sum = (
        result["a_"]
        + result["b_"]
        + result["c_"]
        + result["d_"]
        + result["S_met_"]
    )
    expected = (
        result["a_"]
        * result["d_"]
        / ((result["x_"] / sil_sum) * metal_sum**2)
    )
    assert kd_o == pytest.approx(expected)


def test_strict_solver_rejects_nonconverged_coupled_state():
    params = SulfurKDParams(
        Fe_t=0.3,
        Ni_t=0.08,
        Si_t=0.2,
        O_L=0.1,
        KD_Ni=20.0,
        KD_Si=0.01,
        u=0.4,
        m=0.05,
        n=0.03,
        S_t=0.15,
    )
    with pytest.raises(SulfurPartitionNotReachableError):
        ForwardKDOSolverSulfur(
            params,
            P=20.0,
            T=2500.0,
            sulfur_model="boujibar",
            coupled_max_iter=1,
            coupled_tol=1e-12,
            enforce_z_box=True,
            enforce_d_nonneg=True,
        ).forward_solve(0.2)


def test_diagnostic_mode_returns_flagged_nonconverged_state():
    params = SulfurKDParams(
        Fe_t=0.3,
        Ni_t=0.08,
        Si_t=0.2,
        O_L=0.1,
        KD_Ni=20.0,
        KD_Si=0.01,
        u=0.4,
        m=0.05,
        n=0.03,
        S_t=0.15,
    )
    result = ForwardKDOSolverSulfur(
        params,
        P=20.0,
        T=2500.0,
        sulfur_model="boujibar",
        coupled_max_iter=1,
        sulfur_max_iter=1,
        coupled_tol=1e-12,
        physical_mode="diagnostic",
    ).forward_solve(0.2)
    assert not result["coupled_converged"]
    assert not result["overall_converged"]
    assert result["failure_reason"]


@pytest.mark.parametrize("model", ["rose_weston", "boujibar"])
@pytest.mark.parametrize("sulfur_total", [0.0, 0.005, 0.02])
def test_forward_states_pass_independent_physical_audit(model, sulfur_total):
    params = sulfur_params(sulfur_total)
    result = solver(params, model=model).forward_solve(0.2)
    report = audit_physical_state(result, params)
    assert report.physical
    assert report.max_abs < 1e-10
    assert result["overall_converged"]


def test_outer_solver_rejects_nonconverged_sulfur_state():
    params = sulfur_params()
    factory = make_sulfur_solver_factory(
        P=20.0,
        T=2500.0,
        sulfur_model="rose_weston",
        sulfur_max_iter=1,
        sulfur_tol=1e-30,
    )
    outer = SulfurOLSolver(params, 0.01, -2.0, factory, precision="fast")
    assert math.isinf(outer._residual(params.O_L))
    with pytest.raises(RuntimeError, match="No cached result"):
        outer.get_final_result()


def test_disappearing_metal_endpoint_is_reported_as_unreachable():
    params = SulfurKDParams(
        Fe_t=1.0,
        Ni_t=0.0,
        Si_t=0.0,
        O_L=1.0,
        KD_Ni=1.0,
        KD_Si=1.0,
        u=1.0,
        S_t=0.1,
    )
    with pytest.raises(OxygenKDNotReachableError, match="No positive FeO interval"):
        ForwardKDOSolverSulfur(
            params,
            P=1.0,
            T=2000.0,
            sulfur_model="rose_weston",
            enforce_d_nonneg=True,
        ).solve_x_for_KD_O(0.01, grid_N=10)


def test_extreme_rose_weston_partition_uses_log_sulfur_coordinate():
    params = SulfurKDParams(
        Fe_t=0.30,
        Ni_t=0.02,
        Si_t=0.20,
        O_L=0.10,
        KD_Ni=7.0,
        KD_Si=0.01,
        u=0.40,
        m=0.05,
        n=0.03,
        S_t=0.01,
    )
    result = ForwardKDOSolverSulfur(
        params,
        P=1000.0,
        T=1000.0,
        sulfur_model="rose_weston",
        sulfur_grid_N=20,
    ).forward_solve(0.20)

    assert 0.0 < result["S_sil_"] < params.S_t * 1e-200
    assert result["S_met_"] + result["S_sil_"] == pytest.approx(
        params.S_t,
        rel=1e-15,
    )
    assert abs(result["sulfur_residual"]) <= 1e-12
    assert result["sulfur_iterations"] < 10


def test_finite_sulfur_root_beyond_previous_fixed_q_limit_is_recovered():
    params = SulfurKDParams(
        Fe_t=0.30,
        Ni_t=0.02,
        Si_t=0.20,
        O_L=0.10,
        KD_Ni=7.0,
        KD_Si=0.01,
        u=0.40,
        m=0.05,
        n=0.03,
        S_t=0.01,
    )
    result = ForwardKDOSolverSulfur(
        params,
        P=1362.0,
        T=1000.0,
        sulfur_model="rose_weston",
        sulfur_grid_N=200,
        root_search="thorough",
    ).forward_solve(0.20)
    q_value = math.log(result["S_met_"] / result["S_sil_"])

    # The physical root must remain reachable beyond the former hard-coded
    # q=700 ceiling.  Do not tie this regression test to a duplicate atomic
    # mass table: the exact q shifts slightly with the canonical Elements data.
    assert q_value > 700.0
    assert abs(result["sulfur_residual"]) <= 1e-12
    assert q_value > 700.0
    assert math.isfinite(result["D_S_mole"])
    assert result["S_sil_"] > 0.0
    assert abs(result["sulfur_residual"]) <= 1e-12


def test_warm_start_cannot_change_physical_acceptance_of_same_state():
    params = SulfurKDParams(
        Fe_t=0.46768124266489997,
        Ni_t=0.027771645230251153,
        Si_t=0.13871567379671562,
        O_L=0.030811329295888472,
        KD_Ni=44.05941012885012,
        KD_Si=2.287838579062453e-7,
        u=0.526168238482919,
        m=0.060367940467362574,
        n=0.0391300650178483,
        S_t=0.07968593744228278,
    )
    kwargs = dict(
        P=8.842982041790137,
        T=1700.1288501764031,
        sulfur_model="rose_weston",
        sulfur_grid_N=16,
        ni_grid_N=20,
    )
    target_x = 0.1479847832988515
    cold = ForwardKDOSolverSulfur(params, **kwargs).forward_solve(target_x)

    warm_solver = ForwardKDOSolverSulfur(params, **kwargs)
    neighbor = warm_solver.forward_solve(0.11692031066622499)
    assert neighbor["overall_converged"]
    warm = warm_solver.forward_solve(target_x)

    for field in ("y_", "z_", "S_met_", "S_sil_", "D_S_mole"):
        assert warm[field] == pytest.approx(cold[field], rel=1e-10, abs=1e-15)
    assert warm["overall_converged"]
    assert abs(warm["residuals"]["Si_equilibrium"]) <= 1e-10


def test_thorough_connected_search_uses_full_subdivision_budget():
    params = sulfur_params()

    def artificial_sample(coordinate):
        if not 0.49 <= coordinate <= 0.51:
            return None
        state = {"coordinate": coordinate}
        return coordinate, coordinate - 0.503, state

    fast = solver(params, root_search="fast")
    exact, brackets, valid = fast._discover_connected_roots(
        artificial_sample,
        0.0,
        1.0,
        tolerance=1e-12,
        subdivisions=20,
    )
    assert not exact
    assert not brackets
    assert valid

    thorough = solver(params, root_search="thorough")
    exact, brackets, valid = thorough._discover_connected_roots(
        artificial_sample,
        0.0,
        1.0,
        tolerance=1e-12,
        subdivisions=1000,
    )
    assert exact or brackets
    assert len(valid) > 10


def test_search_mode_uses_and_propagates_actual_budget():
    params = sulfur_params()
    fast = solver(params, root_search="fast", ni_grid_N=1000)
    thorough = solver(params, root_search="thorough", ni_grid_N=1000)
    assert fast._effective_search_subdivisions(1000) == 64
    assert thorough._effective_search_subdivisions(1000) == 1000

    driver = SulfurDifferentiationEventDriver(sulfur_grid_N=100)
    outer = driver.make_ol_solver(
        params,
        event_index=0,
        P=20.0,
        T=2500.0,
        KD_O_target=0.01,
        IW_target=-0.5,
        root_search="thorough",
        outer_grid_N=100,
    )
    assert outer.root_search == "thorough"
    assert outer.solver_factory(params).root_search == "thorough"


def test_unreachable_oxygen_error_reports_search_budget():
    with pytest.raises(OxygenKDNotReachableError, match="search budget"):
        solver(sulfur_params(), root_search="fast").solve_x_for_KD_O(
            1e-20,
            grid_N=200,
        )
