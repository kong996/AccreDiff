# Sulfur differentiation refactor: 15-step checklist

This checklist reconstructs the agreed refactor from the scientific and
interface requirements recorded during the `ForwardKDSolver_sulfur_V2`
review.  The original implementations remain independent references.

| Step | Requirement | Status | Implementation evidence |
|---:|---|---|---|
| 1 | Preserve `differentiation.py`, `CoreMantleSolver`, and the original split solvers as independent references. | Complete | No changes to `src/accrediff/differentiation.py` or `trace_test/definition.py`. |
| 2 | Put the production sulfur calculation in an independent module. | Complete | `src/accrediff/differentiation_sulfur.py` does not import or inherit solver/data classes from `differentiation.py`. |
| 3 | Preserve the original input/output field layout without inheritance and without changing the meanings of `x_`–`d_`. | Complete | Standalone `SulfurKDParams` and `SulfurKDResult`; `c_` remains metallic Si and `d_` metallic O. |
| 4 | Keep the external sulfur calculation on a molar basis. | Complete | Both public `D_S` functions accept mole fractions and return `D_S_mole`. |
| 5 | Implement Rose-Weston Eq. 16 / Brennan Eq. S8 for primordial differentiation. | Complete | `D_S_model_rose_weston_molar`. |
| 6 | Implement Boujibar Eqs. 6 and 11 / Brennan Eq. S7 for later events. | Complete | `D_S_model_molar`. |
| 7 | Exclude sulfur from the Eq. 6 oxide normalization. | Complete | `oxide_sum` contains only oxide components. |
| 8 | Record deferred TiO2, Na2O, and K2O sulfide-capacity terms. | Complete | `TODO(science)` beside Eq. 6. |
| 9 | Convert Eq. 11 FeO and metal compositions to mass fractions internally, then convert `D_S_weight` back to `D_S_mole`. | Complete | `D_S_model_molar`. |
| 10 | Solve sulfur allocation by an inner root iteration for every trial major-element state. | Complete | `ForwardKDOSolverSulfur.forward_solve`. |
| 11 | Include sulfur in metal/silicate phase totals and in the Si and O equilibrium definitions. | Complete | `_f_z_with_sulfur` and `compute_KD_O_from_x`. |
| 12 | Make the Fe-S-Ni activity correction an intrinsic part of the standalone model. | Complete | `SulfurKDCalculator` is always active; the redundant compatibility switch and function wrapper were removed. |
| 13 | Defer the sulfur correction to oxygen partitioning and use the supplementary activity treatment. | Complete | `gamma_FeO_silicate=1.7` is fixed from Supplementary Methods S1; metal Fe/Ni/S activities use Ma (2001); 1873 K epsilon values are explicitly recorded as unscaled. |
| 14 | Use S8 for the primordial event and S7 for subsequent accretion events, with an explicit event-stage driver and no hidden model switch inside one solve. | Complete | `SulfurDifferentiationEventDriver`, `sulfur_model_for_event`, and `sulfur_params_from_composition`. |
| 15 | Expose the standalone API and provide detailed reference comparisons, convergence diagnostics, and regression tests. | Complete | Public exports, `tests/test_differentiation_sulfur.py`, and the regression calculations in `trace_test/comparison.ipynb`. |

## Event-level rule

The event driver, rather than the inner nonlinear solver, owns model
selection:

```text
event 0 (primordial differentiation) -> Rose-Weston / S8
event 1..N (subsequent accretion)     -> Boujibar / S7
```

Within any one event, the selected sulfur model remains fixed while the
major-element and sulfur allocations iterate to convergence.
