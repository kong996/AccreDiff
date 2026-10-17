"""Reusable gas-drag/Type-I time-scale model.

Structure mirrors AccreDiff's gas_migration.py: configuration, physical
models, comparison layer, diagnostics, and a factory function.
"""
from dataclasses import dataclass
import numpy as np

from .gas_migration import Gas_ModelConfig, Gas_build_model

__all__ = [
    "ComparisonConfig",
    "TypeIAdapter",
    "AerodynamicDragModel",
    "TimescaleComparison",
    "build_comparison",
]

YEAR_SECONDS = 365.25*86400.0


# =========================================================
# 1. Configuration
# =========================================================
@dataclass
class ComparisonConfig:
    radius_au: float = 2.0
    time_myr: float = 0.0
    mass_min_mearth: float = 1e-8
    mass_max_mearth: float = 1e-1
    n_mass: int = 500
    target_mass_mearth: float = 1e-4
    eccentricities: tuple = (0.001, 0.01, 0.1)
    inclination_over_e: float = 0.5
    solid_density_cgs: float = 2.0
    drag_coefficient: float = 0.5
    sigma0_cgs: float = 1700.0
    tau_decay_myr: float = 1.0
    a_alpha1: float = 1.0
    a_alpha2: float = 0.5
    beta: float = 0.25
    h1: float = 0.0335743
    gamma_eff: float = 1.4
    alpha_visc: float = 1e-4
    prandtl: float = 1.0
    use_saturation: bool = True

    @property
    def mass_grid_mearth(self):
        return np.logspace(np.log10(self.mass_min_mearth),
                           np.log10(self.mass_max_mearth), self.n_mass)


# =========================================================
# 2. Shared AccreDiff model factory
# =========================================================
class TypeIAdapter:
    def __init__(self, cfg):
        self.cfg = cfg

    def build_model(self, mass_mearth):
        p = self.cfg
        gas_cfg = Gas_ModelConfig(
            particle_mass_mearth=mass_mearth, a0=p.radius_au,
            Sigma0_cgs=p.sigma0_cgs, tau_decay_myr=p.tau_decay_myr,
            a_alpha1=p.a_alpha1, a_alpha2=p.a_alpha2, beta=p.beta,
            def_h_1=p.h1, gamma_eff=p.gamma_eff,
            alpha_visc=p.alpha_visc, prandtl=p.prandtl,
            use_saturation=p.use_saturation,
        )
        return Gas_build_model(gas_cfg)

    def rates(self, mass_mearth, eccentricities):
        model = self.build_model(mass_mearth)
        units, disk = model["units"], model["disk"]
        integrator = model["integrator"]
        t_code = float(units.myr_to_code(self.cfg.time_myr))
        a = np.array([self.cfg.radius_au])
        adot_torque = float(integrator.da_dt(t_code, a)[0])
        adot_total = {}
        for e in eccentricities:
            inc = self.cfg.inclination_over_e*e
            adot_total[e] = float(integrator.da_dt_eccentric(
                t_code, a, np.array([e]), np.array([inc]))[0])
        return model, adot_torque, adot_total


# =========================================================
# 3. Aerodynamic drag, using the same local disk
# =========================================================
class AerodynamicDragModel:
    def __init__(self, cfg):
        self.cfg = cfg

    def semimajor_axis_timescale_year(self, mass_mearth, e, disk, t_code):
        p, a = self.cfg, self.cfg.radius_au
        inc = p.inclination_over_e*e
        gas_cfg = disk.cfg
        mass_g = mass_mearth*gas_cfg.Mearth_to_Msun*gas_cfg.Msun_in_g
        radius_cm = (3*mass_g/(4*np.pi*p.solid_density_cgs))**(1/3)
        h = float(disk.aspect_ratio(a))
        sigma = float(disk.sigma_cgs(a, t_code))
        rho_g = sigma/(np.sqrt(2*np.pi)*h*a*gas_cfg.AU_in_cm)
        code_time_seconds = YEAR_SECONDS/(2*np.pi)
        vk = float(disk.omega(a))*gas_cfg.AU_in_cm/code_time_seconds
        tau0 = (2*mass_g/(np.pi*p.drag_coefficient*radius_cm**2*rho_g*vk)
                / YEAR_SECONDS)
        pressure_slope = float(disk.sigma_slope_p(a)) + 2.0 - p.beta
        eta = 0.5*h**2*pressure_slope
        qrel = np.sqrt(5*e**2/8 + inc**2/2 + eta**2)
        return tau0/(2*qrel*(eta + 17*e**2/16 + inc**2/8))


# =========================================================
# 4. Time-mass comparison
# =========================================================
class TimescaleComparison:
    def __init__(self, cfg, type1, aerodynamic):
        self.cfg, self.type1, self.aerodynamic = cfg, type1, aerodynamic

    @staticmethod
    def _timescale_year(a, adot, year_per_code):
        return np.inf if adot == 0 else a/abs(adot)*year_per_code

    def profile(self):
        p, masses = self.cfg, self.cfg.mass_grid_mearth
        gas = {e: [] for e in p.eccentricities}
        total = {e: [] for e in p.eccentricities}
        total_sign = {e: [] for e in p.eccentricities}
        torque, torque_sign = [], []
        for mass in masses:
            model, adot_torque, adot_total = self.type1.rates(mass, p.eccentricities)
            units, disk = model["units"], model["disk"]
            t_code = float(units.myr_to_code(p.time_myr))
            torque.append(self._timescale_year(p.radius_au, adot_torque, units.YEAR_PER_CODE))
            torque_sign.append(np.sign(adot_torque))
            for e in p.eccentricities:
                gas[e].append(self.aerodynamic.semimajor_axis_timescale_year(
                    mass, e, disk, t_code))
                total[e].append(self._timescale_year(
                    p.radius_au, adot_total[e], units.YEAR_PER_CODE))
                total_sign[e].append(np.sign(adot_total[e]))
        return {
            "mass_mearth": masses,
            "gas_ta_year": {e: np.asarray(v) for e, v in gas.items()},
            "type1_torque_ta_year": np.asarray(torque),
            "type1_torque_sign": np.asarray(torque_sign),
            "cn08_total_ta_year": {e: np.asarray(v) for e, v in total.items()},
            "cn08_total_sign": {e: np.asarray(v) for e, v in total_sign.items()},
        }

    @staticmethod
    def crossings(mass, first, second):
        z = np.log(first/second)
        ok = np.isfinite(z[:-1]) & np.isfinite(z[1:]) & (z[:-1]*z[1:] <= 0)
        roots = []
        for j in np.where(ok)[0]:
            roots.append(np.exp(np.log(mass[j])-z[j]*np.log(mass[j+1]/mass[j])
                                /(z[j+1]-z[j])))
        return roots


# =========================================================
# 5. Public factory
# =========================================================
def build_comparison(cfg=None):
    cfg = cfg or ComparisonConfig()
    type1 = TypeIAdapter(cfg)
    aerodynamic = AerodynamicDragModel(cfg)
    comparison = TimescaleComparison(cfg, type1, aerodynamic)
    return {"config": cfg, "type1": type1, "aerodynamic": aerodynamic,
            "comparison": comparison}
