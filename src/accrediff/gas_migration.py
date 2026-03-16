"""
gas_migration.py
────────────────
行星气盘迁移计算模块，基于 Cresswell & Nelson (2008)。

用法示例
--------
from wisdom.gas_migration import ModelConfig, build_model

cfg = ModelConfig(particle_mass_mearth=0.1, a0=1.0, tau_decay_myr=1.0)
model = build_model(cfg)
result = model["integrator"].integrate_full(e0=0.05)
"""

from dataclasses import dataclass
import numpy as np # type ignore


# =========================================================
# 1. ModelConfig
# =========================================================

@dataclass
class Gas_ModelConfig:
    # ---- GENGA-like code units ----
    G: float = 1.0
    Mstar: float = 1.0
    Mearth_to_Msun: float = 3.003489e-6

    # ---- particle ----
    particle_mass_mearth: float = 0.1
    a0: float = 0.5             # AU

    # ---- disk structure ----
    a_alpha1: float = 1.0       # alpha(r) = a_alpha1 * ln(r) + a_alpha2
    a_alpha2: float = 0.0
    beta: float = 0.25          # H/r = def_h_1 * r^beta
    gamma_eff: float = 1.4
    def_h_1: float = 0.0335743

    # ---- gas surface density (Msun/AU^2) ----
    #Sigma0_code: float = 1.800743117510361e-4
    AU_in_cm: float = 1.495978707e13
    Msun_in_g: float = 1.98847e33
    Sigma0_cgs: float = 1600.0  # g/cm^2
    Sigma0_code = Sigma0_cgs * AU_in_cm**2 / Msun_in_g

    # ---- gas disk decay ----
    tau_decay_myr: float = 2.0

    # ---- integration ----
    t_max_myr: float = 3.0
    n_steps: int = 20000
    a_min_stop: float = 0.05
    a_max_stop: float = 10.0


# =========================================================
# 2. UnitSystem
# =========================================================

class Gas_UnitSystem:
    """GENGA-like 单位换算（G=1, Mstar=1, AU=1）。"""

    def __init__(self):
        self.YEAR_PER_CODE = 1.0 / (2.0 * np.pi)
        self.MYR_PER_CODE  = self.YEAR_PER_CODE / 1.0e6

    def code_to_year(self, t_code):
        return np.asarray(t_code) * self.YEAR_PER_CODE

    def code_to_myr(self, t_code):
        return np.asarray(t_code) * self.MYR_PER_CODE

    def myr_to_code(self, t_myr):
        return np.asarray(t_myr) / self.MYR_PER_CODE


# =========================================================
# 3. DiskModel
# =========================================================

class Gas_DiskModel:
    """
    气盘结构模型：
      Sigma(r,t) = Sigma0 * exp[-(a_alpha1*ln(r)+a_alpha2)*ln(r)] * exp(-t/tau)
      H/r = def_h_1 * r^beta
    """

    def __init__(self, cfg: Gas_ModelConfig, units: Gas_UnitSystem | None = None):
        self.cfg   = cfg
        self.units = units or Gas_UnitSystem()
        self.tau_decay_code = self.units.myr_to_code(cfg.tau_decay_myr)

    def sigma_slope_p(self, r):
        """p(r) = -dlnSigma/dlnr"""
        r = np.asarray(r)
        return 2.0 * self.cfg.a_alpha1 * np.log(r) + self.cfg.a_alpha2

    def temp_slope_q(self, r):
        """q = -dlnT/dlnr = 1 - 2*beta"""
        r = np.asarray(r)
        return np.full_like(r, 1.0 - 2.0 * self.cfg.beta, dtype=float)

    def entropy_slope_xi(self, r):
        p = self.sigma_slope_p(r)
        q = self.temp_slope_q(r)
        return q - (self.cfg.gamma_eff - 1.0) * p

    def aspect_ratio(self, r):
        """H/r"""
        r = np.asarray(r)
        return self.cfg.def_h_1 * r ** self.cfg.beta

    def scale_height(self, r):
        """H(r)"""
        r = np.asarray(r)
        return self.cfg.def_h_1 * r * r ** self.cfg.beta

    def sigma(self, r, t_code):
        r      = np.asarray(r)
        t_code = np.asarray(t_code)
        return (
            self.cfg.Sigma0_code
            * np.exp(-(self.cfg.a_alpha1 * np.log(r) + self.cfg.a_alpha2) * np.log(r))
            * np.exp(-t_code / self.tau_decay_code)
        )

    def omega(self, r):
        r = np.asarray(r)
        return np.sqrt(self.cfg.G * self.cfg.Mstar / r ** 3)


# =========================================================
# 4. TorqueModel
# =========================================================

class Gas_TorqueModel:
    """
    未饱和静态力矩模型（Paardekooper et al. 2010/2011）：
      Gamma_tot = Gamma_L + Gamma_C_baro + Gamma_C_ent
    """

    def __init__(self, cfg: Gas_ModelConfig, disk: Gas_DiskModel):
        self.cfg  = cfg
        self.disk = disk
        self.m_particle = cfg.particle_mass_mearth * cfg.Mearth_to_Msun

    @property
    def qplanet(self):
        return self.m_particle / self.cfg.Mstar

    def gamma0(self, r, t_code):
        r      = np.asarray(r)
        t_code = np.asarray(t_code)
        h      = self.disk.aspect_ratio(r)
        sigma  = self.disk.sigma(r, t_code)
        omega  = self.disk.omega(r)
        return (self.qplanet / h) ** 2 * sigma * r ** 4 * omega ** 2

    def gamma_l(self, r, t_code):
        p = self.disk.sigma_slope_p(r)
        q = self.disk.temp_slope_q(r)
        return (-2.5 - 1.7 * q + 0.1 * p) / self.cfg.gamma_eff * self.gamma0(r, t_code)

    def gamma_c_baro(self, r, t_code):
        p = self.disk.sigma_slope_p(r)
        return 1.1 * (1.5 - p) / self.cfg.gamma_eff * self.gamma0(r, t_code)

    def gamma_c_ent(self, r, t_code):
        xi = self.disk.entropy_slope_xi(r)
        return 7.9 * xi / self.cfg.gamma_eff ** 2 * self.gamma0(r, t_code)

    def gamma_tot(self, r, t_code):
        return self.gamma_l(r, t_code) + self.gamma_c_baro(r, t_code) + self.gamma_c_ent(r, t_code)

    def profile(self, r_array, t_code=0.0):
        r_array = np.asarray(r_array)
        return {
            "r":            r_array,
            "gamma_l":      self.gamma_l(r_array, t_code),
            "gamma_c_baro": self.gamma_c_baro(r_array, t_code),
            "gamma_c_ent":  self.gamma_c_ent(r_array, t_code),
            "gamma_tot":    self.gamma_tot(r_array, t_code),
        }


# =========================================================
# 5. MigrationIntegrator — Cresswell & Nelson (2008)
# =========================================================

class Gas_MigrationIntegrator:
    """
    联立积分 a(t), e(t), inc(t)。
    阻尼公式来自 Cresswell & Nelson (2008) [CN08]。
    """

    def __init__(self, cfg: Gas_ModelConfig, torque_model: Gas_TorqueModel, units: Gas_UnitSystem | None = None):
        self.cfg    = cfg
        self.torque = torque_model
        self.units  = units or Gas_UnitSystem()
        self.m_particle = cfg.particle_mass_mearth * cfg.Mearth_to_Msun

    # ---- 基础力矩驱动 ----

    def da_dt(self, t_code, a):
        a    = np.asarray(a, dtype=float)
        out  = np.zeros_like(a)
        mask = a > 0.0
        out[mask] = (
            2.0 * self.torque.gamma_tot(a[mask], t_code) * np.sqrt(a[mask])
            / (self.m_particle * np.sqrt(self.cfg.G * self.cfg.Mstar))
        )
        return float(out) if (np.isscalar(a) or a.shape == ()) else out

    def migration_timescale_myr(self, t_code, a):
        a    = np.asarray(a, dtype=float)
        adot = np.asarray(self.da_dt(t_code, a), dtype=float)
        out  = np.full_like(a, np.inf)
        mask = np.abs(adot) > 0.0
        out[mask] = a[mask] / np.abs(adot[mask])
        out = out * self.units.MYR_PER_CODE
        return float(out) if (np.isscalar(a) or a.shape == ()) else out

    # ---- CN08 eq.7：tau_wave ----

    def tau_wave(self, t_code, a):
        a     = np.asarray(a, dtype=float)
        h     = self.torque.disk.aspect_ratio(a)
        sigma = self.torque.disk.sigma(a, t_code)
        omega = self.torque.disk.omega(a)
        q_p   = self.m_particle / self.cfg.Mstar
        return (1.0 / q_p) * (self.cfg.Mstar / (sigma * a ** 2)) * h ** 4 / omega

    # ---- CN08 eq.11：tau_e ----

    def eccentricity_damping_timescale(self, t_code, a, e, inc):
        a, e, inc = np.asarray(a, dtype=float), np.asarray(e, dtype=float), np.asarray(inc, dtype=float)
        h  = self.torque.disk.aspect_ratio(a)
        tw = self.tau_wave(t_code, a)
        eh, ih = e / h, inc / h
        bracket = np.maximum(1.0 - 0.14 * eh**2 + 0.06 * eh**3 + 0.18 * eh * ih**2, 1e-10)
        return tw / (0.780 * bracket)

    # ---- CN08 eq.12：tau_i ----

    def inclination_damping_timescale(self, t_code, a, e, inc):
        a, e, inc = np.asarray(a, dtype=float), np.asarray(e, dtype=float), np.asarray(inc, dtype=float)
        h  = self.torque.disk.aspect_ratio(a)
        tw = self.tau_wave(t_code, a)
        eh, ih = e / h, inc / h
        bracket = np.maximum(1.0 - 0.30 * ih**2 + 0.24 * ih**3 + 0.14 * eh**2 * ih, 1e-10)
        return tw / (0.544 * bracket)

    # ---- de/dt, di/dt ----

    def de_dt(self, t_code, a, e, inc):
        e     = np.asarray(e, dtype=float)
        tau_e = self.eccentricity_damping_timescale(t_code, a, e, inc)
        out   = np.zeros_like(e)
        mask  = (tau_e > 0.0) & np.isfinite(tau_e)
        out[mask] = -e[mask] / tau_e[mask]
        return float(out) if (np.isscalar(e) or e.shape == ()) else out

    def di_dt(self, t_code, a, e, inc):
        inc   = np.asarray(inc, dtype=float)
        tau_i = self.inclination_damping_timescale(t_code, a, e, inc)
        out   = np.zeros_like(inc)
        mask  = (tau_i > 0.0) & np.isfinite(tau_i)
        out[mask] = -inc[mask] / tau_i[mask]
        return float(out) if (np.isscalar(inc) or inc.shape == ()) else out

    # ---- CN08 eq.9：带偏心率修正的 da/dt ----

    def da_dt_eccentric(self, t_code, a, e, inc):
        a, e, inc = np.asarray(a, dtype=float), np.asarray(e, dtype=float), np.asarray(inc, dtype=float)
        eh      = e / self.torque.disk.aspect_ratio(a)
        P_e     = (1.0 + 0.5 * eh**5) / (1.0 + 0.25 * eh**5)
        da_circ = np.asarray(self.da_dt(t_code, a), dtype=float)
        de      = np.asarray(self.de_dt(t_code, a, e, inc), dtype=float)
        denom   = np.where(1.0 - e**2 > 1e-6, 1.0 - e**2, 1e-6)
        result  = da_circ / P_e + 2.0 * a * e * de / denom
        return float(result) if (np.isscalar(a) or a.shape == ()) else result

    # ---- 内部公共方法：生成时间轴 ----

    def _make_tcode(self):
        t_max_code = self.units.myr_to_code(self.cfg.t_max_myr)
        return np.linspace(0.0, t_max_code, self.cfg.n_steps)

    # ---- 圆轨道积分 ----

    def integrate(self):
        """圆轨道 RK4 积分。"""
        t_code = self._make_tcode()
        dt = t_code[1] - t_code[0]

        a_arr  = np.zeros(self.cfg.n_steps)
        gL_arr = np.zeros(self.cfg.n_steps)
        gB_arr = np.zeros(self.cfg.n_steps)
        gE_arr = np.zeros(self.cfg.n_steps)
        gT_arr = np.zeros(self.cfg.n_steps)
        a_arr[0]   = self.cfg.a0
        stop_index = self.cfg.n_steps - 1

        for i in range(self.cfg.n_steps - 1):
            t_now, a_now = t_code[i], a_arr[i]
            gL_arr[i] = self.torque.gamma_l(a_now, t_now)
            gB_arr[i] = self.torque.gamma_c_baro(a_now, t_now)
            gE_arr[i] = self.torque.gamma_c_ent(a_now, t_now)
            gT_arr[i] = self.torque.gamma_tot(a_now, t_now)

            if a_now <= self.cfg.a_min_stop or a_now >= self.cfg.a_max_stop:
                stop_index = i
                a_arr[i+1:] = a_now
                gL_arr[i+1:] = gL_arr[i]; gB_arr[i+1:] = gB_arr[i]
                gE_arr[i+1:] = gE_arr[i]; gT_arr[i+1:] = gT_arr[i]
                break

            k1 = self.da_dt(t_now,          a_now)
            k2 = self.da_dt(t_now + 0.5*dt, a_now + 0.5*dt*k1)
            k3 = self.da_dt(t_now + 0.5*dt, a_now + 0.5*dt*k2)
            k4 = self.da_dt(t_now + dt,      a_now + dt*k3)
            a_next = a_now + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

            if a_next <= 0.0:
                a_next = self.cfg.a_min_stop
                stop_index = i + 1
                a_arr[i+1] = a_next; a_arr[i+2:] = a_next
                break
            a_arr[i+1] = a_next

        sl = slice(0, stop_index + 1)
        return {
            "t_code":       t_code[sl],
            "t_myr":        self.units.code_to_myr(t_code[sl]),
            "a":            a_arr[sl],
            "gamma_l":      gL_arr[sl],
            "gamma_c_baro": gB_arr[sl],
            "gamma_c_ent":  gE_arr[sl],
            "gamma_tot":    gT_arr[sl],
        }

    # ---- 联立积分 a, e ----

    def integrate_eccentric(self, e0: float = 0.0, inc0: float = 0.0):
        """联立积分 a(t), e(t)，CN08 完整公式。"""
        t_code = self._make_tcode()
        dt = t_code[1] - t_code[0]

        a_arr  = np.zeros(self.cfg.n_steps)
        e_arr  = np.zeros(self.cfg.n_steps)
        gT_arr = np.zeros(self.cfg.n_steps)
        a_arr[0] = self.cfg.a0; e_arr[0] = e0
        stop_index = self.cfg.n_steps - 1

        def derivs(t, a, e):
            return self.da_dt_eccentric(t, a, e, inc0), self.de_dt(t, a, e, inc0)

        for i in range(self.cfg.n_steps - 1):
            t_now, a_now, e_now = t_code[i], a_arr[i], e_arr[i]
            gT_arr[i] = self.torque.gamma_tot(a_now, t_now)

            if a_now <= self.cfg.a_min_stop or a_now >= self.cfg.a_max_stop:
                stop_index = i
                a_arr[i+1:] = a_now; e_arr[i+1:] = e_now; gT_arr[i+1:] = gT_arr[i]
                break

            ka1,ke1 = derivs(t_now,          a_now,           e_now)
            ka2,ke2 = derivs(t_now + 0.5*dt, a_now+0.5*dt*ka1, e_now+0.5*dt*ke1)
            ka3,ke3 = derivs(t_now + 0.5*dt, a_now+0.5*dt*ka2, e_now+0.5*dt*ke2)
            ka4,ke4 = derivs(t_now + dt,     a_now+dt*ka3,     e_now+dt*ke3)

            a_next = a_now + (dt/6.0)*(ka1 + 2*ka2 + 2*ka3 + ka4)
            e_next = np.clip(e_now + (dt/6.0)*(ke1 + 2*ke2 + 2*ke3 + ke4), 0.0, 0.999)

            if a_next <= 0.0:
                stop_index = i + 1
                a_arr[i+1] = self.cfg.a_min_stop; e_arr[i+1] = e_next
                a_arr[i+2:] = self.cfg.a_min_stop; e_arr[i+2:] = e_next
                break
            a_arr[i+1] = a_next; e_arr[i+1] = e_next

        sl = slice(0, stop_index + 1)
        return {
            "t_code":    t_code[sl],
            "t_myr":     self.units.code_to_myr(t_code[sl]),
            "a":         a_arr[sl],
            "e":         e_arr[sl],
            "gamma_tot": gT_arr[sl],
        }

    # ---- 完整积分 a, e, inc ----

    def integrate_full(self, e0: float = 0.0, inc0: float = 0.0):
        """
        联立积分 a(t), e(t), inc(t)，CN08 完整公式。

        Parameters
        ----------
        e0   : 初始偏心率
        inc0 : 初始倾角（rad）
        """
        t_code = self._make_tcode()
        dt = t_code[1] - t_code[0]

        a_arr   = np.zeros(self.cfg.n_steps)
        e_arr   = np.zeros(self.cfg.n_steps)
        inc_arr = np.zeros(self.cfg.n_steps)
        gT_arr  = np.zeros(self.cfg.n_steps)
        a_arr[0] = self.cfg.a0; e_arr[0] = e0; inc_arr[0] = inc0
        stop_index = self.cfg.n_steps - 1

        def derivs(t, a, e, inc):
            return (self.da_dt_eccentric(t, a, e, inc),
                    self.de_dt(t, a, e, inc),
                    self.di_dt(t, a, e, inc))

        for i in range(self.cfg.n_steps - 1):
            t_now   = t_code[i]
            a_now, e_now, inc_now = a_arr[i], e_arr[i], inc_arr[i]
            gT_arr[i] = self.torque.gamma_tot(a_now, t_now)

            if a_now <= self.cfg.a_min_stop or a_now >= self.cfg.a_max_stop:
                stop_index = i
                a_arr[i+1:] = a_now; e_arr[i+1:] = e_now
                inc_arr[i+1:] = inc_now; gT_arr[i+1:] = gT_arr[i]
                break

            ka1,ke1,ki1 = derivs(t_now,          a_now,            e_now,            inc_now)
            ka2,ke2,ki2 = derivs(t_now+0.5*dt,   a_now+0.5*dt*ka1, e_now+0.5*dt*ke1, inc_now+0.5*dt*ki1)
            ka3,ke3,ki3 = derivs(t_now+0.5*dt,   a_now+0.5*dt*ka2, e_now+0.5*dt*ke2, inc_now+0.5*dt*ki2)
            ka4,ke4,ki4 = derivs(t_now+dt,        a_now+dt*ka3,     e_now+dt*ke3,     inc_now+dt*ki3)

            a_next   = a_now   + (dt/6.0)*(ka1 + 2*ka2 + 2*ka3 + ka4)
            e_next   = np.clip(e_now   + (dt/6.0)*(ke1 + 2*ke2 + 2*ke3 + ke4), 0.0, 0.999)
            inc_next = np.clip(inc_now + (dt/6.0)*(ki1 + 2*ki2 + 2*ki3 + ki4), 0.0, np.pi)

            if a_next <= 0.0:
                stop_index = i + 1
                a_arr[i+1] = self.cfg.a_min_stop
                e_arr[i+1] = e_next; inc_arr[i+1] = inc_next
                a_arr[i+2:] = self.cfg.a_min_stop
                e_arr[i+2:] = e_next; inc_arr[i+2:] = inc_next
                break
            a_arr[i+1] = a_next; e_arr[i+1] = e_next; inc_arr[i+1] = inc_next

        sl = slice(0, stop_index + 1)
        return {
            "t_code":    t_code[sl],
            "t_myr":     self.units.code_to_myr(t_code[sl]),
            "a":         a_arr[sl],
            "e":         e_arr[sl],
            "inc_rad":   inc_arr[sl],
            "inc_deg":   np.degrees(inc_arr[sl]),
            "gamma_tot": gT_arr[sl],
        }

    def fixed_radius_decay(self, r_fixed, t_code_array):
        t_code_array = np.asarray(t_code_array)
        return {
            "r_fixed":   r_fixed,
            "t_code":    t_code_array,
            "t_myr":     self.units.code_to_myr(t_code_array),
            "sigma":     self.torque.disk.sigma(r_fixed, t_code_array),
            "gamma_tot": self.torque.gamma_tot(r_fixed, t_code_array),
        }


# =========================================================
# 6. Factory
# =========================================================

def Gas_build_model(cfg: Gas_ModelConfig | None = None) -> dict:
    """
    一键构建完整模型。

    Returns
    -------
    dict with keys: config, units, disk, torque, integrator
    """
    cfg        = cfg or Gas_ModelConfig()
    units      = Gas_UnitSystem()
    disk       = Gas_DiskModel(cfg, units=units)
    torque     = Gas_TorqueModel(cfg, disk=disk)
    integrator = Gas_MigrationIntegrator(cfg, torque_model=torque, units=units)
    return {
        "config":     cfg,
        "units":      units,
        "disk":       disk,
        "torque":     torque,
        "integrator": integrator,
    }