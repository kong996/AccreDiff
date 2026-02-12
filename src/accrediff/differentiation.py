# accrediff/differentiation.py
import pandas as pd # type: ignore
import numpy as np # type: ignore
from dataclasses import dataclass, replace
from typing import Literal, Optional, Tuple, Dict, Callable
import math
from .utils import power_law
from scipy.optimize import minimize_scalar 

# =========================
# 数据结构
# =========================
@dataclass
class KD_Params:
    Fe_t: float
    Ni_t: float
    Si_t: float        # alpha
    O_L: float       # constant oxygen Loss
    KD_Ni: float
    KD_Si: float
    u: float = 0.0
    m: float = 0.0
    n: float = 0.0

@dataclass
class KD_Result:
    x_: float
    a_: float
    y_: float
    b_: float
    z_: float
    c_: float
    d_: float
    KD_O_: float
    residual: float

#**************************************************************************************************************************************
def Early_pressure(mass, a=112.06, b=0.37):
    """
    计算Early phase impacts equilibrium pressure (Gu et al. 2023)
    :param mass: 质量
    :param a: 参数a
    :param b: 参数b
    :return: equilibrium pressure
    """
    P_emb = a*mass + b 
    return P_emb
#**************************************************************************************************************************************
def CMB_pressure(mass, a=136.50, b=0.91):
        """
        Calculate the pressure at the CMB for a given mass.
        
        Parameters:
        mass (float): Mass in Earth masses (ME).
        
        Returns:
        float: Pressure in GPa.
        """
        P = power_law(mass, a, b) 
        return P
#**************************************************************************************************************************************
def equil_pressure(mass, ratio):
    """
    计算平衡压力
    :param mass: 质量
    :param ratio: 比例
    :return: 平衡压力
    """
    P_equil = CMB_pressure(mass) * ratio
    return P_equil
#**************************************************************************************************************************************
def P_to_T(P):
    """Convert pressure in GPa to temperature in Kelvin.Rubie et al. 2015"""
    if P < 24:
        T = 1874 + 55.43*P - 1.74*(P**2) + 0.0193*(P**3)
    else:
        T = 1249 + 58.28*P - 0.395*(P**2) + 0.0011*(P**3)
    return T
#**************************************************************************************************************************************
class EarlyComUpdater:
    """
    用于批量更新early_com DataFrame中各粒子的化学参数。
    依赖自定义库 wisdom (import wisdom as wd)。
    """
    def __init__(self, df, dict_com, keys_list, accrediff, m_ratio=True):
        self.df = df
        self.dict_com  = dict_com
        self.keys_list = keys_list
        self.accrediff = accrediff
        self.m_ratio = m_ratio
        # 自动获取所有成分类型（如['IW15','IW35']），无需硬编码
        #self.com_types = [k for k in dict_IW.keys() if isinstance(dict_IW[k], dict) and k in df.columns]
        self.com_types = list(dict_com.keys())
        #print(self.com_types)
        self.key_map = {
            'x_': 'FeO',
            'y_': 'NiO',
            'z_': 'SiO2',
            'a_': 'Fe',
            'b_': 'Ni',
            'c_': 'Si',
            'd_': 'O'
        }
        self.update_cols = ['FeO', 'SiO2', 'NiO', 'Fe', 'Ni', 'Si', 'O', 'MgO', 'Al2O3', 'CaO']

    def update_row(self, pid):
        row = self.df[self.df['i'] == pid].iloc[0]
        P, T = row['P'], row['T']
        # 动态读取所有成分类型的摩尔分数
        if self.m_ratio == True:
        # 只选取 ratio 列中，且在 dict_com（dict_IW）中存在的成分类型
            ratio_cols = [col for col in self.df.columns if col.endswith('_ratio') and col.replace('_ratio', '') in self.com_types]
            m_dict = {col.replace('_ratio', ''): row[col] for col in ratio_cols}
        else:
            m_dict = {ctype: row[ctype] for ctype in self.com_types}     
        # 计算 bulk        
        bulk = {k: sum(self.dict_com[ctype][k] * m_dict[ctype] for ctype in self.com_types) for k in self.keys_list}
        #print(f'Particle ID {pid} bulk composition:', bulk)
        Mg_total = bulk['MgO']
        Al_total = bulk['Al2O3']*2  # 转化为 AlO1.5
        Ca_total = bulk['CaO']
        Fe_total = bulk['FeO'] + bulk['Fe']
        Ni_total = bulk['NiO'] + bulk['Ni']
        Si_total = bulk['SiO2'] + bulk['Si']
        O_L = bulk['Fe'] + bulk['Ni'] + 2 * bulk['Si'] - bulk['O']
        # 计算 KD
        kd_calc = self.accrediff.KDCalculator()
        KD_Ni = kd_calc.get_KD('Ni', P, T)
        KD_Si = kd_calc.get_KD('Si', P, T)
        KD_O = kd_calc.get_KD('O', P, T)
        # 构建参数
        params_dict = {
            "Fe_t": Fe_total,
            "Ni_t": Ni_total,
            "Si_t": Si_total,
            "O_L": O_L,
            "KD_Ni": KD_Ni,
            "KD_Si": KD_Si,
            "u": Mg_total,
            "m": Al_total,
            "n": Ca_total,
        }
        p0 = self.accrediff.KD_Params(**params_dict)
        solver = self.accrediff.ForwardKDOSolver(
            p0,
            nonneg='clip',
            enforce_z_box=True,
            enforce_d_nonneg=True,
        )
        res = solver.solve_x_for_KD_O(KD_O, tol=1e-12, max_iter=300, grid_N=400)
        res_dict = vars(res)
        res_dict_new = {self.key_map[k]: v for k, v in res_dict.items() if k in self.key_map}
        res_dict_new['MgO'] = Mg_total
        res_dict_new['Al2O3'] = Al_total/2 # 转化回 Al2O3
        res_dict_new['CaO'] = Ca_total
        row_idx = self.df[self.df['i'] == pid].index[0]
        for col in self.update_cols:
            if col in res_dict_new:
                self.df.at[row_idx, col] = res_dict_new[col]

    def batch_update(self):
        for pid in self.df['i']:
            self.update_row(pid)

#**************************************************************************************************************************************
class ForwardKDOSolver:
    """
    Forward + 反向（匹配 KD_O）求解器。
    """

    def __init__(
        self,
        params: KD_Params,
        *,
        nonneg: Literal["clip", "softplus", "none"] = "clip",
        enforce_z_box: bool = False,
        enforce_d_nonneg: bool = False

    ) -> None:
        self.p = params
        self.nonneg = nonneg
        self.enforce_z_box = enforce_z_box
        self.enforce_d_nonneg = enforce_d_nonneg

    # ---------------- 公共接口 ----------------
    def forward_solve(self, x_init: float) -> Dict[str, float]:
        """
        给定 x_init，按依赖顺序求解 (x_, a_, y_, b_, z_, c_, d_)。
        支持 enforce_z_box 与 enforce_d_nonneg 开关。
        """
        p = self.p
        x_ = min(max(x_init, 0.0), p.Fe_t)

        a_ = self._project_nonneg(self._f_a(x_), self.nonneg)
        y_ = self._project_nonneg(self._f_y(x_, a_), self.nonneg)

        b_raw = self._f_b(y_)
        b_ = self._project_nonneg(min(max(b_raw, 0.0), p.Ni_t), self.nonneg)

        z_raw = self._f_z(x_, y_, a_, b_, enforce_box=self.enforce_z_box)
        if self.enforce_z_box:
            z_raw = min(max(z_raw, 0.0), p.Si_t)
        z_ = self._project_nonneg(z_raw, self.nonneg)

        c_ = self._project_nonneg(self._f_c(z_), self.nonneg)

        d_raw = self._f_d(x_, y_, z_)
        if self.enforce_d_nonneg:
            d_raw = max(d_raw, 0.0)
        d_ = self._project_nonneg(d_raw, self.nonneg)

        return dict(x_=x_, a_=a_, y_=y_, b_=b_, z_=z_, c_=c_, d_=d_)

    def compute_KD_O_from_x(self, x_: float) -> Tuple[float, Dict[str, float]]:
        p = self.p
        sol = self.forward_solve(x_)
        denom = sol["x_"] + sol["y_"] + sol["z_"] + p.u + p.m + p.n
        X_Sil_FeO = self._safe_div(sol["x_"], denom)
        X_MW_FeO = 1.148 * X_Sil_FeO + 1.319 * (X_Sil_FeO ** 2)
        metal_sum = sol["a_"] + sol["b_"] + sol["c_"] + sol["d_"]
        KD_O_ = self._safe_div(sol["a_"] * sol["d_"], X_MW_FeO * (metal_sum ** 2))
        return KD_O_, sol

    def solve_x_for_KD_O(
        self,
        KD_O_target: float,
        *,
        tol: float = 1e-12,
        max_iter: int = 200,
        grid_N: int = 200
    ) -> KD_Result:
        p = self.p
        a, b = 0.0, p.Fe_t

        xs = [a + (b - a) * i / grid_N for i in range(grid_N + 1)]
        vals = [(self.compute_KD_O_from_x(xi)[0] - KD_O_target) for xi in xs]

        bracket = None
        for i in range(grid_N):
            if vals[i] == 0.0:
                KD_best, sol_best = self.compute_KD_O_from_x(xs[i])
                return KD_Result(KD_O_=KD_best, residual=0.0, **sol_best)
            if vals[i] * vals[i+1] < 0:
                bracket = (xs[i], xs[i+1], vals[i], vals[i+1])
                break

        if bracket is None:
            idx = min(range(len(xs)), key=lambda k: abs(vals[k]))
            x_best = xs[idx]
            KD_best, sol_best = self.compute_KD_O_from_x(x_best)
            return KD_Result(KD_O_=KD_best, residual=KD_best - KD_O_target, **sol_best)

        a, b, fa, fb = bracket
        mid = 0.5 * (a + b)
        KD_mid, sol_mid = self.compute_KD_O_from_x(mid)
        fm = KD_mid - KD_O_target
        for _ in range(max_iter):
            mid = 0.5 * (a + b)
            KD_mid, sol_mid = self.compute_KD_O_from_x(mid)
            fm = KD_mid - KD_O_target
            if abs(fm) < tol or abs(b - a) < tol:
                return KD_Result(KD_O_=KD_mid, residual=fm, **sol_mid)
            if fa * fm <= 0:
                b, fb = mid, fm
            else:
                a, fa = mid, fm

        return KD_Result(KD_O_=KD_mid, residual=fm, **sol_mid)

    # ---------------- 私有工具 ----------------
    def _project_nonneg(self, v: float, mode: Literal["clip", "softplus", "none"]) -> float:
        if mode == "clip":
            return v if v >= 0.0 else 0.0
        elif mode == "softplus":
            return math.log1p(math.exp(v))
        else:
            return v

    @staticmethod
    def _safe_div(num: float, den: float, eps: float = 1e-15) -> float:
        return num / (den if abs(den) > eps else (eps if den >= 0 else -eps))

    @staticmethod
    def _quadratic_roots(A: float, B: float, C: float) -> Optional[tuple]:
        if abs(A) < 1e-18:
            if abs(B) < 1e-18:
                return None
            z = ForwardKDOSolver._safe_div(C, B)
            return (z, z)
        D = B * B - 4.0 * A * C
        if D < 0:
            return None
        sqrtD = math.sqrt(D)
        r1 = (B + sqrtD) / (2.0 * A)
        r2 = (B - sqrtD) / (2.0 * A)
        return (r1, r2)

    def _choose_physical_root(self, candidates, alpha: float, *, enforce_box: bool = False) -> float:
        cand = [z for z in candidates if z is not None]
        if not enforce_box:
            in_box = [z for z in cand if 0.0 <= z <= alpha]
            if in_box:
                return max(in_box)
            nonneg = [z for z in cand if z >= 0.0]
            if nonneg:
                return max(nonneg)
            return 0.0
        in_box = [z for z in cand if 0.0 <= z <= alpha]
        if in_box:
            return max(in_box)
        if not cand:
            return 0.0
        def dist_to_box(z):
            if z < 0.0:
                return -z
            if z > alpha:
                return z - alpha
            return 0.0
        z_near = min(cand, key=dist_to_box)
        if z_near < 0.0:
            return 0.0
        if z_near > alpha:
            return alpha
        return z_near

    def _total_O(self) -> float:
        p = self.p
        return p.Fe_t + p.Ni_t + 2.0 * p.Si_t - p.O_L

    # ------- 方程 -------
    def _f_a(self, x_: float) -> float:
        return self.p.Fe_t - x_

    def _f_y(self, x_: float, a_: float) -> float:
        den = a_ * self.p.KD_Ni + x_
        return self._safe_div(x_ * self.p.Ni_t, den)

    def _f_b(self, y_: float) -> float:
        return self.p.Ni_t - y_

    def _f_z(self, x_: float, y_: float, a_: float, b_: float, *, enforce_box: bool = False) -> float:
        p = self.p
        alpha = p.Si_t
        Ot = self._total_O()
        gamma = a_ + b_ + Ot + p.Si_t - x_ - y_
        sigma = x_ + y_ + p.u + p.m + p.n
        A = 3.0 * x_ * x_ - (a_ * a_) * p.KD_Si
        B = gamma * x_ * x_ + 3.0 * alpha * x_ * x_ + (a_ * a_) * sigma * p.KD_Si
        C = alpha * gamma * x_ * x_
        roots = self._quadratic_roots(A, B, C)
        if roots is None:
            return 0.0
        return self._choose_physical_root(roots, alpha=alpha, enforce_box=enforce_box)

    def _f_c(self, z_: float) -> float:
        return self.p.Si_t - z_

    def _f_d(self, x_: float, y_: float, z_: float) -> float:
        return self._total_O() - x_ - y_ - 2.0 * z_
#**************************************************************************************************************************************    
class OLSolver:
    """
    O_L 反向求解器：寻找使 IW 与目标值一致的最佳 O_L。
    使用 ForwardKDOSolver 作为内核。
    """

    def __init__(
        self,
        base_params: KD_Params,
        KD_O_target: float,
        IW_target: float,
        solver_factory: Callable[[KD_Params], ForwardKDOSolver],
        nonneg: Literal["clip", "softplus", "none"] = "clip",
        enforce_z_box: bool = False,
        enforce_d_nonneg: bool = False,
    ):
        self.base_params = base_params
        self.KD_O_target = KD_O_target
        self.IW_target = IW_target
        self.solver_factory = solver_factory
        self.nonneg = nonneg
        self.enforce_z_box = enforce_z_box
        self.enforce_d_nonneg = enforce_d_nonneg

    def _compute_IW(self, sol: Dict[str, float]) -> float:
        """
        根据正向求解得到的元素分布，计算 IW = 2 * log10(X_FeO / X_Fe)
        """
        sil_sum = sol["x_"] + sol["y_"] + sol["z_"] + self.base_params.u + self.base_params.m + self.base_params.n
        metal_sum = sol["a_"] + sol["b_"] + sol["c_"] + sol["d_"]

        eps = 1e-12  # 避免 log(0)
        X_FeO = max(sol["x_"] / sil_sum, eps) if sil_sum > 0 else eps
        X_Fe = max(sol["a_"] / metal_sum, eps) if metal_sum > 0 else eps

        return 2.0 * math.log10(X_FeO / X_Fe)

    def _residual(self, O_L_guess: float) -> float:
        """
        计算当前猜测的 O_L 对应的 IW 残差
        """
        try:
            new_params = replace(self.base_params, O_L=O_L_guess)
            solver = self.solver_factory(new_params)
            kd_result = solver.solve_x_for_KD_O(self.KD_O_target)
            IW_model = self._compute_IW(kd_result.__dict__)
            return abs(IW_model - self.IW_target)
        except Exception:
            return 1e6  # 若计算失败，返回大残差

    def solve(
        self,
        O_bounds: Tuple[float, float] = (-10.0, 10.0),
        tol: float = 1e-6,
        max_iter: int = 200
    ) -> Tuple[float, float]:
        """
        使用 scalar 最优化方法搜索最佳 O_L 值
        返回：(最优 O_L, IW 残差)
        """
        result = minimize_scalar(
            self._residual,
            bounds=O_bounds,
            method="bounded",
            options={"xatol": tol, "maxiter": max_iter}
        )
        return result.x, result.fun
#**************************************************************************************************************************************