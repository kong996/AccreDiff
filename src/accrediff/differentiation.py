# accrediff/differentiation.py
import pandas as pd # type: ignore
import numpy as np # type: ignore
from dataclasses import dataclass, replace
from typing import Literal, Optional, Tuple, Dict, Callable, Any
import math
import copy
#from .utils import power_law
from scipy.optimize import minimize_scalar  # type: ignore

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
def power_law(x, a, b):
    """Power law function."""
    return a * (x ** b)
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
    def __init__(self, df, dict_com, keys_list, accrediff, m_ratio=True, 
                 KD_oxygen: Literal["Rubie", "Fischer"] = "Rubie"  # 修改：加 Literal 类型注解
                 ):
        self.df = df
        self.dict_com  = dict_com
        self.keys_list = keys_list
        self.accrediff = accrediff
        self.m_ratio = m_ratio
        self.KD_oxygen = KD_oxygen
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
            KD_oxygen=self.KD_oxygen,
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
    金属–硅酸盐分异正向求解器 + KD_O 反向匹配。

    物理背景
    --------
    在行星吸积过程中，金属相（核）与硅酸盐相（幔）之间的元素分配受
    分配系数（KD）控制。本求解器给定总组成和 KD_Ni、KD_Si 后，以
    硅酸盐中 FeO 含量 x_ 为自由变量，依次求解所有组分在两相间的分布，
    再通过网格搜索 + 二分法找到使计算的 KD_O 匹配目标值的 x_。

    变量符号约定（摩尔量）
    ----------------------
    硅酸盐相:  x_ = FeO,  y_ = NiO,  z_ = SiO2
    金属相:    a_ = Fe,   b_ = Ni,   c_ = Si,   d_ = O
    惰性组分:  u = MgO,   m = AlO1.5, n = CaO  （全部在硅酸盐相）

    守恒关系
    --------
    Fe_t = x_ + a_      （铁守恒）
    Ni_t = y_ + b_      （镍守恒）
    Si_t = z_ + c_      （硅守恒，alpha = Si_t）
    O_total = x_ + y_ + 2z_ + d_  （参与分配的氧守恒，不含 u/m/n 携带的氧）
    其中 O_total = Fe_t + Ni_t + 2*Si_t - O_L

    参数
    ----
    params : KD_Params
        包含总组成 (Fe_t, Ni_t, Si_t)、氧损失 (O_L)、
        分配系数 (KD_Ni, KD_Si) 和惰性组分 (u, m, n)。
    nonneg : "clip" | "softplus" | "none"
        非负约束方式。"clip" 将负值截断为 0。
    enforce_z_box : bool
        是否强制 z_ 在 [0, Si_t] 范围内。
    enforce_d_nonneg : bool
        是否强制 d_（金属中氧）非负。
    """

    def __init__(
        self,
        params: KD_Params,
        *,
        nonneg: Literal["clip", "softplus", "none"] = "clip",
        enforce_z_box: bool = False,
        enforce_d_nonneg: bool = False,
        KD_oxygen: Literal["Rubie", "Fischer"] = "Rubie"

    ) -> None:
        self.p = params
        self.nonneg: Literal["clip", "softplus", "none"] = nonneg
        self.enforce_z_box = enforce_z_box
        self.enforce_d_nonneg = enforce_d_nonneg
        self.KD_oxygen = KD_oxygen

    # ======================== 公共接口 ========================

    def forward_solve(self, x_init: float) -> Dict[str, float]:
        """
        给定 x_init（硅酸盐中 FeO 摩尔量），按依赖链求解所有组分。

        求解顺序（每一步仅依赖前面已求得的量）：
          1. x_ = clamp(x_init, 0, Fe_t)
          2. a_ = Fe_t - x_                          （铁守恒）
          3. y_ = f(x_, a_)  由 KD_Ni 平衡方程求解    （镍在硅酸盐中）
          4. b_ = Ni_t - y_                           （镍守恒）
          5. z_ = f(x_, y_, a_, b_)  由 KD_Si 二次方程求解（硅在硅酸盐中）
          6. c_ = Si_t - z_                           （硅守恒）
          7. d_ = O_total - x_ - y_ - 2*z_           （氧守恒）

        Returns
        -------
        dict : {"x_", "a_", "y_", "b_", "z_", "c_", "d_"}
        """
        p = self.p
        # Step 1: 限制 x_ 在物理范围 [0, Fe_t]
        x_ = min(max(x_init, 0.0), p.Fe_t)

        # Step 2: 金属中 Fe = 总 Fe - 硅酸盐中 FeO
        a_ = self._project_nonneg(self._f_a(x_), self.nonneg)

        # Step 3: 由 KD_Ni 平衡求硅酸盐中 NiO
        y_ = self._project_nonneg(self._f_y(x_, a_), self.nonneg)

        # Step 4: 金属中 Ni = 总 Ni - 硅酸盐中 NiO
        b_raw = self._f_b(y_)
        b_ = self._project_nonneg(min(max(b_raw, 0.0), p.Ni_t), self.nonneg)

        # Step 5: 由 KD_Si 平衡（二次方程）求硅酸盐中 SiO2
        z_raw = self._f_z(x_, y_, a_, b_, enforce_box=self.enforce_z_box)
        if self.enforce_z_box:
            z_raw = min(max(z_raw, 0.0), p.Si_t)
        z_ = self._project_nonneg(z_raw, self.nonneg)

        # Step 6: 金属中 Si = 总 Si - 硅酸盐中 SiO2
        c_ = self._project_nonneg(self._f_c(z_), self.nonneg)

        # Step 7: 金属中 O 由氧守恒确定
        d_raw = self._f_d(x_, y_, z_)
        if self.enforce_d_nonneg:
            d_raw = max(d_raw, 0.0)
        d_ = self._project_nonneg(d_raw, self.nonneg)

        return dict(x_=x_, a_=a_, y_=y_, b_=b_, z_=z_, c_=c_, d_=d_)

    def compute_KD_O_from_x(self, x_: float) -> Tuple[float, Dict[str, float]]:
        """
        给定 x_，正向求解后计算 KD_O。

        KD_O 定义 (Fischer et al. 2017):
            KD_O = X_Fe_met * X_O_met / X_FeO_sil_MW

        其中:
            X_Fe_met  = a_ / metal_sum
            X_O_met   = d_ / metal_sum
            X_FeO_sil = x_ / sil_sum       （理想摩尔分数）
            X_FeO_MW  = 1.148 * X_FeO_sil + 1.319 * X_FeO_sil²  （Margules 非理想修正活度）

        合并后:
            KD_O = (a_ * d_) / (X_FeO_MW * metal_sum²)

        Returns
        -------
        (KD_O_, sol) : (float, dict)
        """
        p = self.p
        sol = self.forward_solve(x_)

        # 硅酸盐相总摩尔数 = FeO + NiO + SiO2 + MgO + AlO1.5 + CaO
        denom = sol["x_"] + sol["y_"] + sol["z_"] + p.u + p.m + p.n
        X_Sil_FeO = self._safe_div(sol["x_"], denom)

        # Rubie et al. 2011: FeO 的活度系数修正
        if self.KD_oxygen == "Rubie":
            X_MW_FeO = 1.148 * X_Sil_FeO + 1.319 * (X_Sil_FeO ** 2)
        elif self.KD_oxygen == "Fischer":
            X_MW_FeO = X_Sil_FeO  # 先用理想近似 (Fischer et al. 2017 中 KD_O 定义即为理想近似)，后续可替换为非理想修正
        else:
            raise ValueError(f"Unsupported KD_oxygen model: {self.KD_oxygen}")

        # 金属相总摩尔数 = Fe + Ni + Si + O
        metal_sum = sol["a_"] + sol["b_"] + sol["c_"] + sol["d_"]

        # KD_O = X_Fe * X_O / activity_FeO = (a_*d_) / (X_MW_FeO * metal_sum²)
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
        """
        搜索使 compute_KD_O_from_x(x_) == KD_O_target 的 x_ 值。

        算法
        ----
        1. 在 [0, Fe_t] 上均匀撒 grid_N+1 个点，计算 KD_O(x) - target
        2. 寻找第一个变号区间（bracket）
        3. 若找到 bracket，用二分法精确求解
        4. 若无 bracket，返回网格上残差最小的点（近似解）

        Parameters
        ----------
        KD_O_target : float
            目标 KD_O 值。
        tol : float
            收敛容差（残差或区间宽度）。
        max_iter : int
            二分法最大迭代次数。
        grid_N : int
            初始网格点数。

        Returns
        -------
        KD_Result : 包含 x_, a_, y_, b_, z_, c_, d_, KD_O_, residual
        """
        p = self.p
        a, b = 0.0, p.Fe_t

        # --- 第 1 步: 网格搜索 ---
        xs = [a + (b - a) * i / grid_N for i in range(grid_N + 1)]
        vals = [(self.compute_KD_O_from_x(xi)[0] - KD_O_target) for xi in xs]

        # --- 第 2 步: 寻找变号区间 ---
        bracket = None
        for i in range(grid_N):
            if vals[i] == 0.0:
                # 恰好命中目标
                KD_best, sol_best = self.compute_KD_O_from_x(xs[i])
                return KD_Result(KD_O_=KD_best, residual=0.0, **sol_best)
            if vals[i] * vals[i+1] < 0:
                # 找到第一个变号区间
                bracket = (xs[i], xs[i+1], vals[i], vals[i+1])
                break

        # --- 无变号区间: 返回网格上最接近的点 ---
        if bracket is None:
            idx = min(range(len(xs)), key=lambda k: abs(vals[k]))
            x_best = xs[idx]
            KD_best, sol_best = self.compute_KD_O_from_x(x_best)
            return KD_Result(KD_O_=KD_best, residual=KD_best - KD_O_target, **sol_best)

        # --- 第 3 步: 二分法精确求解 ---
        a, b, fa, fb = bracket
        KD_mid, sol_mid, fm = 0.0, {}, 0.0  # 初始化，防止 max_iter=0 时未绑定
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

    # ======================== 私有工具方法 ========================

    def _project_nonneg(self, v: float, mode: Literal["clip", "softplus", "none"]) -> float:
        """非负投影: clip→截断, softplus→光滑近似, none→不处理。"""
        if mode == "clip":
            return v if v >= 0.0 else 0.0
        elif mode == "softplus":
            return math.log1p(math.exp(v))
        else:
            return v

    @staticmethod
    def _safe_div(num: float, den: float, eps: float = 1e-15) -> float:
        """安全除法，避免除以零。当 |den| < eps 时用 ±eps 替代。"""
        return num / (den if abs(den) > eps else (eps if den >= 0 else -eps))

    @staticmethod
    def _quadratic_roots(A: float, B: float, C: float) -> Optional[tuple]:
        """
        求解一元二次方程 Az² + Bz + C = 0。

        退化情况:
          - A ≈ 0 且 B ≈ 0 → 无解 (None)
          - A ≈ 0 → 线性方程 z = -C/B
        正常情况:
          - D < 0 → 无实根 (None)
          - D >= 0 → 返回 (r1, r2)，r1 = (-B+√D)/(2A), r2 = (-B-√D)/(2A)
        """
        if abs(A) < 1e-18:
            if abs(B) < 1e-18:
                return None
            z = -ForwardKDOSolver._safe_div(C, B)
            return (z, z)
        D = B * B - 4.0 * A * C
        if D < 0:
            return None
        sqrtD = math.sqrt(D)
        r1 = (-B + sqrtD) / (2.0 * A)
        r2 = (-B - sqrtD) / (2.0 * A)
        return (r1, r2)

    def _choose_physical_root(self, candidates, alpha: float, *, enforce_box: bool = False) -> float:
        """
        从二次方程的候选根中选择物理合理的 z_ 值。

        选择策略:
          1. 优先选 [0, alpha] 区间内的最大根（更多 Si 在硅酸盐相）
          2. 若无 box 内的根且 enforce_box=False: 选非负最大根，否则返回 0
          3. 若 enforce_box=True: 选离 [0, alpha] 最近的根并 clip 到边界
        """
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
            """计算 z 到 [0, alpha] 区间的距离。"""
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

    # ======================== 守恒与平衡方程 ========================

    def _total_O(self) -> float:
        """
        参与分配的总氧摩尔数（不含 MgO/AlO1.5/CaO 中的氧）。

        O_total = Fe_t + Ni_t + 2*Si_t - O_L

        物理含义: 原始组成中与 Fe/Ni/Si 配对的氧，减去损失量 O_L。
        这些氧分布为: 硅酸盐中 (x_ + y_ + 2z_) + 金属中 (d_)。
        """
        p = self.p
        return p.Fe_t + p.Ni_t + 2.0 * p.Si_t - p.O_L

    def _f_a(self, x_: float) -> float:
        """Fe 守恒: a_ (金属中 Fe) = Fe_t - x_ (硅酸盐中 FeO)"""
        return self.p.Fe_t - x_

    def _f_y(self, x_: float, a_: float) -> float:
        """
        由 KD_Ni 平衡求硅酸盐中 NiO (y_)。

        KD_Ni 定义 (Fischer et al. 2017):
            KD_Ni = (X_Ni_met / X_NiO_sil) * (X_FeO_sil / X_Fe_met)
                  = (b_ / y_) * (x_ / a_)    （简化，假设相总量约消）

        化简为:  KD_Ni = (Ni_t - y_) * x_ / (y_ * a_)
        解出:    y_ = x_ * Ni_t / (a_ * KD_Ni + x_)
        """
        den = a_ * self.p.KD_Ni + x_
        return self._safe_div(x_ * self.p.Ni_t, den)

    def _f_b(self, y_: float) -> float:
        """Ni 守恒: b_ (金属中 Ni) = Ni_t - y_ (硅酸盐中 NiO)"""
        return self.p.Ni_t - y_

    def _f_z(self, x_: float, y_: float, a_: float, b_: float, *, enforce_box: bool = False) -> float:
        """
        由 KD_Si 平衡求硅酸盐中 SiO2 (z_)，需解二次方程。

        KD_Si 定义 (from KDCalculator, Fischer et al. 2017):
            KD_Si = X_Si_met * X_FeO_sil² / (X_SiO2_sil * X_Fe_met²)

        其中摩尔分数为:
            X_Si_met  = c_ / metal_sum = (alpha - z_) / (gamma - 3z_)
            X_Fe_met  = a_ / metal_sum = a_ / (gamma - 3z_)
            X_FeO_sil = x_ / sil_sum   = x_ / (sigma + z_)
            X_SiO2_sil= z_ / sil_sum   = z_ / (sigma + z_)

        辅助量:
            alpha = Si_t                                    （总 Si）
            gamma = a_ + b_ + O_total + Si_t - x_ - y_     （= metal_sum + 3z_）
            sigma = x_ + y_ + u + m + n                    （= sil_sum - z_）

        代入 KD_Si 定义并交叉相乘:
            (alpha - z_) · x_² · (gamma - 3z_) = KD_Si · z_ · a_² · (sigma + z_)

        展开整理为标准二次方程 Az² + Bz + C = 0:
            A = 3·x² - KD·a²
            B = -((gamma + 3·alpha)·x² + KD·a²·sigma)
            C = alpha·gamma·x²

        物理验证:
          - 当 KD_Si → 0（Si 不亲铁，低压条件）: 两个根为 alpha 和 gamma/3，
            选 max → z_ ≈ alpha，即几乎所有 Si 留在硅酸盐相 ✓
          - 当 KD_Si → ∞（Si 极度亲铁）: z_ → 0，Si 全部进入金属相 ✓
        """
        p = self.p
        alpha = p.Si_t
        Ot = self._total_O()
        gamma = a_ + b_ + Ot + p.Si_t - x_ - y_
        sigma = x_ + y_ + p.u + p.m + p.n

        a2 = a_ * a_
        x2 = x_ * x_
        KD = p.KD_Si
        A = 3.0 * x2 - KD * a2
        B = -(( gamma + 3.0 * alpha) * x2 + KD * a2 * sigma)
        C = alpha * gamma * x2

        roots = self._quadratic_roots(A, B, C)
        if roots is None:
            return 0.0
        return self._choose_physical_root(roots, alpha=alpha, enforce_box=enforce_box)

    def _f_c(self, z_: float) -> float:
        """Si 守恒: c_ (金属中 Si) = Si_t - z_ (硅酸盐中 SiO2)"""
        return self.p.Si_t - z_

    def _f_d(self, x_: float, y_: float, z_: float) -> float:
        """氧守恒: d_ (金属中 O) = O_total - x_ - y_ - 2·z_"""
        return self._total_O() - x_ - y_ - 2.0 * z_
#**************************************************************************************************************************************    

class OLSolver:
    """改进版本：精度参数简化 + 结果缓存 + 递归修复"""
    
    def __init__(
        self,
        base_params: KD_Params,
        KD_O_target: float,
        IW_target: float,
        solver_factory: Callable[[KD_Params], ForwardKDOSolver],
        precision: Literal["fast", "normal", "high"] = "normal",
        outer_tol: Optional[float] = None,
        max_outer_iter: int = 300,
    ):
        """
        Parameters
        ----------
        precision : str
            精度预设级别: "fast", "normal", "high"
        outer_tol : float, optional
            覆盖 precision 预设的外层精度
        """
        self.base_params = base_params
        self.KD_O_target = KD_O_target
        self.IW_target = IW_target
        self.solver_factory = solver_factory
        self.max_outer_iter = max_outer_iter
        
        # ── 精度预设 ──
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
                "inner_tol": 1e-14,
                "inner_max_iter": 500,
                "inner_grid_N": 400,
                "outer_tol": 1e-8,
            },
        }
        
        preset = precision_presets[precision]
        self.inner_tol = preset["inner_tol"]
        self.inner_max_iter = preset["inner_max_iter"]
        self.inner_grid_N = preset["inner_grid_N"]
        self.outer_tol = outer_tol if outer_tol is not None else preset["outer_tol"]
        
        # ── 结果缓存 ──
        self._cached_result: Optional[Tuple[Dict, float]] = None
        self._best_residual: float = float('inf')  # ← 追踪最小残差

    def _residual(self, O_L_guess: float) -> float:
        """
        计算残差，并缓存最优解。
        
        ⚠️ 不再在此调用自身！避免无限递归。
        """
        try:
            # ── 创建新参数（使用 dataclass replace 或深拷贝）
            # 假设 KD_Params 是 dataclass
            from copy import deepcopy
            new_params = deepcopy(self.base_params)
            new_params.O_L = O_L_guess
            
            # ── 内层求解
            solver = self.solver_factory(new_params)
            kd_result = solver.solve_x_for_KD_O(
                self.KD_O_target,
                tol=self.inner_tol,
                max_iter=self.inner_max_iter,
                grid_N=self.inner_grid_N,
            )
            res_dict = vars(kd_result) if hasattr(kd_result, '__dict__') else kd_result
            
            # ── 确保 res_dict 是字典类型
            if not isinstance(res_dict, dict):
                # 如果 KD_Result 有转换方法，使用它
                res_dict = vars(res_dict)
            # ── 计算 IW
            IW_model = self._compute_IW(res_dict)
            residual = abs(IW_model - self.IW_target)
            
            # ── 无条件缓存所有计算结果（不做递归比较！）
            if residual < self._best_residual:
                self._best_residual = residual
                self._cached_result = (res_dict, IW_model)
            
            return residual
            
        except RecursionError as e:
            print(f"  ERROR: RecursionError at O_L={O_L_guess}: {e}")
            return 1e6
        except Exception as e:
            print(f"  Warning: residual computation failed at O_L={O_L_guess}: {type(e).__name__}")
            return 1e6

    def _compute_IW(self, res_dict: Dict[str, float]) -> float:
        """
        计算 IW 指标（FeO-Fe 缓冲）。
        
        IW = 2 * log10(X_FeO_sil / X_Fe_met)
        """
        p = self.base_params
        x_ = res_dict.get('x_', 0.0)
        y_ = res_dict.get('y_', 0.0)
        z_ = res_dict.get('z_', 0.0)
        a_ = res_dict.get('a_', 0.0)
        b_ = res_dict.get('b_', 0.0)
        c_ = res_dict.get('c_', 0.0)
        d_ = res_dict.get('d_', 0.0)
        
        # 硅酸盐相和金属相总量
        sil_sum = x_ + y_ + p.u + p.m + p.n + z_
        metal_sum = a_ + b_ + c_ + d_
        
        if sil_sum <= 1e-15 or metal_sum <= 1e-15:
            return 1e10
        
        X_FeO_sil = x_ / sil_sum
        X_Fe_met = a_ / metal_sum
        
        if X_FeO_sil <= 1e-15 or X_Fe_met <= 1e-15:
            return 1e10
        
        import math
        try:
            IW = 2.0 * math.log10(X_FeO_sil / X_Fe_met)
            return IW
        except (ValueError, ZeroDivisionError):
            return 1e10

    def solve(self, O_bounds: Tuple[float, float] = (-10.0, 10.0)) -> Tuple[float, float]:
        """
        求解最优 O_L，使 IW 匹配目标值。
        
        Returns
        -------
        best_O_L : float
            最优的 O_L 值
        min_residual : float
            最小残差
        """
        result = minimize_scalar(
            self._residual,
            bounds=O_bounds,
            method="bounded",
            options={
                "xatol": self.outer_tol,
                "maxiter": self.max_outer_iter
            }
        )
        
        # ── 如果优化失败，使用缓存的最优解
        if self._cached_result is None:
            print("  Warning: Optimization did not produce valid results. "
                  "Check O_bounds or parameters.")
            return result.x, result.fun
        
        return result.x, result.fun

    def get_final_result(self) -> Tuple[Dict[str, float], float]:
        """
        获取最优 O_L 下的求解结果。
        
        Returns
        -------
        res_dict : dict
            求解结果 (x_, y_, z_, a_, b_, c_, d_, 等)
        IW_final : float
            最终的 IW 值
            
        Raises
        ------
        RuntimeError
            如果 solve() 尚未被调用或无有效缓存
        """
        if self._cached_result is None:
            raise RuntimeError(
                "No cached result. Call solve() first to compute the optimal O_L. "
                "If solve() was called but failed, check error messages above."
            )
        res_dict, IW_final = self._cached_result
        return res_dict, IW_final
#**************************************************************************************************************************************