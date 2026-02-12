# accrediff/accretion.py
import pandas as pd # type: ignore
import numpy as np # type: ignore
from scipy.optimize import curve_fit # type: ignore
import contextlib
import io
from .utils import normalize_max
from .differentiation import equil_pressure, P_to_T
from .melt_model import MeltScalingModel
#**************************************************************************************************************************************
class CollisionTracer:
    def __init__(self, df, columns=None):
        """
        初始化碰撞追踪器，支持自定义字段映射。

        参数：
        df : pd.DataFrame
            包含碰撞事件的数据表格
        columns : dict 或 None
            字段名映射，例如：
            {
                'time': 'Time',
                'indexi': 'MainID',
                'm_ei': 'MainMass',
                'indexj': 'ImpactorID',
                'm_ej': 'ImpactorMass'
            }
        """
        self.df = df.copy()
        default_cols = {
            'time': 'time',
            'indexi': 'indexi',
            'm_ei': 'm_ei',
            'indexj': 'indexj',
            'm_ej': 'm_ej'
        }
        self.columns = columns or default_cols

        # 校验字段是否存在，并重命名为标准格式
        for key, col in self.columns.items():
            if col not in self.df.columns:
                raise ValueError(f"缺少必要列: '{col}'")
        self.df.rename(columns={v: k for k, v in self.columns.items()}, inplace=True)

        # 为每条碰撞记录添加真实产物 ID（即碰撞后保留的粒子编号）
        #self.df["product_id"] = self.df.apply(self._resolve_product_id, axis=1)
        self.df["product_id"] = self.df.apply(lambda row: int(self._resolve_product_id(row)), axis=1)
    def _resolve_product_id(self, row):
        """
        根据质量和编号规则判断每次碰撞产物的 ID：
        - 质量大者保留；
        - 若质量相等，取 index 较小者。
        """
        if row["m_ei"] > row["m_ej"]:
            return row["indexi"]
        elif row["m_ej"] > row["m_ei"]:
            return row["indexj"]
        else:
            return min(row["indexi"], row["indexj"])

    def trace_full_history(self, particle_id):
        """
        获取粒子的完整碰撞谱系历史：向前（祖先）和向后（后代）

        参数：
        particle_id : int
            要追踪的粒子编号

        返回：
        pd.DataFrame
            含所有相关碰撞记录（按时间排序）
        """
        visited = set()
        full_history = []

        def trace_ancestors(pid):
            if pid in visited:
                return
            visited.add(pid)

            # 查找所有产物是 pid 的碰撞记录
            rows = self.df[self.df["product_id"] == pid]
            full_history.extend(rows.to_dict("records"))

            for _, row in rows.iterrows():
                # 追踪这次碰撞中的两个来源粒子
                trace_ancestors(row["indexi"])
                trace_ancestors(row["indexj"])

        def trace_descendants(pid):
            if pid in visited:
                return
            visited.add(pid)

            # 查找所有该粒子参与、且产物中保留的事件
            rows = self.df[
                ((self.df["indexi"] == pid) | (self.df["indexj"] == pid)) &
                (self.df["product_id"] == pid)
            ]
            full_history.extend(rows.to_dict("records"))

            for _, row in rows.iterrows():
                trace_descendants(row["product_id"])

        # 启动双向谱系追踪
        visited.clear()
        trace_ancestors(particle_id)

        visited.clear()
        trace_descendants(particle_id)

        return pd.DataFrame(full_history).drop_duplicates().sort_values(by="time").reset_index(drop=True)
#**************************************************************************************************************************************
def resolve_product_id(row):
    if row['mi'] > row['mj']:
        return row['indexi']
    elif row['mi'] < row['mj']:
        return row['indexj']
    else:
        return min(row['indexi'], row['indexj'])
#**************************************************************************************************************************************    
class GrowthModel:
    @staticmethod
    def mars_growth(t, tau):
        return np.tanh(t / tau) ** 3

    @staticmethod
    def earth_growth(t, tau):
        return 1 - np.exp(-t / tau)

    @staticmethod
    def mars_growth_with_error(t, tau, err_ratio=0.5):
        err = tau * err_ratio
        upper = GrowthModel.mars_growth(t, tau + err)
        lower = GrowthModel.mars_growth(t, tau - err)
        return upper, lower

    @staticmethod
    def earth_growth_with_error(t, tau, err_ratio=0.2):
        err = tau * err_ratio
        upper = GrowthModel.earth_growth(t, tau + err)
        lower = GrowthModel.earth_growth(t, tau - err)
        return upper, lower
#**************************************************************************************************************************************
class GrowthFitter:
    """
    胚胎生长曲线拟合工具类。
    """

    def __init__(self):
        pass

    def fit_earth_growth_tau(self, E_archive, E_cols):
        """
        对E_archive中每个E_cols列进行按最大值归一化，并用curve_fit拟合wd.GrowthModel.earth_growth，返回每个胚胎的tau_fit列表。
        """
        tau_fits = []
        t_val = E_archive['Time']
        for col in E_cols:
            m_val = E_archive[col]
            m_norm = normalize_max(m_val)
            popt, _ = curve_fit(GrowthModel.earth_growth, t_val, m_norm, p0=[10])
            tau_fit = popt[0]
            tau_fits.append(tau_fit)
        return tau_fits

    def fit_mars_growth_tau(self, M_archive, M_cols):
        """
        对M_archive中每个M_cols列进行按最大值归一化，并用curve_fit拟合GrowthModel.mars_growth，返回每个胚胎的tau_fit列表。
        """
        tau_fits = []
        t_val = M_archive['Time']
        for col in M_cols:
            m_val = M_archive[col]
            m_norm = normalize_max(m_val)
            popt, _ = curve_fit(GrowthModel.mars_growth, t_val, m_norm, p0=[10])
            tau_fit = popt[0]
            tau_fits.append(tau_fit)
        return tau_fits
#**************************************************************************************************************************************    
def build_pl_source_dict(pl_history_dict, df_initial):
    """
    根据每个行星的碰撞历史，构建其来源物质分布（按初始 a 排序，并给出 m_norm）。
    返回: {planet_key: DataFrame[a, e, inc, m_e, m_norm]}
    """
    out = {}
    for key, pl_history in pl_history_dict.items():
        # 收集参与碰撞的粒子ID（i/j 两侧）
        C_id_set = set(pl_history['indexi'].tolist()) | set(pl_history['indexj'].tolist())
        if not C_id_set:
            out[key] = pd.DataFrame(columns=['a','e','inc','m_e','m_norm'])
            continue
        mask = df_initial.index.isin(list(C_id_set))
        pl_source = df_initial.loc[mask, ['a','e','inc','m_e']].sort_values(by='a', ascending=True).copy()
        M_total = pl_source['m_e'].sum()
        if M_total <= 0:
            pl_source['m_norm'] = 0.0
        else:
            pl_source['m_norm'] = pl_source['m_e'] / M_total
        out[key] = pl_source
    return out
#**************************************************************************************************************************************

class ImpactEventProcessor:
    """
    用于处理碰撞事件的分析与标记，包括碰撞平衡条件构建和大事件的冲击模型计算。
    """

    def __init__(self, df, p_ratio=0.6, Mtotal=0.1, ratio=[0.1, 0.5]):
        """
        初始化并构建碰撞平衡条件的DataFrame（df_CE）
        参数:
            df: 原始DataFrame
            p_ratio: 平衡压力的CMB比例
            Mtotal: 质量阈值
            ratio: 冲击比阈值范围
        """
        self.df = df
        self.p_ratio = p_ratio
        self.Mtotal = Mtotal
        self.ratio = ratio
        self.df_CE = self.build_collision_equilibrium()

    def build_collision_equilibrium(self):
        df_CE = pd.DataFrame()
        df_CE[['Time', 'target_id']] = self.df[['Time', 'product_id']].copy()
        df_CE['impactor_id'] = self.df.apply(
            lambda row: int(row['indexi']) if row['indexi'] != row['product_id'] else int(row['indexj']), axis=1
        )
        df_CE['m_target'] = self.df.apply(
            lambda row: max(row['m_ei'], row['m_ej']), axis=1
        )
        df_CE['m_impactor'] = self.df.apply(
            lambda row: min(row['m_ei'], row['m_ej']), axis=1
        )
        df_CE['Mass'] = self.df['Mass']
        df_CE['Impact ratio'] = df_CE['m_impactor'] / df_CE['Mass']
        df_CE['events'] = 'small'
        df_CE.loc[
            (df_CE['Mass'] >= self.Mtotal) & 
            (df_CE['Impact ratio'] >= self.ratio[0]) & 
            (df_CE['Impact ratio'] < self.ratio[1]),
            'events'
        ] = 'global'
        df_CE['P_equil'] = df_CE['m_target'].apply(lambda mass: equil_pressure(mass, self.p_ratio))
        df_CE['T_equil'] = df_CE['P_equil'].apply(lambda P: P_to_T(P))
        df_CE['k_mantle'] = np.nan
        return df_CE

    def update_global_melting_process(self, duration=5, PEF=0.5):
        """
        对df_CE中的'global'事件进行冲击模型计算，并更新相关small事件的属性。
        参数:
            duration: 时间窗口（Myr），默认±5 Myr
        返回:
            df_CE: 更新后的DataFrame
        """
        Global_index = self.df_CE[self.df_CE['events'] == 'global'].index.tolist()
        for i in Global_index:
            Mass, Impact_ratio = self.df_CE.loc[i, ['Mass', 'Impact ratio']].values
            impact = MeltScalingModel(
                Mtotal=10 * Mass,
                gamma=Impact_ratio,
                vel=1.0,
                entropy0=1100,
                impact_angle=45
            )
            with contextlib.redirect_stdout(io.StringIO()):
                impact_res = impact.run_model()
            P_Global = impact_res['max pressure (global model)'][0] * PEF
            k_mantle = impact_res['melt fraction']
            self.df_CE.loc[i, ['P_equil', 'k_mantle']] = [P_Global, k_mantle]

            Timing = self.df_CE.loc[i, 'Time']
            mask = (
                (self.df_CE['target_id'] == self.df_CE.loc[i, 'target_id']) &
                (self.df_CE['Time'].between(Timing, Timing + duration)) &
                (self.df_CE['events'] == 'small')
            )
            '''
            mask = (
                (self.df_CE['target_id'] == self.df_CE.loc[i, 'target_id']) &
                (self.df_CE['Time'].between(Timing - duration, Timing + duration)) &
                (self.df_CE['events'] == 'small')
            )
            '''
            self.df_CE.loc[mask, 'events'] = 'global'
            self.df_CE.loc[mask, 'P_equil'] = P_Global
            self.df_CE.loc[mask, 'T_equil'] = P_to_T(P_Global)
            self.df_CE.loc[mask, 'k_mantle'] = k_mantle
        return self.df_CE

# 用法示例
# processor = ImpactEventProcessor(df)
# df_CE = processor.df_CE
# df_CE = processor.update_global_melting_process()
#**************************************************************************************************************************************
def update_partial_melting_process(df_CE, Mtotal=0.1, ratio=[0.1, 0.5], events='partial'):
    """
    标记满足条件的事件为'partial'，并计算其P_equil和T_equil
    """
    # 满足条件[Mass>=0.1, 0.1<=Impact ratio<0.5的标记为'partial'
    df_CE.loc[(df_CE['Mass'] >= Mtotal) & (df_CE['Impact ratio'] >= ratio[0]) & (df_CE['Impact ratio'] < ratio[1]), 'events'] = events
    Partial_index = df_CE[df_CE['events'] == events].index.tolist()
    for i in Partial_index:
        Mass, Impact_ratio = df_CE.loc[i, ['Mass', 'Impact ratio']].values
        impact = MeltScalingModel(
            Mtotal=10 * Mass,
            gamma=Impact_ratio,
            vel=1.0,
            entropy0=1100,
            impact_angle=45
        )
        # 屏蔽 impact.run_model() 的输出
        with contextlib.redirect_stdout(io.StringIO()):
            impact_res = impact.run_model()
        P_Partial = impact_res['max pressure (melt pool model)'][0]
        df_CE.loc[i, ['P_equil', 'T_equil']] = [P_Partial, P_to_T(P_Partial)]
    return df_CE
# 用法示例：
# df_partial= update_partial_melting_process(df_CE)

#**************************************************************************************************************************************