# accrediff/chemistry.py
import numpy as np # type: ignore
import re # type: ignore
from typing import Optional 

from .constant import Elements
#**************************************************************************************************************************************
class MolarMassCalculator:
    """简单化合物摩尔质量计算器（不含括号/点号）"""

    # 项目唯一的默认原子量数据源。实例仍持有副本，避免调用者修改
    # ``mass_table`` 时污染 ``constant.Elements``。
    DEFAULT_MASS_TABLE = Elements

    _TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")  # 元素符号+可选数字

    def __init__(self, mass_table: Optional[dict] = None):
        """
        mass_table: dict 可选，元素符号 → 摩尔质量。
        如果不给，则使用 DEFAULT_MASS_TABLE。
        """
        if mass_table is not None:
            self.mass_table = mass_table
        else:
            self.mass_table = self.DEFAULT_MASS_TABLE.copy()

    def molar_mass(self, formula: str) -> float:
        """返回化学式的摩尔质量"""
        formula = formula.strip()
        if not formula:
            raise ValueError("空化学式")

        total = 0.0
        i = 0
        while i < len(formula):
            m = self._TOKEN.match(formula, i)
            if not m:
                raise ValueError(f"无法解析: '{formula[i:]}' in '{formula}'")

            el, num = m.groups()
            if el not in self.mass_table:
                raise ValueError(f"未知元素: {el}")

            n = int(num) if num else 1
            total += self.mass_table[el] * n
            i = m.end()

        return total

    def breakdown(self, formula: str) -> dict:
        """返回化学式的元素个数分解，例如 SiO2 -> {"Si":1, "O":2}"""
        formula = formula.strip()
        if not formula:
            raise ValueError("空化学式")

        counts = {}
        i = 0
        while i < len(formula):
            m = self._TOKEN.match(formula, i)
            if not m:
                raise ValueError(f"无法解析: '{formula[i:]}' in '{formula}'")

            el, num = m.groups()
            if el not in self.mass_table:
                raise ValueError(f"未知元素: {el}")

            n = int(num) if num else 1
            counts[el] = counts.get(el, 0) + n
            i = m.end()

        return counts
#**************************************************************************************************************************************    
class KDCalculator:
    def __init__(self):
        # 参数格式：{元素: [(a, b, c) for P<5, (a, b, c) for P>=5]}
        raw_params = {
            'Si': [(1.3, -11400, -430), (1.3, -13500, 0)],  # Fisher et al. (2017)
            'O':  [(0.6, -2500, -240), (0.6, -3800, 22)],   # Fisher et al. (2017)
            'Ni': [(0.46, 3400, -200), (0.46, 2700, -61)],  # Fisher et al. (2017)
            'Co': [(0.36, 1800, -92), (0.36, 1500, -33)],   # Fisher et al. (2017)
            'Cr': [(-0.3, -2900, 130), (-0.3, -2200, 0)],   # Fisher et al. (2017)
            'V': [(-1.5, -2200, -17), (-1.5, -2300, 9)],    # Fisher et al. (2017)
            'Nb': [(0, -10314, -67), (0, -10314, -67)],     # Huang et al. (2020)
            'Ta': [(-2.71, -9270, 0), (-2.71, -9270, 0)],   # Huang et al. (2020)
            'Mo': [(4.1, -5563, -215), (4.1, -5563, -215)], # Huang et al. (2021)
            'W': [(1.8, -6942, -84), (1.8, -6942, -84)]     # Huang et al. (2021)

        }
        # 转换为 float
        self.params = {k: [(float(a), float(b), float(c)), (float(d), float(e), float(f))]
                       for k, [(a, b, c), (d, e, f)] in raw_params.items()}

    def add_element(self, element, params_lowP, params_highP):

        """添加新元素的参数（自动转为float）"""
        if len(params_lowP) != 3 or len(params_highP) != 3:
            raise ValueError("参数必须为3个元素的元组")
        self.params[element] = [
            (float(params_lowP[0]), float(params_lowP[1]), float(params_lowP[2])),
            (float(params_highP[0]), float(params_highP[1]), float(params_highP[2]))
        ]

    def get_KD(self, element, P, T):
        if element not in self.params:
            raise ValueError(f"未找到元素 {element} 的参数")
        idx = 0 if P < 5 else 1
        a, b, c = self.params[element][idx]
        log_KD = a + b / T + c * P / T
        return 10 ** log_KD
#**************************************************************************************************************************************    
def compute_core_ratio(**kwargs):
    """
    计算地核/地幔质量比，支持用关键字参数（如**res_dict, Mg_total=..., Al_total=..., Ca_total=..., Elements=...）传参
    """
    Elements = kwargs['Elements']
    Mg_total = kwargs['Mg_total']
    Al_total = kwargs['Al_total']
    Ca_total = kwargs['Ca_total']
    # 其余参数从kwargs获取
    x_ = kwargs['x_']
    y_ = kwargs['y_']
    z_ = kwargs['z_']
    a_ = kwargs['a_']
    b_ = kwargs['b_']
    c_ = kwargs['c_']
    d_ = kwargs['d_']

    mass_m = (
        x_ * Elements['FeO'] +
        y_ * Elements['NiO'] +
        z_ * Elements['SiO2'] +
        Mg_total * Elements['MgO'] +
        (Al_total/2) * Elements['Al2O3'] +
        Ca_total * Elements['CaO']
    )
    mass_c = (
        a_ * Elements['Fe'] +
        b_ * Elements['Ni'] +
        c_ * Elements['Si'] +
        d_ * Elements['O']
    )
    return mass_c / (mass_m + mass_c)
#**************************************************************************************************************************************
def compute_IW(x_=None, a_=None, y_=None, b_=None, z_=None, c_=None, d_=None, u=None, m=None, n=None, **kwargs):
    """
    计算IW值，支持通过res_dict等字典方式传参
    """
    # 支持从kwargs补充参数
    x_ = x_ if x_ is not None else kwargs.get('x_')
    a_ = a_ if a_ is not None else kwargs.get('a_')
    y_ = y_ if y_ is not None else kwargs.get('y_')
    b_ = b_ if b_ is not None else kwargs.get('b_')
    z_ = z_ if z_ is not None else kwargs.get('z_')
    c_ = c_ if c_ is not None else kwargs.get('c_')
    d_ = d_ if d_ is not None else kwargs.get('d_')
    u  = u  if u  is not None else kwargs.get('u')
    m  = m  if m  is not None else kwargs.get('m')
    n  = n  if n  is not None else kwargs.get('n')

    X_FeO  = x_ / (x_ + y_ + z_ + u + m + n)
    #X_MW_FeO = 1.148 * X_FeO + 1.319 * (X_FeO ** 2)
    X_Fe = a_ / (a_ + b_ + c_ + d_)
    return 2 * np.log10(X_FeO / X_Fe)
#**************************************************************************************************************************************
def compute_IW_mol(IW_15, Elements):
    """
    计算IW_15条件下各组分的摩尔分率（以MgO为基准），返回两个字典：
    - IW_15_mol: 各组分摩尔数
    - IW_15_mol_mg: 各组分摩尔分率（以MgO为基准）
    """
    IW_15_mol, IW_15_mol_mg = {}, {}
    Silicate = ['FeO', 'NiO', 'SiO2', 'MgO', 'Al2O3', 'CaO']
    Metal = ['Fe', 'Ni', 'Si', 'O']
    core = IW_15['core']  # core占bulk的质量比

    for i in IW_15.keys():
        if i in Silicate:
            IW_15_mol[i] = IW_15[i] / Elements[i] * (1 - core)
        elif i in Metal:
            IW_15_mol[i] = IW_15[i] / Elements[i] * core
        else:
            pass  # 如果元素不在列表中，可以选择跳过或处理

    mg_mol = IW_15_mol.get('MgO', 1)  # 防止MgO缺失导致除零
    for i in IW_15_mol.keys():
        if i in Elements:
            IW_15_mol_mg[i] = IW_15_mol[i] / mg_mol
        else:
            print(f"元素 {i} 不在元素列表中，无法计算摩尔分率。")
    return IW_15_mol, IW_15_mol_mg
#**************************************************************************************************************************************
