# accrediff/utils.py
import math
import numpy as np # type: ignore
import pandas as pd  # type: ignore
#**************************************************************************************************************************************
def normalize_max(data):
    """
    对输入数据按最大值归一化，返回归一化后的结果。
    支持 numpy 数组、pandas Series 或 DataFrame。
    """
    
    if isinstance(data, pd.DataFrame):
        return data / data.max()
    elif isinstance(data, pd.Series):
        return data / data.max()
    elif isinstance(data, np.ndarray):
        return data / np.max(data)
    else:
        raise TypeError('仅支持 numpy.ndarray, pandas.Series 或 pandas.DataFrame 类型')
#**************************************************************************************************************************************
def power_law(x, a, b):
    """Power law function."""
    return a * (x ** b)
#**************************************************************************************************************************************
def error_rate(x, y):
    """
    计算误差率
    :param x: 实际值
    :param y: 预测值
    :return: 误差率 %
    """
    return np.abs((x - y) / x) * 100
#**************************************************************************************************************************************
class WeightedECDF:
    """
    加权经验分布函数 F(x)=P(A<=x)，权重需非负。
    若 weights 归一化，F(x) 在 [0,1]；否则自动按总权重归一化。
    """
    def __init__(self, a, weights):
        a = np.asarray(a, dtype=float)
        w = np.asarray(weights, dtype=float)
        # 去除 NaN 和 非有限值
        mask = np.isfinite(a) & np.isfinite(w)
        a, w = a[mask], w[mask]
        if a.size == 0:
            raise ValueError("输入为空")
        if np.any(w < 0):
            raise ValueError("权重必须非负")
        order = np.argsort(a)
        self.a_sorted = a[order]
        self.w_sorted = w[order]
        total = self.w_sorted.sum()
        if total <= 0:
            raise ValueError("权重和必须大于 0")
        self.cdf = np.cumsum(self.w_sorted) / total  # 归一化为 [0,1]

    def __call__(self, x):
        return self.cdf_at(x)

    def cdf_at(self, x):
        x_arr = np.atleast_1d(x).astype(float)
        idx = np.searchsorted(self.a_sorted, x_arr, side='right') - 1
        idx = np.clip(idx, -1, len(self.cdf)-1)
        y = np.where(idx >= 0, self.cdf[idx], 0.0)
        y = np.where(x_arr >= self.a_sorted[-1], 1.0, y)
        return y if np.ndim(x) else float(y[0])

    def quantile(self, q):
        q_arr = np.atleast_1d(q).astype(float)
        if np.any((q_arr < 0) | (q_arr > 1)):
            raise ValueError("q 必须在 [0,1] 区间内")
        idx = np.searchsorted(self.cdf, q_arr, side='left')
        idx = np.clip(idx, 0, len(self.a_sorted)-1)
        res = self.a_sorted[idx]
        return res if np.ndim(q) else float(res[0])
#**************************************************************************************************************************************