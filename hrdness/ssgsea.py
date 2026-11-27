# NOTICE: Python port of the Broad Institute ssGSEA projection code (R; Pablo Tamayo, Chet Birger et al.; ssGSEA module for GenePattern,
# BSD-style license, Copyright (c) 2003-2020 The Broad Institute, Inc., Massachusetts Institute of Technology and Regents of the
# University of California). The full license text and disclaimer are reproduced in NOTICE.md at the repository root.
# Method: Barbie et al., Nature 2009;462:108-112. このモジュールは上記コードの処理を Python に移植したもの。
"""出力1: ssGSEA による HRDness score (UP − DN)。元解析の R ライブラリ (ssGSEAProjection.Library) の Python 移植。

設定は元解析と同一: sample.norm.type="rank", weight=0.75, combine.mode="combine.add" (SIG = ES_UP − ES_DN), min.overlap=5。
サンプル内ランクにのみ依存するので単調変換 (TPM / FPKM / log2) に不変。ただし遺伝子セットが全遺伝子中のランクで
評価されるため、249 遺伝子だけでなく**全遺伝子の行列**が必要。
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .core import load_signature

WEIGHT = 0.75
MIN_OVERLAP = 5


def _es(order: np.ndarray, rank_vals: np.ndarray, in_set: np.ndarray, weight: float) -> float:
    """1 サンプル・1 遺伝子セットの ES (Project.to.GeneSet と同一)。order: 発現降順の行 index。"""
    n = len(order)
    x = rank_vals[order]
    ranked = (x - x.mean()) / x.std(ddof=1)
    ind = np.flatnonzero(in_set[order])                      # 0-based
    nh = len(ind); nm = n - nh
    r = np.abs(ranked[ind]) ** weight
    up = r / r.sum()
    ind1 = ind + 1                                           # R の 1-based
    gaps = np.append(ind1 - 1, n) - np.insert(ind1, 0, 0)
    down = gaps / nm
    res = np.cumsum(np.append(up, up[-1]) - down)
    valleys = res[:nh] - up
    gaps = gaps + 1
    area = np.append(valleys, 0) * gaps + 0.5 * (np.insert(res[:nh], 0, 0) - np.append(valleys, 0)) * gaps
    return float(area.sum())


def collapse_genes(expr: pd.DataFrame) -> pd.DataFrame:
    """重複 Gene Symbol は平均発現が最大の行を採用 (元解析の gs_max 相当)。NA は行中央値で補完。"""
    expr = expr[expr.index.notna() & (expr.index.astype(str) != "")].astype(float)
    if expr.index.duplicated().any():
        o = expr.mean(axis=1).sort_values(ascending=False, kind="stable").index
        expr = expr.loc[o]
        expr = expr[~expr.index.duplicated(keep="first")]
    if expr.isna().any().any():
        expr = expr.apply(lambda r: r.fillna(r.median()), axis=1)
    return expr


def ssgsea_hrdness(expr: pd.DataFrame, weight: float = WEIGHT) -> pd.DataFrame:
    """expr: 遺伝子(Symbol) × サンプル、全遺伝子。返り値: サンプル × [UP, DN, SIG]。"""
    expr = collapse_genes(expr)
    genes = expr.index.to_numpy()
    up_g, dn_g = load_signature()
    ng = len(genes)
    sets = {}
    for k, g in {"UP": up_g, "DN": dn_g}.items():
        m = np.isin(genes, g)
        if m.sum() < MIN_OVERLAP:
            raise ValueError(f"{k} シグネチャ遺伝子の一致が {m.sum()} 個 (< {MIN_OVERLAP})。遺伝子名(Symbol)を確認")
        sets[k] = m
    out = {k: np.zeros(expr.shape[1]) for k in sets}
    M = expr.to_numpy(dtype=float)
    for j in range(M.shape[1]):
        rv = 10000.0 * rankdata(M[:, j], method="average") / ng
        order = np.argsort(-rv, kind="stable")
        for k, m in sets.items():
            out[k][j] = _es(order, rv, m, weight)
    res = pd.DataFrame(out, index=expr.columns)
    res["SIG"] = res["UP"] - res["DN"]
    res.attrs["n_up"] = int(sets["UP"].sum()); res.attrs["n_dn"] = int(sets["DN"].sum())
    return res
