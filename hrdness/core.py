"""出力2: サンプル内 249 遺伝子 z-score → TCGA 学習済み分類器 (LR + RF) → HRDness probability。

- 特徴量はサンプルごとの z-score (入力遺伝子 ∩ 参照遺伝子の範囲で標準化) → 1 サンプルから計算可能。
- 学習は TCGA-OV の原発 (-01) サンプルのみ、HRD-high/low の極端群 (症例単位ラベル)。RNA-seq (GDC log2(FPKM+1)) のみ対応。
- probability = LR と RF(5 seed 平均) の確率の平均。
"""
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "data/reference_TCGA"
# rnaseq: GDC log2(FPKM+1) / affy: Xena HT_HG-U133A (RMA 系 log2, 値域 ≈2.6–13)。
# Agilent (TCGA_OV_Agilent_249genes.tsv) は外部での絶対確率のずれが大きく保留 (追加すれば動くが未検証)。
PLATFORMS = {"rnaseq": "TCGA_OV_RNAseq_log2fpkm_249genes.tsv", "affy": "TCGA_OV_Affy_249genes.tsv"}


def load_signature(sig_dir: Path = ROOT / "data/signature"):
    up = [l.strip() for l in open(sig_dir / "TCGA_OV_HRDness_SIG_UP.txt") if l.strip()]
    dn = [l.strip() for l in open(sig_dir / "TCGA_OV_HRDness_SIG_DN.txt") if l.strip()]
    return up, dn


def sample_z(X: pd.DataFrame) -> pd.DataFrame:
    """X: サンプル × 遺伝子 → サンプルごと (行ごと) に遺伝子間で z-score。"""
    return X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=1), axis=0)


def gene_z(Z: pd.DataFrame) -> pd.DataFrame:
    """サンプル × 遺伝子 → 遺伝子ごと (列ごと) に集団内で z-score。"""
    return (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)


@lru_cache(maxsize=None)
def _reference(platform: str):
    ref = pd.read_csv(REF / PLATFORMS[platform], sep="\t", index_col=0)
    lab = pd.read_csv(REF / "TCGA_OV_HRD_label_rnaseq_case.tsv", sep="\t", index_col=0)
    c2s = {}
    for c in sorted(x for x in ref.columns if x[13:15] == "01"):
        c2s.setdefault(c[:12], c)
    cases = [c for c in lab.index if c in c2s]
    X = ref[[c2s[c] for c in cases]].T
    y = lab.loc[cases, "HRD_label"].astype(int).values
    return X, y


@lru_cache(maxsize=None)
def fit_reference(platform: str, genes: tuple, adjust: bool = False):
    """指定遺伝子 (入力と参照の共通集合) で TCGA 参照から LR / RF を学習。

    adjust=True: サンプル内 z の後、遺伝子ごとに TCGA 内で再標準化した特徴量で学習 (cohort 補正モード用)。
    """
    X, y = _reference(platform)
    Xz = sample_z(X[list(genes)])
    if adjust:
        Xz = gene_z(Xz)
    lr = LogisticRegression(max_iter=1000, random_state=3).fit(Xz, y)
    rfs = [RandomForestClassifier(random_state=s).fit(Xz, y) for s in range(3, 8)]
    return lr, rfs


def cohort_stats(expr: pd.DataFrame, platform: str = "rnaseq") -> pd.DataFrame:
    """cohort の遺伝子ごと統計 (サンプル内 z の後の平均 / SD)。同一パイプラインの後続サンプルの補正用に保存する。"""
    ref_genes = list(_reference(platform)[0].columns)
    genes = [g for g in ref_genes if g in expr.index and expr.loc[g].notna().all()]
    Z = sample_z(expr.loc[genes].T.astype(float))
    return pd.DataFrame({"mean": Z.mean(axis=0), "sd": Z.std(axis=0, ddof=1)})


def predict_probability(expr: pd.DataFrame, platform: str = "rnaseq", min_genes: int = 200, cohort_adjust: bool = False,
                        stats: pd.DataFrame = None) -> pd.DataFrame:
    """expr: 遺伝子 × サンプル (log スケール; RNA-seq は log2(TPM/FPKM+1))。

    cohort_adjust=True: 入力 cohort 内の遺伝子統計で補正 (n≥50 推奨)。stats を渡すと、保存済みの別 cohort 統計
    (cohort_stats の出力; 同一パイプライン) で補正する → 1 サンプルでも補正できる。
    返り値: サンプル × [LR_prob, RF_prob, HRDness_probability, n_genes]。
    """
    ref_genes = list(_reference(platform)[0].columns)
    genes = tuple(g for g in ref_genes if g in expr.index and expr.loc[g].notna().all()
                  and (stats is None or g in stats.index))
    cohort_adjust = cohort_adjust or stats is not None
    if len(genes) < min_genes:
        raise ValueError(f"シグネチャ遺伝子の一致が {len(genes)} 個 (< {min_genes})。遺伝子名(Symbol)/行列の向きを確認")
    lr, rfs = fit_reference(platform, genes, cohort_adjust)
    Xz = sample_z(expr.loc[list(genes)].T.astype(float))
    if stats is not None:                   # 保存済み cohort 統計で補正 (単一サンプル可)
        Xz = (Xz - stats.loc[list(genes), "mean"]) / stats.loc[list(genes), "sd"]
    elif cohort_adjust:                     # 入力 cohort 内で遺伝子ごとに再標準化 (n が十分ある場合のみ使うこと)
        Xz = gene_z(Xz)
    p_lr = lr.predict_proba(Xz)[:, 1]
    p_rf = np.mean([m.predict_proba(Xz)[:, 1] for m in rfs], axis=0)
    return pd.DataFrame({"LR_prob": p_lr, "RF_prob": p_rf, "HRDness_probability": (p_lr + p_rf) / 2,
                         "n_genes": len(genes)}, index=expr.columns)


def annotate(prob: pd.Series, hi: float, lo: float) -> pd.Series:
    return pd.Series(np.where(prob >= hi, "HRD-high", np.where(prob <= lo, "HRD-low", "undetermined")), index=prob.index)
