"""CLI: python3 -m hrdness <expr.tsv> <out_prefix> [--platform rnaseq|affy] [--log]

入力: 1列目=Gene Symbol、1行目=サンプルID の発現行列 (遺伝子 × サンプル、全遺伝子)。HGSOC を想定、推奨 50 例以上。
  rnaseq (既定): TPM / FPKM (既に log2 の場合は --log)
  affy: RMA 系の log2 発現 (MAS5 でない方; 既に log2 なので変換しない)。Agilent は未対応
出力1 (ssGSEA, cohort 内 z): ssGSEA_UP / ssGSEA_DN / HRDness_score (=UP−DN) / HRDness_score_z  (50 例未満は参考)
出力2 (基本 = 「2回 z」): log2(FPKM/TPM+1) の 249 遺伝子 → サンプル内 z → 遺伝子ごと z (cohort 内) → TCGA 学習済み LR+RF
        HRDness_probability / HRDness_annotation (HRD-high ≥ 0.65 / HRD-low ≤ 0.35 / それ以外 undetermined)
        n≥50: cohort 内の遺伝子 z (2Z-cohort) / --stats 指定: 保存済み統計の遺伝子 z (2Z-stats, 1 サンプルでも可) /
        どちらも無い (n<50): サンプル内 z のみ (1Z-single, 参考値)。HRDness_mode 列に記録。HRDness_probability_1Z は常に参考出力。
--save-stats: n≥50 の cohort で遺伝子統計を保存 (同一パイプラインの後続サンプル用)。
"""
import argparse
import sys

import numpy as np
import pandas as pd

from .core import annotate, cohort_stats, predict_probability
from .ssgsea import collapse_genes, ssgsea_hrdness

HIGH, LOW = 0.65, 0.35
MIN_COHORT = 50


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hrdness", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("expr"); ap.add_argument("out_prefix")
    ap.add_argument("--platform", choices=["rnaseq", "affy"], default="rnaseq",
                    help="rnaseq: TPM/FPKM → log2(x+1) / affy: RMA 系 log2 をそのまま使う (遺伝子 Symbol に集約済みの行列)")
    ap.add_argument("--log", action="store_true", help="入力が既に log スケール (RNA-seq の log2(x+1) を適用しない)")
    ap.add_argument("--save-stats", metavar="FILE", help="(n≥50) この cohort の遺伝子統計を保存。同一パイプラインの後続サンプル用")
    ap.add_argument("--stats", metavar="FILE", help="保存済み cohort 統計で補正した確率 (*_refadj) も出す。1 サンプルでも可")
    a = ap.parse_args(argv)

    expr = collapse_genes(pd.read_csv(a.expr, sep="\t", index_col=0))
    n = expr.shape[1]
    warns = []
    if n < MIN_COHORT:
        warns.append(f"n={n} < {MIN_COHORT}: cohort 内 z・cohort 補正は参考程度 (補正版は出力しない)")
    plat = a.platform
    if plat == "affy":
        a.log = True                                    # array は既に log2 スケール
        if expr.max().max() > 50:
            sys.exit("ERROR: --platform affy は RMA 系の log2 発現 (値域 ≈2–15) を想定。最大値が 50 超: 非 log の値では？")
        if expr.min().min() < 1.0:
            warns.append(f"最小値 {expr.min().min():.2f}: MAS5 の log2 (負値あり) の可能性。TCGA 参照は RMA 系 (≈2.6–13)。"
                         "2Z で影響は小さい(GSE26712: Spearman 0.94)が、基本は RMA 系を使うこと")
    if a.log and expr.max().max() > 50:
        warns.append("--log 指定だが最大値が 50 を超える: 非 log の TPM/FPKM なら --log を外す")

    vals = expr.to_numpy(dtype=float)
    if plat == "rnaseq" and np.allclose(vals, np.round(vals)) and vals.max() > 1000 and not a.log:
        sys.exit("ERROR: 整数値の発現量 (raw counts?)。遺伝子長で補正された TPM または FPKM を入力してください "
                 "(counts は /rnaseq-normalize スキル等で TPM/FPKM に変換)")
    out = pd.DataFrame(index=expr.columns)
    if expr.shape[0] >= 5000:                                     # 出力1 (全遺伝子行列が必要)
        s = ssgsea_hrdness(expr)
        out["ssGSEA_UP"], out["ssGSEA_DN"], out["HRDness_score"] = s["UP"], s["DN"], s["SIG"]
        out["HRDness_score_z"] = (s["SIG"] - s["SIG"].mean()) / s["SIG"].std(ddof=1) if n > 1 else np.nan
        out["HRDness_score_reliable"] = n >= MIN_COHORT
        if n >= MIN_COHORT:
            out["HRDness_score_group"] = np.where(s["SIG"] >= s["SIG"].median(), "high", "low")
    else:
        warns.append(f"遺伝子数 {expr.shape[0]}: 全遺伝子行列でないため出力1 (ssGSEA) をスキップ")

    e2 = expr if a.log else np.log2(expr.clip(lower=0) + 1)   # 出力2: log2(FPKM/TPM+1) → 249 遺伝子
    p1 = predict_probability(e2, plat)                     # 1Z: サンプル内 z のみ (単一サンプル可・絶対値は cohort でずれる)
    if a.stats:                                                # 2Z (保存済み統計): サンプル内 z → 保存済み cohort の遺伝子 z
        st = pd.read_csv(a.stats, sep="\t", index_col=0)
        p, mode = predict_probability(e2, plat, stats=st), "2Z-stats"
        if p["n_genes"].iloc[0] < 240:
            warns.append(f"stats と一致した遺伝子が {p['n_genes'].iloc[0]}/249")
    elif n >= MIN_COHORT:                                      # 2Z (cohort): サンプル内 z → cohort 内の遺伝子 z  (基本)
        p, mode = predict_probability(e2, plat, cohort_adjust=True), "2Z-cohort"
    else:
        p, mode = p1, "1Z-single"
        warns.append(f"n={n} < {MIN_COHORT} で --stats も無い: 1Z (サンプル内 z のみ) の参考値。cohort でずれる可能性があるため "
                     "同一パイプラインの 50 例以上で --save-stats した統計を --stats で渡すこと")
    out["LR_prob"], out["RF_prob"] = p["LR_prob"], p["RF_prob"]
    out["HRDness_probability"] = p["HRDness_probability"]
    out["HRDness_annotation"] = annotate(p["HRDness_probability"], HIGH, LOW)
    out["HRDness_mode"] = mode
    out["HRDness_probability_1Z"] = p1["HRDness_probability"]
    out["n_signature_genes"] = p["n_genes"]
    mp = p1["HRDness_probability"].mean()
    if n >= 10 and (mp < 0.3 or mp > 0.7):
        warns.append(f"1Z の平均 probability={mp:.2f} が極端 (前処理/プラットフォーム差)。主出力は 2Z なので影響は小さいが 1Z 列は参考にしない")
    if a.save_stats:
        if n < MIN_COHORT:
            warns.append(f"n={n} < {MIN_COHORT}: 統計は不安定なので --save-stats を実行しない")
        else:
            cohort_stats(e2, plat).to_csv(a.save_stats, sep="\t"); print(f"saved cohort stats: {a.save_stats}")
    out.index.name = "sample"
    out.round(4).to_csv(f"{a.out_prefix}_HRDness.tsv", sep="\t")
    cnt = out["HRDness_annotation"].value_counts().to_dict()
    print(f"n={n} genes={p['n_genes'].iloc[0]}/249 | annotation: {cnt}")
    for w in warns:
        print("WARNING:", w, file=sys.stderr)


if __name__ == "__main__":
    main()
