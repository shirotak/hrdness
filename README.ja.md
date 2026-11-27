# hrdness

[English](README.md) | **日本語**

TCGA-OV 由来の 249 遺伝子 HRDness signature を使い、遺伝子発現から高異型度漿液性卵巣癌（HGSOC）の
**HRDness**（相同組換え修復欠損様の発現パターン）を推定します。

> Signature: Takamatsu *et al.*, *Br J Cancer* 2023;128:1095–1104, doi:[10.1038/s41416-022-02122-9](https://doi.org/10.1038/s41416-022-02122-9)
> （UP 173 遺伝子 + DN 75 遺伝子）。`hrdness` は、この signature を新規 cohort に適用するための再実装です。

**学術研究目的（非商用）に限ります。臨床診断には使用できません。**

## できること

| 出力 | 手法 | 用途 |
|---|---|---|
| **1. HRDness score**（連続値） | UP / DN 遺伝子セットの ssGSEA（rank 正規化, weight 0.75）。score = UP − DN。元の R ライブラリの Python 移植（R の出力と一致）。 | **cohort 内**でのサンプル間比較。n < 50 は参考値。 |
| **2. HRDness probability**（連続値）+ **annotation**（`HRD-high` / `HRD-low` / `undetermined`） | 「2回 z」（2Z）: 249 遺伝子の log2 発現 → **サンプルごと**に z-score → **遺伝子ごと**（サンプル間）に z-score → TCGA-OV（HRD-high vs HRD-low の極端群）で学習したロジスティック回帰 + ランダムフォレスト。 | cohort に対する相対的なサンプルごとの判定。保存した cohort 統計を使えば **1 サンプル**でも可能。 |

annotation の閾値: probability ≥ 0.65 → `HRD-high`、≤ 0.35 → `HRD-low`、それ以外は `undetermined`。

対応プラットフォーム: **RNA-seq**（既定）と **Affymetrix**（`--platform affy`）。Agilent は未対応です。

## 判定を使う前に必ずお読みください

2Z 法では、annotation は**入力した cohort に対する相対的な判定**になります。

- `HRD-high` は「*この* cohort の他のサンプルより HRD 的」という意味で、「HRD 陽性」ではありません。high の割合は TCGA の学習データの混合比（約 50%）に引き寄せられ、
  **有病率の推定値ではありません**。HRD スコアの閾値（例: 42）や BRCA / HRD 検査の代わりにもなりません。
- **入力は HGSOC に限ります。** 他の組織型、正常組織、腫瘍純度が非常に低い・間質成分が多いサンプルは、遺伝子ごとの統計を歪め、判定がずれます。
- **50 サンプル以上を推奨します**（実務上の目安で、実験的に決めた閾値ではありません）。それ未満では、サンプル内 z のみのモード（`1Z-single`）になります。
  このモードの絶対確率は cohort やプラットフォームによってずれるため、参考値としてください。
- **cohort は対象集団を代表している必要があります。** BRCA 変異例のみ、HRD-high のみ、奏効例のみなどの cohort は、無理に約 50% ずつに割り振られます。
  複数のバッチ・施設が混在する場合は、別々に解析するか補正してください。
- 1 サンプルの判定は、**同一パイプライン**で処理した 50 サンプル以上から保存した統計に対してのみ可能です
  （[1 サンプルの判定](#1-サンプルの判定)を参照）。別のプラットフォームやパイプラインの統計は流用できません。

## インストール

Python 3 と `numpy`、`pandas`、`scipy`、`scikit-learn` が必要です（R も ComBat も不要）。
Python 3.11、numpy 2.4、pandas 2.3、scipy 1.17、scikit-learn 1.8 で開発・動作確認しています。

```bash
git clone <このリポジトリ>
cd hrd_gexp_signatures
pip install -r requirements.txt
python3 -m hrdness --help
```

TCGA の参照テーブルを `data/` から読み込むため、リポジトリ内から（またはリポジトリと同じ場所に置いて）実行してください。

## 入力

タブ区切りの行列（**遺伝子 × サンプル**）:

- 1 列目 = HGSOC の**遺伝子シンボル**（hg38 / HGNC）、1 行目 = サンプル ID
- 249 遺伝子だけでなく**全遺伝子**（ssGSEA に必要。5000 遺伝子未満の場合は出力 1 をスキップ）
- 重複するシンボルは平均発現が最大の行に集約し、欠損値は行の中央値で補完

| プラットフォーム | 入力するもの |
|---|---|
| RNA-seq（`--platform rnaseq`, 既定） | **TPM または FPKM**（線形スケール）。ツールが log2(x + 1) を適用します。すでに log2(x + 1) の場合は `--log` を指定。raw counts（整数）は拒否されます。先に正規化してください。 |
| Affymetrix（`--platform affy`） | 遺伝子シンボルに集約済みの **RMA 系 log2** 発現。log2(x + 1) は**適用しないでください**。MAS5 の値（負値あり）は警告、50 超の値は拒否されます。 |

サンプル内 z-score を取るため、TPM と FPKM の違い、定数倍、pseudocount（0.01〜1）の影響は小さいです。

## 使い方

```bash
python3 -m hrdness <expr.tsv> <out_prefix> [--platform rnaseq|affy] [--log] [--save-stats FILE] [--stats FILE]
```

例（RNA-seq TPM、100 サンプル）:

```bash
python3 -m hrdness my_cohort_tpm.tsv results/my_cohort --save-stats results/my_cohort_stats.tsv
# -> results/my_cohort_HRDness.tsv, results/my_cohort_stats.tsv
```

### 1 サンプルの判定

1. 50 サンプル以上の cohort を `--save-stats` 付きで実行します（同じプラットフォーム、同じ定量パイプライン）。
2. 新しいサンプルを、その統計に対して判定します:

```bash
python3 -m hrdness new_sample_tpm.tsv results/new_sample --stats results/my_cohort_stats.tsv
```

結果は、そのサンプルを参照 cohort の中に入れて実行した場合と同一です（`HRDness_mode = 2Z-stats`）。
判定は参照 cohort に対する相対的なものです。

## 出力: `<out_prefix>_HRDness.tsv`

| 列 | 意味 |
|---|---|
| `ssGSEA_UP`, `ssGSEA_DN` | UP / DN 遺伝子セットの ssGSEA enrichment score |
| `HRDness_score` | UP − DN（出力 1） |
| `HRDness_score_z` | `HRDness_score` の cohort 内 z-score |
| `HRDness_score_reliable` | n ≥ 50 のとき `True` |
| `HRDness_score_group` | cohort の中央値で二分した `high` / `low`（n ≥ 50 のみ） |
| `LR_prob`, `RF_prob` | ロジスティック回帰 / ランダムフォレスト（5 seed の平均）の確率 |
| `HRDness_probability` | `LR_prob` と `RF_prob` の平均（出力 2） |
| `HRDness_annotation` | `HRD-high`（≥ 0.65）/ `HRD-low`（≤ 0.35）/ `undetermined` |
| `HRDness_mode` | `2Z-cohort`（n ≥ 50）、`2Z-stats`（`--stats`）、`1Z-single`（n < 50 で統計なし; 参考値） |
| `HRDness_probability_1Z` | サンプル内 z のみで求めた確率。参考値で、cohort によってずれる |
| `n_signature_genes` | 入力に含まれていた signature 遺伝子数（249 のうち） |

警告（n が小さい、スケールが不自然、遺伝子の一致数が少ない、1Z の平均確率が極端）は stderr に出力されます。

## 仕組み

1. 入力と TCGA 参照の両方にある 249 遺伝子を取り出します（200 遺伝子以上が必要）。
2. それらの遺伝子間で、サンプルごとに z-score を取ります。
3. 遺伝子ごとに、サンプル間で z-score を取ります（cohort の統計、または `--stats` で保存した統計）。
4. TCGA-OV の原発腫瘍（`-01` サンプル、1 症例 1 サンプル; HRD-high vs HRD-low の極端群）を**同じ変換**にかけて学習した LR + RF で判定します。
   RNA-seq の参照は GDC の log2(FPKM + 1)、Affy の参照は Xena の `HT_HG-U133A` です。

参照の HRD ラベル: HRD-high = HRD スコアと SigMA（signature 3）がともに高い、HRD-low = ともに低い。不一致の症例は除外
（閾値は BRCA1/2 LOH に対する ROC から決定）。詳細は論文を参照してください。

## 性能と注意点

- **TCGA 内部 CV（5×5）**: accuracy ≈ 0.88、AUC ≈ 0.94。これは**楽観的**な値です。249 遺伝子は、同じ TCGA の HRD-high/low サンプルから DEG として選ばれているためです。
- **独立 cohort（RNA-seq、4 cohort、TCGA のみで学習）**: cohort 内の HRD スコアとの Spearman ≈ 0.53、HRD スコアが中央値以上かどうかの AUC ≈ 0.80。識別力は中程度です。
  low 判定でも腫瘍が HR 正常であるとは限りません（これらの cohort では `HRD-low` 判定の 25〜45% が HRD スコア < 42 でした）。
- 元の ComBat-to-TCGA + RF/SVM/KNN/LR 方式（再実装）と比べて、2Z はこれらの cohort で同等以上で、`undetermined` が少なくなりました。
  探索的な比較です（4 cohort、同一の TCGA ラベル定義）。
- 分類器（SVM、KNN、MLP、勾配ブースティング）を追加しても、外部での性能は改善しませんでした。
- **Affymetrix** は動作確認のみです（RMA / MAS5 入力、1 サンプルでの再現）。独立した HRD 指標に対する精度は**検証していません**。
  TCGA の Affy 参照の Xena での正規化方法は、値域（≈ 2.6〜13）からの推定です。
- 統計はプラットフォーム固有です: `--stats` ファイルは、作成したプラットフォームとパイプラインに固有です。

## リポジトリ構成

```
hrdness/              パッケージ（core.py: 2Z + LR/RF, ssgsea.py: ssGSEA, __main__.py: CLI）
data/signature/       249 遺伝子 signature（UP / DN のリスト, gmt）
data/reference_TCGA/  TCGA-OV の参照発現（249 遺伝子）と HRD ラベル
requirements.txt      Python の依存パッケージ
NOTICE.md             サードパーティの表示（ssGSEA の license、TCGA データ）
```

## 引用

このツールを使用する場合は、signature の論文を引用してください:

> Takamatsu S, Yoshihara K, Baba T, Shimada M, Yoshida H, Kajiyama H, Oda K, Mandai M, Okamoto A, Enomoto T, Matsumura N.
> Prognostic relevance of HRDness gene expression signature in ovarian high-grade serous carcinoma; JGOG3025-TR2 study.
> *Br J Cancer* 2023;128(6):1095–1104. doi:[10.1038/s41416-022-02122-9](https://doi.org/10.1038/s41416-022-02122-9)

## AI の利用

このソフトウェアは、Claude Code を通じて Claude（Anthropic）の支援を受けて開発しました（コードの作成、評価スクリプト、ドキュメントの草稿を含む）。
作者（Shiro Takamatsu）が作業を指揮して設計上の判断を行い、公開内容に責任を負います。
Claude の支援を受けたコミットには `Co-Authored-By` トレーラーを付けています。

## ライセンス

**学術・非商用利用に限ります。** コードは
[PolyForm Noncommercial License 1.0.0](LICENSE) で公開しています。研究、教育、および非営利団体・公的研究機関・教育機関・政府機関による利用は許可されますが、
**著作権者の別途書面による許諾なしの商用利用は許可されません。**
これはソース公開型のライセンスで、OSI 承認のオープンソースライセンスではありません。

- `hrdness/ssgsea.py` は Broad Institute の ssGSEA コード（BSD 系 license）の Python 移植です。license の全文は
  [NOTICE.md](NOTICE.md) に掲載しています。手法: Barbie *et al.*, *Nature* 2009;462:108–112。
- `data/reference_TCGA/` の参照データは、TCGA のオープンアクセスデータ（GDC および UCSC Xena）に由来します。参加者を特定しようとしないでください。
  論文等では TCGA / GDC / UCSC Xena をデータ元として明記してください（[NOTICE.md](NOTICE.md) を参照）。
