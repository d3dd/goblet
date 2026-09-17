# GOBLET — Goner-Buffering Lethality Transcriptomic Pipeline

GOBLET is an exploratory computational pipeline for asking whether tumour gene
expression co-varies with the burden of predicted loss-of-function (LOF)
variants in an operational set of putative loss-intolerant genes (Goners).
Candidate expression associations can subsequently be filtered using curated
synthetic-lethal relationships.

## Scientific status

The proof-of-concept analysis uses pancreatic adenocarcinoma (TCGA-PAAD). It
aligns mutation and RNA-expression data for 165 cases and defines each case's
predictor as the number of distinct candidate Goner genes carrying at least one
predicted LOF variant.

No expression association passed the prespecified criteria of Spearman
`rho > 0.3` and Benjamini–Hochberg-adjusted `p < 0.05`. The largest positive
association was PARD6B (`rho = 0.342`, nominal `p = 6.95e-6`, adjusted
`p = 0.121`). Excluding the case with the largest mutation and candidate-LOF
burdens did not produce an FDR-significant result.

BCAT2 and DBT appear in stored output from the legacy exploratory notebook,
but they are **not supported as GOBLET results by the analysis reported in the
Discovery Note**. BCAA metabolism remains a biologically motivated hypothesis
for future testing, not a detected buffering mechanism in this cohort.

GOBLET 1.0 should therefore be understood as an exploratory pipeline and the
TCGA-PAAD analysis as a proof of computational feasibility. It does not yet
constitute a validated detector or establish applicability to other cancers.

## Reproducing the reported PAAD analysis

The analysis and deposited outputs supporting manuscript BIOADV-2026-434 are
in [`bioadv_analysis/`](bioadv_analysis/):

```bash
python bioadv_analysis/run_paad_analysis.py
```

The command uses the deposited inputs in `data/` and `goner_list.txt`. It
normalises raw counts to `log2(CPM + 1)`, tests genes expressed at at least
1 CPM in at least 20% of the aligned cases, applies Benjamini–Hochberg
correction, and repeats the screen after excluding the case with the largest
candidate-LOF burden.

The required Python packages are NumPy, pandas and SciPy. The full legacy
notebook has additional dependencies listed in `requirements.txt`.

## Repository structure

```text
goblet/
├── bioadv_analysis/              # Analysis and outputs reported in the Note
│   ├── run_paad_analysis.py
│   ├── gencode.v36.genes.sorted.tsv
│   └── outputs/
├── GonerFull24Feb2026.ipynb      # Legacy exploratory notebook; see warning
├── goner_list.txt                # Deposited 1,526-gene candidate set
├── requirements.txt              # Legacy notebook dependencies
└── data/
    ├── paad_expr_matrix.xlsx     # 60,660 genes × 183 RNA aliquots
    ├── paad_maaf.xlsx            # TCGA-PAAD somatic-variant table
    ├── Human_Mouse_Common.csv
    ├── HSIAO_Cleaned_paired_genes.csv
    └── gene_sl_gene.csv
```

The expression workbook contains 183 aliquots representing 178 cases. The
reported analysis retains the first deposited aliquot per case and uses the
165 case identifiers shared with the mutation table.

## Legacy notebook

[`GonerFull24Feb2026.ipynb`](GonerFull24Feb2026.ipynb) is retained to document
the development of GOBLET 1.0. Its stored outputs include the earlier BCAT2,
DBT and BCAA-enrichment results. Those outputs use a different burden
construction and should not be interpreted as findings of the current
TCGA-PAAD analysis. The notebook begins with the same warning and directs
readers to `bioadv_analysis/`.

## Inputs

### Candidate Goner list

`goner_list.txt` contains the operational candidate set derived from DepMap
23Q2 Chronos scores: genes with a score below −0.5 in at least 90% of cell
lines. This pan-cancer threshold may not represent PAAD-specific essentiality.

### Mutation data

Predicted LOF consequences comprise frame-shift insertions/deletions, nonsense
mutations, splice-site variants, translation-start-site variants and nonstop
mutations. Multiple variants in the same candidate gene and case count once;
genes outside the candidate set do not contribute to the burden.

### Expression data

The deposited expression workbook contains raw GDC STAR counts. The accompanying
GENCODE v36 mapping in `bioadv_analysis/` restores gene identifiers and symbols
in the deposited row order.

### Synthetic-lethal data

Curated relationships in `data/gene_sl_gene.csv` are used only after an
expression association passes the statistical screen. No expression
association passed that screen in the reported PAAD analysis.

## Limitations

The predictor is sparse and uneven: 137 of 165 cases have burden zero, whereas
one hypermutated case contributes 93 of the 126 candidate-gene/case LOF
observations. The pan-cancer essentiality threshold may miss PAAD-specific
dependencies. Future methods should account for total mutation burden, driver
context, tumour purity and copy number, and should be validated in independent
cohorts.

A null result from this cross-sectional burden–expression test does not exclude
alteration-specific or non-transcriptional buffering, which may fall outside
the tested signal or be obscured by sparse and uneven burden, hypermutation or
other confounding factors.

## Citation

If you use GOBLET, please cite:

> Boman M et al. (2026). Transcriptomic Buffering Detection: A Computational
> Approach to Identifying Essential Gene Dependencies in Pancreatic Cancer.
> *Bioinformatics Advances*, BIOADV-2026-434 (under revision).

## Licence

MIT
