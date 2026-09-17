# GOBLET 1.0 TCGA-PAAD analysis

This directory reproduces the analysis reported for manuscript
BIOADV-2026-434. It analyses the deposited GOBLET 1.0 inputs and is separate
from the development of GOBLET 2.0.

The analysis:

1. assigns the 60,660 expression rows using the deposited GENCODE v36 mapping;
2. retains the first deposited expression aliquot per case and the 165 case
   identifiers shared with the mutation table;
3. counts distinct predicted LOF genes only within the deposited 1,526-gene
   operational candidate set;
4. transforms raw STAR counts to `log2(CPM + 1)`;
5. tests genes expressed at at least 1 CPM in at least 20% of cases using
   Spearman correlation and Benjamini–Hochberg correction; and
6. repeats the screen after excluding the case with the largest mutation and
   candidate-LOF burdens.

No transcript passes the prespecified FDR threshold of 0.05 in either screen.

## Run

From the GOBLET repository root:

```bash
python bioadv_analysis/run_paad_analysis.py
```

Required deposited inputs are:

- `goner_list.txt`
- `data/paad_expr_matrix.xlsx`
- `data/paad_maaf.xlsx`
- `data/gene_sl_gene.csv`

Python dependencies are NumPy, pandas and SciPy. Use `--help` to override the
repository root, output directory or GENCODE mapping path.

## Outputs

- `outputs/analysis_summary.json`: principal dimensions and results.
- `outputs/sample_burdens.csv`: case-level candidate and all-gene LOF burden,
  mutation count and expression library size.
- `outputs/candidate_only_lof_correlations.csv`: complete primary screen.
- `outputs/leave_top_out_correlations.csv`: complete sensitivity screen after
  excluding TCGA-IB-7651.

The GENCODE mapping is retained in this directory so that the analysis does not
depend on a changing remote annotation endpoint.
