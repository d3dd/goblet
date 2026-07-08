# GOBLET — Goner-Buffering Lethality Transcriptomic Pipeline

A computational pipeline for identifying transcriptional buffering candidates
in cancer, based on correlating Goner loss-of-function (LOF) burden with
tumour gene expression data.

## Overview

GOBLET identifies genes whose expression co-varies with the burden of Goner
LOF mutations across a tumour cohort. Goners are mutations in essential genes
that are toxic to tumour cells — present at lower-than-expected frequency in
tumour sequencing data because cells acquiring them are eliminated. When
tumours do carry such mutations, compensatory transcriptional changes
(buffering) may allow survival. GOBLET identifies candidate buffering genes
by correlating Goner LOF burden with expression across TCGA tumour samples,
then filtering for known synthetic lethal (SL) pairs.

The proof-of-concept application is TCGA-PAAD (pancreatic adenocarcinoma,
143 samples), where the top hit is BCAT2 (Spearman ρ = 0.31), a branched-chain
amino acid transaminase with independent experimental support in PDAC
(Ericksen et al. 2019; Zhu et al. 2020).

## Repository structure

\`\`\`
goblet/
├── GonerFull24Feb2026.ipynb   # Main analysis notebook
├── goner_list.txt             # Common-essential gene list (DepMap 23Q2,
│                              # Chronos < −0.5 in ≥ 90% of cell lines)
├── requirements.txt           # Python dependencies
├── data/
│   ├── paad_expr_matrix.xlsx  # TCGA-PAAD expression matrix (143 samples)
│   ├── paad_maaf.xlsx         # TCGA-PAAD mutation annotation flat file
│   ├── Human_Mouse_Common.csv # Human–mouse orthologue pairs
│   ├── HSIAO_Cleaned_paired_genes.csv  # Co-essential gene pairs (Hsiao et al.)
│   └── gene_sl_gene.csv       # Synthetic lethal pairs (SynLethDB)
└── README.md
\`\`\`

## Reproducing the results

### Option A — run in Google Colab (easiest)

The notebook was developed in Google Colab. To run it there:

1. Open [Google Colab](https://colab.research.google.com/)
2. File → Open notebook → GitHub → \`d3dd/goblet\` → \`GonerFull24Feb2026.ipynb\`
3. When prompted by \`google.colab.files.upload()\` cells, upload the
   corresponding files from the \`data/\` directory in this repository
4. Run all cells in order

### Option B — run locally in Jupyter

The notebook uses \`google.colab.files.upload()\` for file input, which does
not work outside Colab. To run locally, replace each upload cell with:

\`\`\`python
# Replace: uploaded = cf.upload()
# With:
import pandas as pd
df = pd.read_excel('data/paad_expr_matrix.xlsx', index_col=0)
\`\`\`

Then install dependencies and run:

\`\`\`bash
git clone https://github.com/d3dd/goblet
cd goblet
pip install -r requirements.txt
jupyter notebook GonerFull24Feb2026.ipynb
\`\`\`

## Downloading fresh data from primary sources

If you wish to reproduce the pipeline with updated data, follow these steps:

### TCGA-PAAD expression matrix (\`paad_expr_matrix.xlsx\`)

1. Go to [cBioPortal](https://www.cbioportal.org/)
2. Select study: **Pancreatic Adenocarcinoma (TCGA, PanCancer Atlas)**
   (\`paad_tcga_pan_can_atlas_2018\`)
3. Download → mRNA expression (RNA-seq V2 RSEM)
4. Filter samples: retain only those with tumour purity ≥ 0.70
   (ABSOLUTE estimates available from the TCGA pan-cancer atlas;
   see Ahn et al. 2021 for rationale)
5. This yields **143 samples** matching the cohort used in the manuscript
6. Save as \`paad_expr_matrix.xlsx\` with genes as rows, sample IDs as columns

### TCGA-PAAD mutation annotation flat file (\`paad_maaf.xlsx\`)

1. From the same cBioPortal study, download → Mutations (MAF format)
2. Filter to non-synonymous somatic mutations only
3. Save as \`paad_maaf.xlsx\`

### Common-essential gene list (\`goner_list.txt\`)

Already deposited. Derived from DepMap 23Q2 Chronos scores: genes with
mean Chronos score < −0.5 in ≥ 90% of 1,019 cell lines. To regenerate:

\`\`\`python
import pandas as pd
chronos = pd.read_csv('CRISPRGeneEffect.csv', index_col=0)
essential = (chronos < -0.5).mean() >= 0.90
goner_list = chronos.columns[essential].tolist()
\`\`\`

Download \`CRISPRGeneEffect.csv\` from
[DepMap 23Q2](https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2).

### Synthetic lethal pairs (\`gene_sl_gene.csv\`)

Downloaded from [SynLethDB 2.0](https://synlethdb.sist.shanghaitech.edu.cn/).
Select: human, experimentally validated pairs only.

### Co-essential pairs (\`HSIAO_Cleaned_paired_genes.csv\`)

Derived from Hsiao et al. (2019). See notebook comments for filtering criteria.

### Human–mouse orthologues (\`Human_Mouse_Common.csv\`)

Downloaded from Ensembl BioMart. Homo sapiens → Mus musculus one-to-one
orthologues, protein-coding genes only.

## Dependencies

\`\`\`
pandas>=1.5.0
numpy>=1.23.0
scipy>=1.9.0
openpyxl>=3.0.10
matplotlib>=3.6.0
seaborn>=0.12.0
statsmodels>=0.13.0
gseapy>=1.0.4
biorosetta>=0.3.0
requests>=2.28.0
jupyter>=1.0.0
\`\`\`

Install with:

\`\`\`bash
pip install -r requirements.txt
\`\`\`

Note: \`google.colab\` is pre-installed in Colab and should not be pip-installed.

## Tumour purity filtering

Following Ahn, Grimes & Datta (2021, *Frontiers in Genetics*), samples with
tumour purity below 0.70 are excluded before expression correlation to avoid
confounding by stromal/immune cell transcription. ABSOLUTE purity estimates
for TCGA-PAAD are available from the TCGA pan-cancer atlas supplementary data.

## PAAD-specific essentiality threshold

The default \`goner_list.txt\` uses a pan-cancer threshold (Chronos < −0.5 in
≥ 90% of all 1,019 DepMap cell lines). A PAAD-specific sensitivity analysis
using only PAAD-lineage cell lines is planned. The PAAD cell lines can be
identified from the DepMap 23Q2 model metadata
(\`Model.csv\`, \`OncotreeLineage == 'Pancreas'\`):

\`\`\`python
import pandas as pd

# Load DepMap 23Q2 model metadata
models = pd.read_csv('Model.csv', index_col=0)
paad_lines = models[models['OncotreeLineage'] == 'Pancreas'].index

# Load Chronos scores and subset to PAAD lines
chronos = pd.read_csv('CRISPRGeneEffect.csv', index_col=0)
chronos_paad = chronos.loc[chronos.index.isin(paad_lines)]

# Apply threshold: essential in ≥ 80% of PAAD lines
# (threshold relaxed from 90% due to smaller n)
essential_paad = (chronos_paad < -0.5).mean() >= 0.80
goner_list_paad = chronos_paad.columns[essential_paad].tolist()
\`\`\`

## Limitations

**Pan-cancer essentiality threshold:** The current Goner list uses a pan-cancer
DepMap threshold. This may include genes that are not essential in PAAD
specifically, or miss PAAD-specific dependencies. A sensitivity analysis
comparing pan-cancer vs PAAD-specific thresholds is in progress.

**Stromal contamination:** TCGA bulk RNA-seq includes stromal and immune cell
transcription. Tumour purity filtering (≥ 0.70, see above) partially addresses
this. Note that for the BCAT2 finding specifically, Zhu et al. (2020) show
that stromal BCAT2 (in cancer-associated fibroblasts) has no effect on PDAC
growth — only tumour-cell BCAT2 matters — which supports the interpretation
of the GOBLET signal as tumour-cell-intrinsic.

**Causal interpretation:** The pipeline identifies correlation between Goner
LOF burden and gene expression — consistent with transcriptional buffering
but not proof of causation. Candidate genes require independent experimental
validation (see Ericksen 2019 for BCAT2).

## Key results

| Gene  | Pathway              | Spearman ρ | FDR-adjusted p | Rank |
|-------|----------------------|-----------|----------------|------|
| BCAT2 | BCAA catabolism      | 0.31      | < 0.05         | 1 (BCAA) |
| DBT   | BCAA catabolism (E2) | 0.30      | < 0.05         | 2 (BCAA) |

The BCAA catabolism pathway (GO:0009083) is the strongest enriched GO term
among the top 25 candidates (adjusted p = 0.0438, 4 enriched terms).
BCAT2 upregulation in PDAC is independently supported by:
- Ericksen et al. (2019): BCAT2 knockdown selectively impairs PDAC but not
  normal pancreatic cell proliferation
- Zhu et al. (2020): tumour-cell BCAT2 (not stromal BCAT2) drives BCAA
  dependency in PDAC

## Citation

If you use GOBLET, please cite:

> Boman M et al. (2026). Transcriptomic Buffering Detection: A Computational
> Approach to Identifying Essential Gene Dependencies in Cancer.
> *Bioinformatics Advances*, BIOADV-2026-298 (under revision).

## Licence

MIT
