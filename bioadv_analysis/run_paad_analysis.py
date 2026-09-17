#!/usr/bin/env python3
"""Run the GOBLET 1.0 TCGA-PAAD analysis reported in BIOADV-2026-434.

The script assigns the deposited expression rows using the sorted GENCODE v36
gene set, forms a candidate-only LOF burden, tests transcriptome-wide Spearman
associations and repeats the screen after excluding the highest-burden case.
"""

from __future__ import annotations

import argparse
import json
import re
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree.ElementTree import iterparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist


SCRIPT_DIR = Path(__file__).resolve().parent

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
LOF_TYPES = {
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonsense_Mutation",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
}


def cell_value(cell):
    value = cell.find(NS + "v")
    if value is not None:
        return value.text or ""
    inline = cell.find(NS + "is")
    if inline is not None:
        return "".join(node.text or "" for node in inline.iter(NS + "t"))
    return ""


def column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    result = 0
    for letter in letters:
        result = result * 26 + ord(letter) - 64
    return result - 1


def read_expression_xlsx(path: Path) -> tuple[list[str], np.ndarray]:
    """Read the single-sheet numeric workbook without materialising XML strings."""
    sample_ids: list[str] | None = None
    matrix: np.ndarray | None = None
    data_row = 0
    with zipfile.ZipFile(path) as archive, archive.open(
        "xl/worksheets/sheet1.xml"
    ) as handle:
        for _, element in iterparse(handle, events=("end",)):
            if element.tag != NS + "row":
                continue
            cells = element.findall(NS + "c")
            row_number = int(element.attrib["r"])
            if row_number == 1:
                sample_ids = [cell_value(cell) for cell in cells]
                matrix = np.zeros((60660, len(sample_ids)), dtype=np.float32)
            else:
                if matrix is None:
                    raise ValueError("Expression workbook data preceded its header")
                for cell in cells:
                    matrix[data_row, column_index(cell.attrib["r"])] = float(
                        cell_value(cell) or 0
                    )
                data_row += 1
            element.clear()
    if sample_ids is None or matrix is None or data_row != matrix.shape[0]:
        raise ValueError(f"Expression workbook parsed as {data_row} data rows")
    return sample_ids, matrix


def read_maf_xlsx(path: Path) -> pd.DataFrame:
    wanted = {"Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"}
    records: list[dict[str, str]] = []
    header: dict[int, str] = {}
    with zipfile.ZipFile(path) as archive, archive.open(
        "xl/worksheets/sheet1.xml"
    ) as handle:
        for _, element in iterparse(handle, events=("end",)):
            if element.tag != NS + "row":
                continue
            row_number = int(element.attrib["r"])
            values = {
                column_index(cell.attrib["r"]): cell_value(cell)
                for cell in element.findall(NS + "c")
            }
            if row_number == 1:
                header = {index: value for index, value in values.items()}
                missing = wanted - set(header.values())
                if missing:
                    raise ValueError(f"Missing MAF columns: {sorted(missing)}")
            else:
                record = {
                    name: values.get(index, "")
                    for index, name in header.items()
                    if name in wanted
                }
                records.append(record)
            element.clear()
    return pd.DataFrame.from_records(records)


def bh_correlations(
    expression: np.ndarray,
    burden: np.ndarray,
    tested: np.ndarray,
) -> pd.DataFrame:
    """Vectorised Spearman correlations and asymptotic two-sided p-values."""
    x = expression[tested]
    ranked_x = rankdata(x, axis=1, method="average")
    ranked_y = rankdata(burden, method="average")
    ranked_x -= ranked_x.mean(axis=1, keepdims=True)
    ranked_y -= ranked_y.mean()
    denominator = np.sqrt(
        np.sum(ranked_x * ranked_x, axis=1) * np.sum(ranked_y * ranked_y)
    )
    rho = np.divide(
        ranked_x @ ranked_y,
        denominator,
        out=np.full(ranked_x.shape[0], np.nan),
        where=denominator != 0,
    )
    degrees = expression.shape[1] - 2
    clipped = np.clip(rho, -0.999999999, 0.999999999)
    statistic = clipped * np.sqrt(degrees / (1.0 - clipped * clipped))
    p_value = 2 * t_dist.sf(np.abs(statistic), degrees)
    finite_p = np.nan_to_num(p_value, nan=1.0)
    order = np.argsort(finite_p)
    ranked = finite_p[order]
    adjusted_ranked = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    p_adjusted = np.empty_like(adjusted_ranked)
    p_adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return pd.DataFrame(
        {
            "row_index": np.flatnonzero(tested),
            "rho": rho,
            "p_value": p_value,
            "p_adjusted": p_adjusted,
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        type=Path,
        default=SCRIPT_DIR.parent,
        help="GOBLET repository root (default: parent of this directory)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIR / "outputs",
        help="Output directory (default: bioadv_analysis/outputs)",
    )
    parser.add_argument(
        "--gencode",
        type=Path,
        default=SCRIPT_DIR / "gencode.v36.genes.sorted.tsv",
        help="Sorted GENCODE v36 gene mapping",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo = args.repo.resolve()
    output = args.output.resolve()
    expression_xlsx = repo / "data" / "paad_expr_matrix.xlsx"
    maf_xlsx = repo / "data" / "paad_maaf.xlsx"
    goners_path = repo / "goner_list.txt"
    sl_path = repo / "data" / "gene_sl_gene.csv"
    output.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(
        args.gencode,
        sep="\t",
        names=["ensembl_gene_id", "gene_symbol", "gene_type"],
        dtype=str,
    )
    if len(genes) != 60660 or genes["ensembl_gene_id"].iloc[0] != "ENSG00000000003.15":
        raise ValueError("GENCODE v36 mapping did not pass row-order checks")

    samples_full, raw_counts = read_expression_xlsx(expression_xlsx)
    maf = read_maf_xlsx(maf_xlsx)
    maf["sample_id"] = maf["Tumor_Sample_Barcode"].str[:12]
    goners = {
        line.strip().upper()
        for line in goners_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }

    # Keep the first aliquot per case, mirroring the deposited notebook.
    sample_ids = [sample[:12] for sample in samples_full]
    first_columns: dict[str, int] = {}
    for column, sample in enumerate(sample_ids):
        first_columns.setdefault(sample, column)

    lof = maf[maf["Variant_Classification"].isin(LOF_TYPES)].copy()
    lof_unique = lof[["sample_id", "Hugo_Symbol"]].drop_duplicates()
    lof_unique["Hugo_Symbol"] = lof_unique["Hugo_Symbol"].str.upper()
    goner_lof = lof_unique[lof_unique["Hugo_Symbol"].isin(goners)]

    maf_samples = set(maf["sample_id"])
    common = sorted(maf_samples & set(first_columns))
    columns = np.array([first_columns[sample] for sample in common])
    raw_counts = raw_counts[:, columns]

    all_lof_counts = (
        lof_unique.groupby("sample_id")["Hugo_Symbol"].nunique().reindex(common, fill_value=0)
    )
    goner_counts = (
        goner_lof.groupby("sample_id")["Hugo_Symbol"].nunique().reindex(common, fill_value=0)
    )
    total_mutations = maf.groupby("sample_id").size().reindex(common, fill_value=0)

    # Library-size normalisation.  log2 is monotone within a gene, but CPM
    # removes the between-sample sequencing-depth component before ranking.
    library_sizes = raw_counts.sum(axis=0)
    cpm = raw_counts / library_sizes[np.newaxis, :] * 1_000_000.0
    expression = np.log2(cpm + 1.0)
    tested = (cpm >= 1.0).sum(axis=1) >= max(1, int(np.ceil(0.20 * len(common))))

    full = bh_correlations(expression, goner_counts.to_numpy(dtype=float), tested)
    full = full.merge(genes.reset_index(names="row_index"), on="row_index", how="left")
    full = full.sort_values(["p_adjusted", "p_value", "rho"], ascending=[True, True, False])

    top_goner_sample = goner_counts.idxmax()
    retained = np.array([sample != top_goner_sample for sample in common])
    sensitivity = bh_correlations(
        expression[:, retained],
        goner_counts.to_numpy(dtype=float)[retained],
        tested,
    )
    sensitivity = sensitivity.merge(
        genes.reset_index(names="row_index"), on="row_index", how="left"
    )
    sensitivity = sensitivity.sort_values(
        ["p_adjusted", "p_value", "rho"], ascending=[True, True, False]
    )

    for name, frame in [
        ("candidate_only_lof_correlations.csv", full),
        ("leave_top_out_correlations.csv", sensitivity),
    ]:
        frame.to_csv(output / name, index=False)

    burden_table = pd.DataFrame(
        {
            "sample_id": common,
            "all_lof_gene_count": all_lof_counts.to_numpy(),
            "goner_lof_gene_count": goner_counts.to_numpy(),
            "total_maf_variant_count": total_mutations.to_numpy(),
            "library_size": library_sizes,
        }
    )
    burden_table.to_csv(output / "sample_burdens.csv", index=False)

    observed_goners = set(goner_lof["Hugo_Symbol"])
    sl = pd.read_csv(sl_path, low_memory=False)
    sl["x_name"] = sl["x_name"].astype(str).str.upper()
    sl["y_name"] = sl["y_name"].astype(str).str.upper()
    partners_of_observed = set(sl.loc[sl["x_name"].isin(observed_goners), "y_name"]) | set(
        sl.loc[sl["y_name"].isin(observed_goners), "x_name"]
    )

    def significant(frame: pd.DataFrame) -> pd.DataFrame:
        return frame[(frame["rho"] > 0.3) & (frame["p_adjusted"] < 0.05)]

    full_sig = significant(full)
    sensitivity_sig = significant(sensitivity)
    full_sl = full_sig[full_sig["gene_symbol"].str.upper().isin(partners_of_observed)]
    sensitivity_sl = sensitivity_sig[
        sensitivity_sig["gene_symbol"].str.upper().isin(partners_of_observed)
    ]
    largest_positive = full.loc[full["rho"].idxmax()]

    summary = {
        "expression_rows": int(raw_counts.shape[0]),
        "expression_aliquots": len(samples_full),
        "unique_expression_cases": len(first_columns),
        "maf_rows": len(maf),
        "lof_rows": len(lof),
        "unique_lof_genes": int(lof["Hugo_Symbol"].nunique()),
        "goner_list_size": len(goners),
        "observed_candidate_genes": len(observed_goners),
        "candidate_gene_case_lof_observations": int(goner_counts.sum()),
        "aligned_cases": len(common),
        "tested_expression_features": int(tested.sum()),
        "goner_burden_distribution": {
            str(key): int(value)
            for key, value in sorted(Counter(goner_counts.to_numpy()).items())
        },
        "top_goner_sample": top_goner_sample,
        "top_goner_burden": int(goner_counts.loc[top_goner_sample]),
        "top_sample_total_maf_variants": int(total_mutations.loc[top_goner_sample]),
        "median_total_maf_variants": float(total_mutations.median()),
        "significant_features": len(full_sig),
        "sl_overlap_features": len(full_sl),
        "leave_top_out_significant_features": len(sensitivity_sig),
        "leave_top_out_sl_overlap_features": len(sensitivity_sl),
        "largest_positive_association": {
            "ensembl_gene_id": largest_positive["ensembl_gene_id"],
            "gene_symbol": largest_positive["gene_symbol"],
            "rho": float(largest_positive["rho"]),
            "p_value": float(largest_positive["p_value"]),
            "p_adjusted": float(largest_positive["p_adjusted"]),
        },
        "significant_genes": full_sig[
            ["ensembl_gene_id", "gene_symbol", "rho", "p_value", "p_adjusted"]
        ].to_dict("records"),
        "leave_top_out_significant_genes": sensitivity_sig[
            ["ensembl_gene_id", "gene_symbol", "rho", "p_value", "p_adjusted"]
        ].to_dict("records"),
    }
    (output / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
