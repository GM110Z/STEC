#!/usr/bin/env python3

import argparse
import os
import re
from collections import Counter

import numpy as np
import pandas as pd


UNCLASSIFIED_VALUES = {
    "",
    "na",
    "nan",
    "none",
    "unknown",
    "unclassified",
    "unclassified caudoviricetes",
}


def clean_id(value):
    """
    Harmonise identifiers such as:

    ADUT01000057.1:133732-152734
    ADUT01000057.1-133732-152734
    ADUT01000057.1-133732-152734.gbff
    """
    value = os.path.basename(str(value).strip())

    value = re.sub(
        r"\.(gbff|gbk|gb|fna|fa|fasta|faa|gff3|gff)$",
        "",
        value,
        flags=re.IGNORECASE,
    )

    # Metadata uses accession:start-end, while filenames may use accession-start-end
    value = re.sub(
        r"^([^:]+):(\d+)-(\d+)$",
        r"\1-\2-\3",
        value,
    )

    return value


def is_classified(value):
    return str(value).strip().lower() not in UNCLASSIFIED_VALUES


def parse_family_genus(value):
    """
    Parse a combined family_genus value such as:

    NewFamily_Marienburgvirus

    The first underscore separates family and genus.
    """
    value = str(value).strip()

    if not is_classified(value):
        return "", ""

    if "_" not in value:
        return value, ""

    return tuple(value.split("_", 1))


def read_similarity_matrix(path):
    matrix = pd.read_csv(
        path,
        sep=None,
        engine="python",
        index_col=0,
    )

    matrix.index = [clean_id(x) for x in matrix.index]
    matrix.columns = [clean_id(x) for x in matrix.columns]

    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    common = [
        genome
        for genome in matrix.index
        if genome in matrix.columns
    ]

    if not common:
        raise SystemExit(
            "The matrix row and column names do not match."
        )

    matrix = matrix.loc[common, common]

    # Convert percentages to proportions when necessary
    finite_values = matrix.to_numpy()
    finite_values = finite_values[np.isfinite(finite_values)]

    if len(finite_values) == 0:
        raise SystemExit("No numeric similarity values found.")

    if finite_values.max() > 1.5:
        print("Similarity values appear to be percentages; dividing by 100.")
        matrix = matrix / 100.0

    # Force symmetry in case of tiny rounding differences
    matrix = (
        matrix
        .combine_first(matrix.T)
        .add(matrix.T, fill_value=0)
        / 2
    )

    np.fill_diagonal(matrix.values, 1.0)

    return matrix


parser = argparse.ArgumentParser(
    description=(
        "Convert a LoVis4u proteomic-similarity matrix and phage metadata "
        "into Cytoscape tables, family-support summaries and a selected "
        "heatmap matrix."
    )
)

parser.add_argument(
    "--matrix",
    required=True,
    help="LoVis4u pairwise similarity matrix",
)

parser.add_argument(
    "--metadata",
    required=True,
    help=(
        "TSV with Genome, Species cluster, Genus cluster and "
        "INPHARED/TaxMyPhage genus columns"
    ),
)

parser.add_argument(
    "--threshold",
    type=float,
    default=0.30,
    help="Minimum similarity for retaining an edge; default 0.30",
)

parser.add_argument(
    "--top-k",
    type=int,
    default=5,
    help=(
        "Retain an edge only if it is among the top-k neighbours of "
        "at least one endpoint; default 5"
    ),
)

parser.add_argument(
    "--heatmap-neighbours",
    type=int,
    default=3,
    help=(
        "Number of classified neighbours retained per unclassified phage "
        "for the reduced heatmap"
    ),
)

parser.add_argument(
    "--output-prefix",
    default="lovis4u",
)

args = parser.parse_args()


# ------------------------------------------------------------
# Read LoVis4u matrix
# ------------------------------------------------------------

matrix = read_similarity_matrix(args.matrix)

print(f"Similarity matrix contains {len(matrix)} phages.")


# ------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------

metadata = pd.read_csv(
    args.metadata,
    sep="\t",
    dtype=str,
).fillna("")

required_columns = {
    "Genome",
    "Species cluster",
    "Genus cluster",
    "INPHARED/TaxMyPhage genus",
}

missing = required_columns - set(metadata.columns)

if missing:
    raise SystemExit(
        "Metadata table is missing: "
        + ", ".join(sorted(missing))
    )

metadata["original_genome"] = metadata["Genome"].str.strip()
metadata["phage"] = metadata["Genome"].map(clean_id)

metadata["species_cluster"] = (
    metadata["Species cluster"]
    .astype(str)
    .str.strip()
)

metadata["genus_cluster"] = (
    metadata["Genus cluster"]
    .astype(str)
    .str.strip()
)

metadata["family_genus"] = (
    metadata["INPHARED/TaxMyPhage genus"]
    .astype(str)
    .str.strip()
)

parsed = metadata["family_genus"].map(parse_family_genus)

metadata["family"] = [x[0] for x in parsed]
metadata["genus"] = [x[1] for x in parsed]

metadata["status"] = np.where(
    metadata["family_genus"].map(is_classified),
    "classified",
    "unclassified",
)

metadata = metadata.drop_duplicates("phage")

metadata_lookup = metadata.set_index("phage").to_dict("index")


# ------------------------------------------------------------
# Match matrix names to metadata
# ------------------------------------------------------------

matrix_phages = set(matrix.index)
metadata_phages = set(metadata["phage"])

matched = matrix_phages & metadata_phages
matrix_only = matrix_phages - metadata_phages
metadata_only = metadata_phages - matrix_phages

print(f"Matched phages: {len(matched)}")
print(f"Matrix phages without metadata: {len(matrix_only)}")
print(f"Metadata phages absent from matrix: {len(metadata_only)}")

pd.DataFrame(
    {"phage": sorted(matrix_only)}
).to_csv(
    f"{args.output_prefix}_matrix_only_names.tsv",
    sep="\t",
    index=False,
)

pd.DataFrame(
    {"phage": sorted(metadata_only)}
).to_csv(
    f"{args.output_prefix}_metadata_only_names.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Cytoscape node table
# ------------------------------------------------------------

node_rows = []

for phage in matrix.index:
    info = metadata_lookup.get(phage, {})

    status = info.get("status", "metadata_missing")
    genus_cluster = info.get("genus_cluster", "")
    family = info.get("family", "")
    genus = info.get("genus", "")
    family_genus = info.get("family_genus", "")

    if status == "unclassified":
        label = (
            f"GC{genus_cluster}"
            if genus_cluster
            else phage
        )
    else:
        label = ""

    node_rows.append({
        "phage": phage,
        "display_label": label,
        "family": family if family else "Unclassified",
        "genus": genus,
        "family_genus": family_genus,
        "species_cluster": info.get("species_cluster", ""),
        "genus_cluster": genus_cluster,
        "status": status,
        "node_shape": (
            "Square"
            if status == "unclassified"
            else "Ellipse"
        ),
    })

nodes = pd.DataFrame(node_rows)

nodes.to_csv(
    f"{args.output_prefix}_cytoscape_nodes.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Rank neighbours for top-k filtering
# ------------------------------------------------------------

top_neighbours = {}

for phage in matrix.index:
    ranked = (
        matrix.loc[phage]
        .drop(index=phage, errors="ignore")
        .dropna()
        .sort_values(ascending=False)
    )

    top_neighbours[phage] = set(
        ranked.head(args.top_k).index
    )


# ------------------------------------------------------------
# Cytoscape edge table
# ------------------------------------------------------------

edge_rows = []
phages = list(matrix.index)

for i, source in enumerate(phages):
    for target in phages[i + 1:]:

        similarity = matrix.at[source, target]

        if pd.isna(similarity):
            continue

        similarity = float(similarity)

        if similarity < args.threshold:
            continue

        top_k_connection = (
            target in top_neighbours[source]
            or source in top_neighbours[target]
        )

        if not top_k_connection:
            continue

        source_info = metadata_lookup.get(source, {})
        target_info = metadata_lookup.get(target, {})

        source_family = source_info.get("family", "")
        target_family = target_info.get("family", "")

        same_family = (
            bool(source_family)
            and source_family == target_family
        )

        edge_rows.append({
            "source": source,
            "target": target,
            "similarity": round(similarity, 4),
            "edge_width": round(1 + similarity * 5, 3),
            "same_known_family": same_family,
            "source_family": source_family,
            "target_family": target_family,
        })

edges = pd.DataFrame(edge_rows)

if not edges.empty:
    edges = edges.sort_values(
        "similarity",
        ascending=False,
    )

edges.to_csv(
    f"{args.output_prefix}_cytoscape_edges.tsv",
    sep="\t",
    index=False,
)

print(f"Retained {len(edges)} network edges.")


# ------------------------------------------------------------
# Closest classified neighbours for every unclassified phage
# ------------------------------------------------------------

classified = set(
    nodes.loc[
        nodes["status"] == "classified",
        "phage",
    ]
)

unclassified = set(
    nodes.loc[
        nodes["status"] == "unclassified",
        "phage",
    ]
)

neighbour_rows = []
support_rows = []
selected_for_heatmap = set(unclassified)

for query in sorted(unclassified):

    if query not in matrix.index:
        continue

    candidates = []

    for reference in classified:

        if reference not in matrix.columns:
            continue

        similarity = matrix.at[query, reference]

        if pd.isna(similarity):
            continue

        reference_info = metadata_lookup.get(reference, {})

        candidates.append({
            "query": query,
            "query_genus_cluster": (
                metadata_lookup
                .get(query, {})
                .get("genus_cluster", "")
            ),
            "reference": reference,
            "reference_family": reference_info.get("family", ""),
            "reference_genus": reference_info.get("genus", ""),
            "reference_family_genus": reference_info.get(
                "family_genus",
                "",
            ),
            "similarity": float(similarity),
        })

    candidates.sort(
        key=lambda row: row["similarity"],
        reverse=True,
    )

    top = candidates[:args.top_k]
    neighbour_rows.extend(top)

    selected_for_heatmap.update(
        row["reference"]
        for row in candidates[:args.heatmap_neighbours]
    )

    if not top:
        support_rows.append({
            "query": query,
            "query_genus_cluster": (
                metadata_lookup
                .get(query, {})
                .get("genus_cluster", "")
            ),
            "best_reference": "",
            "best_similarity": "",
            "dominant_family": "",
            "dominant_family_hits": 0,
            "top_hits_considered": 0,
            "family_consistency": 0,
            "result": "no classified neighbour",
        })
        continue

    family_counts = Counter(
        row["reference_family"]
        for row in top
        if row["reference_family"]
    )

    if family_counts:
        dominant_family, family_hits = family_counts.most_common(1)[0]
        consistency = family_hits / len(top)
    else:
        dominant_family = ""
        family_hits = 0
        consistency = 0

    best = top[0]

    if consistency >= 0.8 and best["similarity"] >= args.threshold:
        result = f"{dominant_family}-associated"
    elif consistency >= 0.6 and best["similarity"] >= args.threshold:
        result = f"possible {dominant_family}-association"
    else:
        result = "ambiguous"

    support_rows.append({
        "query": query,
        "query_genus_cluster": (
            metadata_lookup
            .get(query, {})
            .get("genus_cluster", "")
        ),
        "best_reference": best["reference"],
        "best_reference_family": best["reference_family"],
        "best_reference_genus": best["reference_genus"],
        "best_similarity": round(best["similarity"], 4),
        "dominant_family": dominant_family,
        "dominant_family_hits": family_hits,
        "top_hits_considered": len(top),
        "family_consistency": round(consistency, 4),
        "result": result,
    })


top_neighbours_df = pd.DataFrame(neighbour_rows)

if not top_neighbours_df.empty:
    top_neighbours_df["similarity"] = (
        top_neighbours_df["similarity"].round(4)
    )

top_neighbours_df.to_csv(
    f"{args.output_prefix}_unclassified_top_neighbours.tsv",
    sep="\t",
    index=False,
)

family_support = pd.DataFrame(support_rows)

family_support.to_csv(
    f"{args.output_prefix}_family_support.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Reduced matrix for the clean heatmap
# ------------------------------------------------------------

selected = [
    phage
    for phage in matrix.index
    if phage in selected_for_heatmap
]

reduced_matrix = matrix.loc[selected, selected]

reduced_matrix.to_csv(
    f"{args.output_prefix}_selected_heatmap_matrix.tsv",
    sep="\t",
    index=True,
    index_label="phage",
)

selected_nodes = nodes[
    nodes["phage"].isin(selected)
].copy()

selected_nodes.to_csv(
    f"{args.output_prefix}_selected_heatmap_metadata.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Edge-threshold diagnostic table
# ------------------------------------------------------------

threshold_rows = []

for threshold in np.arange(0.10, 0.81, 0.05):

    count = 0
    involved = set()

    for i, source in enumerate(phages):
        for target in phages[i + 1:]:

            value = matrix.at[source, target]

            if pd.notna(value) and float(value) >= threshold:
                count += 1
                involved.update([source, target])

    threshold_rows.append({
        "threshold": round(float(threshold), 2),
        "edges": count,
        "nodes_with_edges": len(involved),
    })

pd.DataFrame(threshold_rows).to_csv(
    f"{args.output_prefix}_threshold_diagnostics.tsv",
    sep="\t",
    index=False,
)

print("Outputs written:")
print(f"  {args.output_prefix}_cytoscape_nodes.tsv")
print(f"  {args.output_prefix}_cytoscape_edges.tsv")
print(f"  {args.output_prefix}_family_support.tsv")
print(f"  {args.output_prefix}_unclassified_top_neighbours.tsv")
print(f"  {args.output_prefix}_selected_heatmap_matrix.tsv")
print(f"  {args.output_prefix}_selected_heatmap_metadata.tsv")
print(f"  {args.output_prefix}_threshold_diagnostics.tsv")
