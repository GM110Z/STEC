#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def clean_id(x):
    return str(x).strip().replace("\ufeff", "")


def normalise_columns(df):
    df.columns = (
        df.columns.astype(str)
        .str.replace("\ufeff", "", regex=False)
        .str.strip()
    )
    return df


def read_member_file(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df = normalise_columns(df)

    possible_id_cols = [
        "seqid_clean",
        "seqid",
        "Seqid",
        "locus_id",
        "locus id",
        "matrix_id"
    ]

    id_col = None
    for c in possible_id_cols:
        if c in df.columns:
            id_col = c
            break

    if id_col is None:
        raise ValueError(
            f"No usable seqid column found in {path}.\n"
            f"Available columns: {list(df.columns)}"
        )

    cluster_col = None
    for c in ["cluster_clean", "cluster"]:
        if c in df.columns:
            cluster_col = c
            break

    df["matrix_id"] = df[id_col].apply(clean_id)

    if cluster_col:
        df["cluster_for_plot"] = df[cluster_col].apply(clean_id)
    else:
        df["cluster_for_plot"] = "NA"

    # Only representative file has this column
    if "representative_type" in df.columns:
        df["stx_status_for_plot"] = df["representative_type"].replace({
            "stx_positive": "Stx+",
            "stx_negative": "Stx-"
        })
    elif "is_stx_phage" in df.columns:
        df["stx_status_for_plot"] = df["is_stx_phage"].astype(str).replace({
            "True": "Stx+",
            "False": "Stx-",
            "true": "Stx+",
            "false": "Stx-",
            "1": "Stx+",
            "0": "Stx-"
        })
    else:
        df["stx_status_for_plot"] = ""

    return df, id_col, cluster_col


def read_similarity_matrix(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = normalise_columns(df)

    first_col = df.columns[0]
    first_col_numeric = pd.to_numeric(df[first_col], errors="coerce").notna().all()

    if not first_col_numeric:
        df = df.set_index(first_col)

    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()

    df = df.apply(pd.to_numeric, errors="coerce")

    return df


def cluster_tick_positions(order_df, representative_mode=False):
    tmp = order_df.reset_index(drop=True).copy()
    tmp["pos"] = tmp.index

    tick_positions = []
    tick_labels = []
    boundaries = []

    if representative_mode:
        for _, row in tmp.iterrows():
            cluster = row["cluster_for_plot"]
            status = row["stx_status_for_plot"]

            label = f"{cluster} {status}".strip()

            tick_positions.append(row["pos"])
            tick_labels.append(label)
            boundaries.append(row["pos"] + 0.5)

        if boundaries:
            boundaries = boundaries[:-1]

    else:
        for cluster, sub in tmp.groupby("cluster_for_plot", sort=False):
            start = sub["pos"].min()
            end = sub["pos"].max()
            midpoint = (start + end) / 2

            tick_positions.append(midpoint)
            tick_labels.append(str(cluster))
            boundaries.append(end + 0.5)

        if boundaries:
            boundaries = boundaries[:-1]

    return tick_positions, tick_labels, boundaries


def make_heatmap(
    sub,
    order_df,
    out_prefix,
    title,
    cmap,
    vmin=None,
    vmax=None,
    representative_mode=False
):
    n = sub.shape[0]
    fig_size = max(6, min(30, n * 0.28))

    plt.figure(figsize=(fig_size, fig_size))
    plt.imshow(sub.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(label="Proteome similarity")

    tick_positions, tick_labels, boundaries = cluster_tick_positions(
        order_df,
        representative_mode=representative_mode
    )

    plt.xticks(tick_positions, tick_labels, rotation=90, fontsize=7)
    plt.yticks(tick_positions, tick_labels, fontsize=7)

    for b in boundaries:
        plt.axhline(b, linewidth=0.5)
        plt.axvline(b, linewidth=0.5)

    if representative_mode:
        plt.xlabel("Cluster / Stx status")
        plt.ylabel("Cluster / Stx status")
    else:
        plt.xlabel("Cluster")
        plt.ylabel("Cluster")

    plt.title(title)
    plt.tight_layout()

    plt.savefig(f"{out_prefix}.heatmap.png", dpi=300)
    plt.savefig(f"{out_prefix}.heatmap.pdf")
    plt.close()


def make_cluster_mean_heatmap(sub, annotation, out_prefix, cmap, vmin=None, vmax=None):
    annot = (
        annotation
        .drop_duplicates("matrix_id")
        .set_index("matrix_id")["cluster_for_plot"]
    )

    sub2 = sub.copy()
    sub2.index.name = "seqid_1"
    sub2.columns.name = "seqid_2"

    long = sub2.stack(dropna=False).reset_index(name="similarity")

    long["cluster_1"] = long["seqid_1"].map(annot)
    long["cluster_2"] = long["seqid_2"].map(annot)

    long = long.dropna(subset=["cluster_1", "cluster_2", "similarity"])

    if long.empty:
        print(f"Skipping cluster-mean heatmap for {out_prefix}: no valid values")
        return

    cluster_mean = (
        long.groupby(["cluster_1", "cluster_2"], as_index=False)["similarity"]
        .mean()
        .pivot(index="cluster_1", columns="cluster_2", values="similarity")
    )

    try:
        row_order = sorted(cluster_mean.index, key=lambda x: int(float(x)))
        col_order = sorted(cluster_mean.columns, key=lambda x: int(float(x)))
    except Exception:
        row_order = sorted(cluster_mean.index)
        col_order = sorted(cluster_mean.columns)

    cluster_mean = cluster_mean.loc[row_order, col_order]

    cluster_mean.to_csv(f"{out_prefix}.cluster_mean_similarity.tsv", sep="\t")

    n_rows, n_cols = cluster_mean.shape
    fig_w = max(5, min(18, n_cols * 0.55))
    fig_h = max(5, min(18, n_rows * 0.55))

    plt.figure(figsize=(fig_w, fig_h))
    plt.imshow(cluster_mean.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(label="Mean proteome similarity")

    plt.xticks(range(n_cols), cluster_mean.columns, rotation=90, fontsize=8)
    plt.yticks(range(n_rows), cluster_mean.index, fontsize=8)

    plt.xlabel("Cluster")
    plt.ylabel("Cluster")
    plt.title("Mean similarity between clusters")
    plt.tight_layout()

    plt.savefig(f"{out_prefix}.cluster_mean_heatmap.png", dpi=300)
    plt.savefig(f"{out_prefix}.cluster_mean_heatmap.pdf")
    plt.close()


def process_set(member_file, matrix, label, out_prefix, cmap, vmin, vmax, representative_mode=False):
    members, id_col, cluster_col = read_member_file(member_file)

    wanted = list(dict.fromkeys(members["matrix_id"]))
    matrix_ids = set(matrix.index).intersection(set(matrix.columns))

    found = [x for x in wanted if x in matrix_ids]
    missing = [x for x in wanted if x not in matrix_ids]

    if len(found) == 0:
        raise ValueError(
            f"No IDs from {member_file} were found in the similarity matrix.\n"
            f"Check ID format."
        )

    order_cols = ["matrix_id", "cluster_for_plot", "stx_status_for_plot"]

    order_df = (
        members[members["matrix_id"].isin(found)]
        [order_cols]
        .drop_duplicates()
    )

    try:
        order_df["_cluster_sort"] = order_df["cluster_for_plot"].astype(float)
    except Exception:
        order_df["_cluster_sort"] = order_df["cluster_for_plot"]

    # For reps, force Stx+ before Stx- within each cluster
    status_order = {"Stx+": 0, "Stx-": 1, "": 2}
    order_df["_status_sort"] = order_df["stx_status_for_plot"].map(status_order).fillna(2)

    order_df = order_df.sort_values(
        ["_cluster_sort", "_status_sort", "matrix_id"]
    ).drop(columns=["_cluster_sort", "_status_sort"])

    ordered_ids = order_df["matrix_id"].tolist()
    sub = matrix.loc[ordered_ids, ordered_ids]

    prefix = f"{out_prefix}.{label}"

    sub.to_csv(f"{prefix}.subset_matrix.tsv", sep="\t")
    order_df.to_csv(f"{prefix}.seqid_cluster_annotation.tsv", sep="\t", index=False)

    with open(f"{prefix}.missing_from_matrix.txt", "w") as handle:
        for x in missing:
            handle.write(x + "\n")

    make_heatmap(
        sub=sub,
        order_df=order_df,
        out_prefix=prefix,
        title=f"{label}: proteome similarity",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        representative_mode=representative_mode
    )

    # Keep cluster-mean only for full member set.
    # For representative heatmap, cluster/Stx labels are more useful.
    if not representative_mode:
        make_cluster_mean_heatmap(
            sub=sub,
            annotation=order_df,
            out_prefix=prefix,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax
        )

    print()
    print(label)
    print(f"  Input file: {member_file}")
    print(f"  ID column used: {id_col}")
    print(f"  Cluster column used: {cluster_col}")
    print(f"  IDs requested: {len(wanted)}")
    print(f"  IDs found in matrix: {len(found)}")
    print(f"  IDs missing from matrix: {len(missing)}")
    print(f"  Matrix extracted: {sub.shape[0]} x {sub.shape[1]}")
    print(f"  Wrote: {prefix}.subset_matrix.tsv")
    print(f"  Wrote: {prefix}.heatmap.png/pdf")


def main():
    parser = argparse.ArgumentParser(
        description="Plot similarity heatmaps for all Stx cluster members and Stx+/Stx- representatives."
    )

    parser.add_argument(
        "all_members",
        help="stx_coord_clusters.all_members_of_filtered_stx_clusters.tsv"
    )

    parser.add_argument(
        "representatives",
        help="stx_coord_clusters.one_stx_one_nonstx_per_filtered_cluster.tsv"
    )

    parser.add_argument(
        "matrix",
        help="proteome_similarity_matrix.tsv"
    )

    parser.add_argument(
        "--out_prefix",
        default="stx_similarity",
        help="Output prefix"
    )

    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Colour map: viridis, cividis, magma, plasma, inferno, YlGnBu, Blues, turbo"
    )

    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Minimum colour scale value"
    )

    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Maximum colour scale value"
    )

    args = parser.parse_args()

    matrix = read_similarity_matrix(args.matrix)

    print(f"Similarity matrix loaded: {matrix.shape[0]} x {matrix.shape[1]}")

    process_set(
        member_file=args.all_members,
        matrix=matrix,
        label="all_members_filtered_clusters",
        out_prefix=args.out_prefix,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        representative_mode=False
    )

    process_set(
        member_file=args.representatives,
        matrix=matrix,
        label="one_stx_one_nonstx_representatives",
        out_prefix=args.out_prefix,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        representative_mode=True
    )

    print()
    print("Done.")


if __name__ == "__main__":
    main()
