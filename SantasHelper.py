#!/usr/bin/env python3
import sys
import argparse
import pandas as pd


def main():
    p = argparse.ArgumentParser(
        description="Merge operon table with PFAM domains, PADLOC systems, and DefenseFinder (antidefense) annotations."
    )
    p.add_argument("--operon", required=True, help="Operon results table (CSV). Must contain column 'Protein_ID'.")
    p.add_argument("--pfam", required=True, help="PFAM domain table (TSV). Must contain 'target name', 'query name', 'E-value', 'description of target', 'accession'.")
    p.add_argument("--padloc", required=True, help="PADLOC results (TSV). Must contain 'system' and 'target.name'.")
    p.add_argument("--defensefinder", required=True, help="DefenseFinder/antidefense results (TSV). Must contain 'type', 'subtype', and 'sys_beg'.")
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    args = p.parse_args()

    # --- Load inputs ---
    # Operon table with Protein and operon info
    df_operon = pd.read_csv(args.operon, sep=",")
    required_operon = {"Protein_ID", "nuccore_id", "start", "stop", "strand", "operon_number", "product", "accession", "ClusterRep"}
    missing_operon = required_operon - set(df_operon.columns)
    if missing_operon:
        raise SystemExit(f"[operon] Missing required columns: {sorted(missing_operon)}")

    # PFAM table with domain info
    df_pfam = pd.read_csv(args.pfam, sep="\t")
    required_pfam = {"target name", "query name", "E-value", "description of target", "accession"}
    missing_pfam = required_pfam - set(df_pfam.columns)
    if missing_pfam:
        raise SystemExit(f"[pfam] Missing required columns: {sorted(missing_pfam)}")

    # PADLOC results
    df_padloc = pd.read_csv(args.padloc, sep="\t")
    required_padloc = {"system", "target.name"}
    missing_padloc = required_padloc - set(df_padloc.columns)
    if missing_padloc:
        raise SystemExit(f"[padloc] Missing required columns: {sorted(missing_padloc)}")

    # DefenseFinder / antidefense table
    df_def = pd.read_csv(args.defensefinder, sep="\t")
    required_def = {"type", "subtype", "sys_beg"}
    missing_def = required_def - set(df_def.columns)
    if missing_def:
        raise SystemExit(f"[defensefinder] Missing required columns: {sorted(missing_def)}")

    # --- Merge PFAM onto operon by Protein_ID vs 'target name' ---
    merged = pd.merge(
        df_operon,
        df_pfam,
        left_on="Protein_ID",
        right_on="target name",
        how="left"
    )

    # Select columns for final base table
    base_cols = [
        "Protein_ID", "nuccore_id", "start", "stop", "strand",
        "operon_number", "product", "accession", "ClusterRep",
        "query name", "E-value", "description of target"
    ]
    # Some rows may not have PFAM hits; those PFAM columns will be NaN, that's fine.
    final_df = merged[base_cols].copy()

    # --- Merge PADLOC (system by Protein_ID == target.name) ---
    final_df = pd.merge(
        final_df,
        df_padloc[["system", "target.name"]],
        left_on="Protein_ID",
        right_on="target.name",
        how="left"
    ).drop(columns=["target.name"])

    # --- Prepare DefenseFinder 'Antidefense' label and merge by sys_beg == Protein_ID ---
    df_def["Antidefense"] = df_def["type"].astype(str) + "_" + df_def["subtype"].astype(str)
    final_df = pd.merge(
        final_df,
        df_def[["Antidefense", "sys_beg"]],
        left_on="Protein_ID",
        right_on="sys_beg",
        how="left"
    ).drop(columns=["sys_beg"])

    # --- Write output ---
    final_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
