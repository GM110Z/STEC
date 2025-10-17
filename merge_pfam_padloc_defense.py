#!/usr/bin/env python3
import sys
import argparse
import pandas as pd

def read_table(path, sep="\t"):
    return pd.read_csv(path, sep=sep)

def main():
    p = argparse.ArgumentParser(
        description="Merge PFAM, PADLOC, DefenseFinder, and optional AMRFinder by Protein_ID."
    )
    p.add_argument("--pfam", required=True, help="PFAM domtbl/parsed table (TSV). Must contain 'target name' or 'target.name'")
    p.add_argument("--padloc", required=True, help="PADLOC output (TSV). Must contain 'target.name' and 'system'")
    p.add_argument("--antidefense", required=True, help="DefenseFinder/anti-defense table (TSV). Must contain 'type','subtype','sys_beg'")
    p.add_argument("--amrfinder", required=False, help="Optional AMRFinder table (TSV).")
    p.add_argument("-o", "--out", default="Final-R-file.tsv", help="Output TSV (default: Final-R-file.tsv)")
    args = p.parse_args()

    # ---------- PFAM ----------
    df_pfam = read_table(args.pfam, sep="\t")
    # Normalize target id column name -> Protein_ID
    if "target name" in df_pfam.columns:
        df_pfam = df_pfam.rename(columns={"target name": "Protein_ID"})
    elif "target.name" in df_pfam.columns:
        df_pfam = df_pfam.rename(columns={"target.name": "Protein_ID"})
    else:
        sys.exit("PFAM file must have a 'target name' or 'target.name' column.")

    # Keep a common, lean set of PFAM columns if they exist
    pfam_cols = ["Protein_ID"]
    for c in ["query name", "E-value", "description of target"]:
        if c in df_pfam.columns:
            pfam_cols.append(c)
    df_pfam = df_pfam[pfam_cols].drop_duplicates(subset=["Protein_ID"])

    # ---------- PADLOC ----------
    df_padloc = read_table(args.padloc, sep="\t")
    if "target.name" not in df_padloc.columns:
        sys.exit("PADLOC file must have 'target.name' column.")
    if "system" not in df_padloc.columns:
        sys.exit("PADLOC file must have 'system' column.")
    df_padloc = df_padloc.rename(columns={"target.name": "Protein_ID"})
    df_padloc = df_padloc[["Protein_ID", "system"]].drop_duplicates(subset=["Protein_ID"])

    # ---------- DefenseFinder / Anti-defense ----------
    df_anti = read_table(args.antidefense, sep="\t")
    for col in ["type", "subtype", "sys_beg"]:
        if col not in df_anti.columns:
            sys.exit("Anti-defense file must have 'type', 'subtype', and 'sys_beg' columns.")
    df_anti["Antidefense"] = df_anti["type"].astype(str) + "_" + df_anti["subtype"].astype(str)
    # Note: anti-defense maps protein via 'sys_beg' -> Protein_ID
    df_anti = df_anti[["sys_beg", "Antidefense"]].drop_duplicates(subset=["sys_beg"])

    # ---------- AMRFinder (optional) ----------
    if args.amrfinder:
        df_amr = read_table(args.amrfinder, sep="\t")
        if "Protein identifier" in df_amr.columns:
            # Build a compact descriptor if columns exist; otherwise keep as-is.
            parts = [
                ("Sequence name", "Sequence name"),
                ("Scope", "Scope"),
                ("Element type", "Element type"),
                ("Element subtype", "Element subtype"),
                ("Class", "Class"),
                ("Subclass", "Subclass"),
            ]
            if all(k in df_amr.columns for k, _ in parts):
                df_amr["AMRFinder"] = (
                    df_amr["Sequence name"].astype(str) + "_" +
                    df_amr["Scope"].astype(str) + "_" +
                    df_amr["Element type"].astype(str) + "_" +
                    df_amr["Element subtype"].astype(str) + "_" +
                    df_amr["Class"].astype(str) + "_" +
                    df_amr["Subclass"].astype(str)
                )
                df_amr = df_amr[["Protein identifier", "AMRFinder"]]
            else:
                # If not all columns present, just keep the identifier
                df_amr["AMRFinder"] = True
                df_amr = df_amr[["Protein identifier", "AMRFinder"]]
        else:
            sys.exit("AMRFinder file must have 'Protein identifier' column if provided.")
    else:
        df_amr = None

    # ---------- Build a master Protein_ID index (union of all sources) ----------
    protein_ids = set(df_pfam["Protein_ID"]) | set(df_padloc["Protein_ID"]) | set(df_anti["sys_beg"])
    if df_amr is not None:
        protein_ids |= set(df_amr["Protein identifier"])
    master = pd.DataFrame({"Protein_ID": sorted(protein_ids)})

    # ---------- Merge everything onto master ----------
    out_df = master.merge(df_pfam, on="Protein_ID", how="left") \
                   .merge(df_padloc, on="Protein_ID", how="left") \
                   .merge(df_anti, left_on="Protein_ID", right_on="sys_beg", how="left") \
                   .drop(columns=["sys_beg"])

    if df_amr is not None:
        out_df = out_df.merge(df_amr, left_on="Protein_ID", right_on="Protein identifier", how="left") \
                       .drop(columns=["Protein identifier"])

    # ---------- Save ----------
    out_df.to_csv(args.out, index=False, sep="\t")
    print(f"Wrote {len(out_df):,} rows to {args.out}")

if __name__ == "__main__":
    main()
