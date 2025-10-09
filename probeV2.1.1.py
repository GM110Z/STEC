#!/usr/bin/env python3
import os
import sys
import glob
import subprocess as sp
from Bio import SearchIO, SeqIO

# ---------- HMMER ----------
def run_hmmer(model, fasta_file, domtbl_path, threshold=5, log_path=None):
    if log_path is None:
        log_path = os.devnull
    cmd = (
        f"hmmsearch -T {threshold} --incT {threshold} "
        f"-o {log_path} --domtblout {domtbl_path} {model} {fasta_file}"
    )
    sp.run(cmd, shell=True, check=False)

def read_hit_ids_from_domtbl(domtbl_path):
    """Return a list of hit IDs from hmmsearch --domtblout (domain table)."""
    if not os.path.exists(domtbl_path):
        return []
    hits = []
    with open(domtbl_path, "r") as handle:
        for qresult in SearchIO.parse(handle, "hmmsearch3-domtab"):
            for hit in qresult.hits:
                hits.append(hit.id)
    return hits

# ---------- GenBank helpers ----------
def index_genbank(genbank_file):
    """
    Return:
      - rec_by_id: dict[record.id] = SeqRecord
      - cds_index: dict[protein_or_locus_id] = (start, end, record.id)
    Matches by 'protein_id' and 'locus_tag' (version stripped).
    """
    rec_by_id = {}
    cds_index = {}
    for rec in SeqIO.parse(genbank_file, "genbank"):
        rec_by_id[rec.id] = rec
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            quals = feat.qualifiers
            ids = []
            ids.extend(quals.get("protein_id", []))
            ids.extend(quals.get("locus_tag", []))
            # Strip version suffix
            ids = [i.split(".", 1)[0] for i in ids]
            for pid in ids:
                cds_index[pid] = (int(feat.location.start), int(feat.location.end), rec.id)
    return rec_by_id, cds_index

def map_ids_to_coords(ids, cds_index):
    coords = []
    for pid in ids:
        pid_base = pid.split(".", 1)[0]
        if pid_base in cds_index:
            coords.append((pid_base, *cds_index[pid_base]))  # (pid, start, end, chrom)
    return coords

# ---------- Core pairing / fallback logic ----------
def best_pair_interval(coords1, coords2, max_gap=30000):
    """
    Find a model1–model2 pair on same contig where the gap between features <= max_gap,
    and return the interval [min(s1,s2), max(e1,e2)] for the *largest span*.
    Returns (chrom, start, end) or None.
    """
    best_span = -1
    best_coords = None
    # Pre-group model2 by contig for speed
    by_chrom2 = {}
    for _, s2, e2, c2 in coords2:
        by_chrom2.setdefault(c2, []).append((s2, e2))
    for _, s1, e1, c1 in coords1:
        if c1 not in by_chrom2:
            continue
        for (s2, e2) in by_chrom2[c1]:
            span = max(e1, e2) - min(s1, s2)
            # distance between non-overlapping intervals (0 if overlapping)
            dist = 0
            if e1 < s2:
                dist = s2 - e1
            elif e2 < s1:
                dist = s1 - e2
            if dist <= max_gap and span > best_span:
                best_span = span
                best_coords = (c1, min(s1, s2), max(e1, e2))
    return best_coords

def fallback_from_model1_to_end(coords1, coords2, rec_by_id):
    """
    If there’s a model1 hit and *no downstream* model2 hit on that contig,
    return (chrom, start1, contig_length) for the *last* such model1 on that contig.
    Preference: the model1 hit with the greatest start (closest to end).
    Returns (chrom, start, end) or None.
    """
    # Index model2 starts per contig for downstream check
    m2_starts_by_c = {}
    for _, s2, _, c2 in coords2:
        m2_starts_by_c.setdefault(c2, []).append(s2)
    for c in m2_starts_by_c:
        m2_starts_by_c[c].sort()

    # Arrange model1 hits by contig, sorted by start
    m1_by_c = {}
    for _, s1, e1, c1 in coords1:
        m1_by_c.setdefault(c1, []).append((s1, e1))
    for c in m1_by_c:
        m1_by_c[c].sort(key=lambda t: t[0])  # sort by s1

    candidates = []
    for chrom, hits in m1_by_c.items():
        m2_starts = m2_starts_by_c.get(chrom, [])
        for (s1, e1) in hits:
            # Is there any model2 start >= e1 (downstream)?
            has_downstream_m2 = False
            # Binary-ish scan (list is small; linear ok)
            for s2 in m2_starts:
                if s2 >= e1:
                    has_downstream_m2 = True
                    break
            if not has_downstream_m2:
                # No model2 after this model1 on this contig → candidate
                contig_len = len(rec_by_id[chrom].seq)
                candidates.append((chrom, s1, contig_len))

    if not candidates:
        return None
    # Pick the candidate with the largest s1 (furthest to the right / closest to end)
    candidates.sort(key=lambda t: t[1], reverse=True)
    return candidates[0]

# ---------- Extraction ----------
def write_slice(rec, start, end, out_path, note=None):
    sub = rec[start:end]  # Biopython: 0-based, end-exclusive
    sub.id = rec.id
    sub.name = rec.name
    desc_extra = f" | region {start}-{end}"
    if note:
        desc_extra += f" | {note}"
    sub.description = (rec.description or rec.id) + desc_extra
    SeqIO.write(sub, out_path, "genbank")

# ---------- Main ----------
def process_files(model1, model2, protein_directory, genbank_directory, threshold=5, max_gap=30000):
    fasta_files = sorted(glob.glob(os.path.join(protein_directory, "*.faa")))
    genbank_files = glob.glob(os.path.join(genbank_directory, "*.gbff")) + \
                    glob.glob(os.path.join(genbank_directory, "*.gbk"))
    genbank_dict = {os.path.splitext(os.path.basename(f))[0]: f for f in genbank_files}
    print(f"Found {len(fasta_files)} FASTA files and {len(genbank_files)} GenBank files.")

    # Ensure coordinates log exists
    coords_log = "Locus-coordinates.txt"
    if not os.path.exists(coords_log):
        with open(coords_log, "w") as f:
            f.write("contig\tstart_0based\tend_0based\tsource\n")

    for fasta_file in fasta_files:
        base = os.path.splitext(os.path.basename(fasta_file))[0]
        if base not in genbank_dict:
            print(f"[{base}] No matching GenBank file found.")
            continue

        genbank_file = genbank_dict[base]
        print(f"[{base}] Processing {os.path.basename(fasta_file)} with {os.path.basename(genbank_file)}")

        # Unique temp files per sample
        dom1 = f"{base}.model1.domtblout"
        dom2 = f"{base}.model2.domtblout"
        log1 = f"{base}.model1.log"
        log2 = f"{base}.model2.log"

        # Run HMMER
        run_hmmer(model1, fasta_file, dom1, threshold=threshold, log_path=log1)
        run_hmmer(model2, fasta_file, dom2, threshold=threshold, log_path=log2)

        # Read HMM hits
        ids1 = read_hit_ids_from_domtbl(dom1)
        ids2 = read_hit_ids_from_domtbl(dom2)

        if not ids1:
            print(f"[{base}] No hits for model1; skipping.")
            cleanup([dom1, dom2, log1, log2])
            continue

        # Index GenBank & map IDs to coords
        rec_by_id, cds_index = index_genbank(genbank_file)
        coords1 = map_ids_to_coords(ids1, cds_index)
        coords2 = map_ids_to_coords(ids2, cds_index)

        if not coords1:
            print(f"[{base}] model1 hits did not map to GenBank CDS; skipping.")
            cleanup([dom1, dom2, log1, log2])
            continue

        # Try best pair within max_gap
        pair = best_pair_interval(coords1, coords2, max_gap=max_gap)

        if pair:
            chrom, start, end = pair
            out_gbk = f"{base}_{chrom}_{start}_{end}.gbk"
            write_slice(rec_by_id[chrom], start, end, out_gbk, note="model1+model2 interval")
            with open(coords_log, "a") as f:
                f.write(f"{chrom}\t{start}\t{end}\tmodel1+model2\n")
            print(f"[{base}] Extracted pair interval → {out_gbk}")
        else:
            # Fallback: model1 → end (if no downstream model2)
            fb = fallback_from_model1_to_end(coords1, coords2, rec_by_id)
            if fb:
                chrom, start, end = fb
                out_gbk = f"{base}_{chrom}_{start}_{end}.gbk"
                write_slice(rec_by_id[chrom], start, end, out_gbk, note="fallback model1→end")
                with open(coords_log, "a") as f:
                    f.write(f"{chrom}\t{start}\t{end}\tmodel1_to_end\n")
                print(f"[{base}] Extracted fallback interval (model1→end) → {out_gbk}")
            else:
                print(f"[{base}] No valid pair and no fallback (model1→end) candidate found.")

        cleanup([dom1, dom2, log1, log2])

def cleanup(paths):
    for p in paths:
        try:
            if os.path.exists(p):
                os.remove(p)
        except Exception:
            pass

# ---------- CLI ----------
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <model1.hmm> <model2.hmm> <faa_dir> <gb_dir>")
        sys.exit(1)
    model1 = sys.argv[1]
    model2 = sys.argv[2]
    protein_directory = sys.argv[3]
    genbank_directory = sys.argv[4]  # <-- fixed bug
    process_files(model1, model2, protein_directory, genbank_directory, threshold=5, max_gap=30000)

