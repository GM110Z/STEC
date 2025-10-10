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
            ids = [i.split(".", 1)[0] for i in ids]  # strip version
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
def best_pair_intervals(coordsA, coordsB, max_gap=30000):
    """
    Find all A–B pairs on same contig with gap <= max_gap.
    Returns list of (chrom, start, end), de-duplicated.
    """
    if not coordsA or not coordsB:
        return []
    intervals = set()
    by_chromB = {}
    for _, s2, e2, c2 in coordsB:
        by_chromB.setdefault(c2, []).append((s2, e2))
    for _, s1, e1, c1 in coordsA:
        if c1 not in by_chromB:
            continue
        for (s2, e2) in by_chromB[c1]:
            # distance between non-overlapping intervals (0 if overlapping)
            if e1 < s2:
                dist = s2 - e1
            elif e2 < s1:
                dist = s1 - e2
            else:
                dist = 0
            if dist <= max_gap:
                start = min(s1, s2)
                end = max(e1, e2)
                intervals.add((c1, start, end))
    return sorted(intervals)

def fallback_from_model1_to_end(coords1, coords_others, rec_by_id):
    """
    If there’s a model1 hit and *no downstream* hit from ANY other model on that contig,
    return (chrom, start1, contig_length) for the *last* such model1 on that contig.
    """
    # Index all "other" starts per contig
    other_starts_by_c = {}
    for (_pid, s, _e, c) in (coords_others or []):
        other_starts_by_c.setdefault(c, []).append(s)
    for c in other_starts_by_c:
        other_starts_by_c[c].sort()

    # Arrange model1 hits by contig, sorted by start
    m1_by_c = {}
    for _, s1, e1, c1 in (coords1 or []):
        m1_by_c.setdefault(c1, []).append((s1, e1))
    for c in m1_by_c:
        m1_by_c[c].sort(key=lambda t: t[0])

    candidates = []
    for chrom, hits in m1_by_c.items():
        other_starts = other_starts_by_c.get(chrom, [])
        for (s1, e1) in hits:
            has_downstream_other = False
            for s_other in other_starts:
                if s_other >= e1:
                    has_downstream_other = True
                    break
            if not has_downstream_other:
                contig_len = len(rec_by_id[chrom].seq)
                candidates.append((chrom, s1, contig_len))

    if not candidates:
        return None
    candidates.sort(key=lambda t: t[1], reverse=True)
    return candidates[0]

# ---------- Extraction ----------
def write_slice(rec, start, end, out_path, note=None):
    sub = rec[start:end]  # 0-based, end-exclusive
    sub.id = rec.id
    sub.name = rec.name
    desc_extra = f" | region {start}-{end}"
    if note:
        desc_extra += f" | {note}"
    sub.description = (rec.description or rec.id) + desc_extra
    SeqIO.write(sub, out_path, "genbank")

# ---------- Main ----------
def process_files(model1, model2, model3, protein_directory, genbank_directory, threshold=5, max_gap=30000):
    fasta_files = sorted(glob.glob(os.path.join(protein_directory, "*.faa")))
    genbank_files = glob.glob(os.path.join(genbank_directory, "*.gbff")) + \
                    glob.glob(os.path.join(genbank_directory, "*.gbk"))
    genbank_dict = {os.path.splitext(os.path.basename(f))[0]: f for f in genbank_files}
    print(f"Found {len(fasta_files)} FASTA files and {len(genbank_files)} GenBank files.")

    coords_log = "Locus-coordinates.txt"
    if not os.path.exists(coords_log):
        with open(coords_log, "w") as f:
            f.write("sample\tstart_0based\tend_0based\tsource\n")

    for fasta_file in fasta_files:
        base = os.path.splitext(os.path.basename(fasta_file))[0]
        if base not in genbank_dict:
            print(f"[{base}] No matching GenBank file found.")
            continue

        genbank_file = genbank_dict[base]
        print(f"[{base}] Processing {os.path.basename(fasta_file)} with {os.path.basename(genbank_file)}")

        # Unique temp files per sample
        dom1 = f"{base}.model1.domtblout"; log1 = f"{base}.model1.log"
        dom2 = f"{base}.model2.domtblout"; log2 = f"{base}.model2.log"
        run_hmmer(model1, fasta_file, dom1, threshold=threshold, log_path=log1)
        run_hmmer(model2, fasta_file, dom2, threshold=threshold, log_path=log2)

        ids1 = read_hit_ids_from_domtbl(dom1)
        ids2 = read_hit_ids_from_domtbl(dom2)

        # Optional model3
        ids3 = []
        dom3 = log3 = None
        if model3 is not None:
            dom3 = f"{base}.model3.domtblout"; log3 = f"{base}.model3.log"
            run_hmmer(model3, fasta_file, dom3, threshold=threshold, log_path=log3)
            ids3 = read_hit_ids_from_domtbl(dom3)

        if not ids1:
            print(f"[{base}] No hits for model1; skipping.")
            cleanup([p for p in [dom1, dom2, dom3, log1, log2, log3] if p])
            continue

        # Index GenBank & map IDs to coords
        rec_by_id, cds_index = index_genbank(genbank_file)
        coords1 = map_ids_to_coords(ids1, cds_index)
        coords2 = map_ids_to_coords(ids2, cds_index)
        coords3 = map_ids_to_coords(ids3, cds_index) if ids3 else []

        if not coords1:
            print(f"[{base}] model1 hits did not map to GenBank CDS; skipping.")
            cleanup([p for p in [dom1, dom2, dom3, log1, log2, log3] if p])
            continue

        # Find intervals: (1,2) always; (1,3) only if model3 provided & mapped
        intervals_12 = best_pair_intervals(coords1, coords2, max_gap=max_gap)
        intervals_13 = best_pair_intervals(coords1, coords3, max_gap=max_gap) if coords3 else []

        extracted_any = False

        # Extract all 1+2
        for chrom, start, end in intervals_12:
            out_gbk = f"{base}_{chrom}_{start}_{end}_1plus2.gbk"
            write_slice(rec_by_id[chrom], start, end, out_gbk, note="model1+model2 interval")
            with open(coords_log, "a") as f:
                f.write(f"{base}\t{start}\t{end}\tmodel1+model2\n")
            print(f"[{base}] Extracted 1+2 interval → {out_gbk}")
            extracted_any = True

        # Extract all 1+3 (if any)
        for chrom, start, end in intervals_13:
            out_gbk = f"{base}_{chrom}_{start}_{end}_1plus3.gbk"
            write_slice(rec_by_id[chrom], start, end, out_gbk, note="model1+model3 interval")
            with open(coords_log, "a") as f:
                f.write(f"{base}\t{start}\t{end}\tmodel1+model3\n")
            print(f"[{base}] Extracted 1+3 interval → {out_gbk}")
            extracted_any = True

        # Fallback: if nothing extracted, use model1→end when no downstream partner (2 or 3)
        if not extracted_any:
            coords_others = (coords2 or []) + (coords3 or [])
            fb = fallback_from_model1_to_end(coords1, coords_others, rec_by_id)
            if fb:
                chrom, start, end = fb
                out_gbk = f"{base}_{chrom}_{start}_{end}_fallback.gbk"
                write_slice(rec_by_id[chrom], start, end, out_gbk, note="fallback model1→end (no downstream model2/3)")
                with open(coords_log, "a") as f:
                    f.write(f"{base}\t{start}\t{end}\tmodel1_to_end\n")
                print(f"[{base}] Extracted fallback interval (model1→end) → {out_gbk}")
            else:
                print(f"[{base}] No valid intervals and no fallback candidate found.")

        cleanup([p for p in [dom1, dom2, dom3, log1, log2, log3] if p])

def cleanup(paths):
    for p in paths:
        try:
            if os.path.exists(p):
                os.remove(p)
        except Exception:
            pass

# ---------- CLI ----------
if __name__ == "__main__":
    # Accept 4 or 5 positional args after script name:
    # 2-HMMs: script.py m1 m2 faa_dir gb_dir
    # 3-HMMs: script.py m1 m2 m3 faa_dir gb_dir
    if len(sys.argv) == 5:
        model1, model2, protein_directory, genbank_directory = sys.argv[1:5]
        model3 = None
    elif len(sys.argv) == 6:
        model1, model2, model3, protein_directory, genbank_directory = sys.argv[1:6]
    else:
        print("Usage:\n  python script.py <model1.hmm> <model2.hmm> <faa_dir> <gb_dir>\n"
              "  python script.py <model1.hmm> <model2.hmm> <model3.hmm> <faa_dir> <gb_dir>")
        sys.exit(1)

    process_files(model1, model2, model3, protein_directory, genbank_directory, threshold=5, max_gap=45000)

