echo "accession\tassembly_level\tcontig_count\tcontig_n50\tscaffold_count\tscaffold_n50" > assembly_metrics.tsv

while read acc; do
    datasets summary genome accession $acc \
    | jq -r '
    [
      .reports[0].accession,
      .reports[0].assembly_info.assembly_level,
      .reports[0].assembly_stats.number_of_contigs,
      .reports[0].assembly_stats.contig_n50,
      .reports[0].assembly_stats.number_of_scaffolds,
      .reports[0].assembly_stats.scaffold_n50
    ] | @tsv'
done < $1 >> assembly_metrics.tsv
