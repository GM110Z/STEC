python plot_two_stx_similarity_heatmaps.py \
  stx_coord_clusters.all_members_of_filtered_stx_clusters.tsv \
  stx_annotation_reps.representatives.tsv \
  proteome_similarity_matrix.tsv \
  --out_prefix stx_similarity_reds \
  --cmap Reds \
  --vmin 0 --vmax 1
