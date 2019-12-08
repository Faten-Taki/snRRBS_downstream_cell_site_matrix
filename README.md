#this repository has several files that should be used in the following order:

1. merge all CpG_sites_per_cells matrices from one or more single cell sequencing runs using the file "merge_cell_site_matrices.R"

2. get the differentially methylated sites (in this case, I have 2 cell groups/types of interest: Pos and Neg) using the file "get_differential_CpG_sites.R"

3. Visualize the cells relative to their methylation status across sites of interest using the file "visualize_using_tSNE.R"
