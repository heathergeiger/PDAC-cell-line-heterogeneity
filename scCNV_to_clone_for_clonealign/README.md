First, run the wrapper script.

Within that is scCNV_to_clones_for_clonealign.R, which shows an example (on sample MP2B) of converting the long format data (with one line per intersection of node_cnv_calls.bed coordinate+cell with mappable regions) to wide, and using that to get a segment x clone matrix that will work well with clonealign.

The version of this script including plots is available in scCNV_to_clones_for_clonealign.html (output of Rmarkdown version of the R script).

Next, run clonealign. Example of this is in clonealign_example_MP2B.R.

Finally, compare clonealign results for two different sets of segments (one as run for the paper, another including an additional chromosome). This is in MP2B_scRNA_UMAP_and_clustering_vs_clonealign_clones.R and MP2B_scRNA_UMAP_and_clustering_vs_clonealign_clones.html. 

I also included here some small intermediate files such as the coordinates per annotated gene bed file, the segment x clone matrices output by scCNV_to_clones_for_clonealign.R, and the clonealign calls output by clonealign_example_MP2B.R.
