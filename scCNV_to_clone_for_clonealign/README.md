First, run the wrapper script.

Within that is scCNV_to_clones_for_clonealign.R, which shows an example (on sample MP2B) of converting the long format data (with one line per intersection of node_cnv_calls.bed coordinate+cell with mappable regions) to wide, and using that to get a segment x clone matrix that will work well with clonealign.

The version of this script including plots is available in scCNV_to_clones_for_clonealign.html (output of Rmarkdown version of the R script).

Next, run clonealign. Example of this is in clonealign_example_MP2B.R.

Finally, compare clonealign results for two different sets of segments (one as run for the paper, another including an additional chromosome). This is in 
