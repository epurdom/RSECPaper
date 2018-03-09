This is a Github repository to accompany the paper "clusterExperiment and RSEC: A Bioconductor package and framework for clustering of single-cell and other large gene expression datasets"

# Reproduce the manuscript

## From scratch

The provided make file will allow you to recreate the analysis in the paper from scratch. This code includes downloading the data from GEO. 

You should first have R installed and the necessary packages installed: (`BiocInstaller`,`clusterExperiment`,`benchmarkme`,and `scone`)
).

Then type the following commands:

```
make bioarxiv_manuscript.pdf
```

Note that this is quite computationally intensive, and runs on parallel cores (by default the value defined by the environment variable `SLURM_CPUS_PER_TASK` or 6 if such a value is missing). So it should only be run from scratch in a setting appropriate to this.

## Using the provided output (quicker)

To have quicker access to the data and results, we also provide (via `git lfs`) `.rda` files that are the result of the output of downloading the data, and the `RSEC` command (see [https://git-lfs.github.com/](https://git-lfs.github.com/) for instructions on `git lfs`). They are saved under the directory `dataOutput_submitted`. 

To be able to run the above make command, and not have to rerun the computationally intensive commands, you must make the following changes:

1. Move (or copy) the directory  `dataOutput_submitted` into `dataOutput`
2. Edit the file `OEAnalysis.R`  so that the line
	```
	runClus<-FALSE
	```
	
	now reads as 
	
	```
	runClus<-TRUE
	```
After these changes you should be able to run the above make command,
```
make bioarxiv_manuscript.pdf
```

With these settings, the scripts should be able to run on a laptop relatively quickly. 
