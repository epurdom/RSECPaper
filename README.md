This is a Github repository to accompany the paper:

Risso D, Purvis L, Fletcher R, Das D, Ngai J, Dudoit S, Purdom E (2018) "clusterExperiment and RSEC: A Bioconductor package and framework for clustering of single-cell and other large gene expression datasets" PLoS Comput Biol. 2018 Sep 4;14(9):e1006378  http://dx.plos.org/10.1371/journal.pcbi.1006378

# How to Reproduce the analysis

## From scratch

The provided make file will allow you to recreate the analysis in the paper from scratch. This code includes downloading the data from GEO. 

You should first have R installed and the necessary packages installed: (`BiocInstaller`,`clusterExperiment`,`benchmarkme`,and `scone`)
).

Then type the following commands:

```
make OEAnalysis.Rout
```

Note that this is quite computationally intensive, and runs on parallel cores (by default the value defined by the environment variable `SLURM_CPUS_PER_TASK` or 6 if such a value is missing). So it should only be run from scratch in a setting appropriate to this.

## Using the provided intermediate results (quicker)

To have quicker access to the data and results, we also provide (via `git lfs`) `.rda` files that are the result of the output of downloading the data, and the `RSEC` command (see [https://git-lfs.github.com/](https://git-lfs.github.com/) for instructions on `git lfs` and how to download these files). They are saved under the directory `dataOutput_submitted`. 

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
make OEAnalysis.Rout
```

With these settings, the scripts should be able to run on a laptop relatively quickly. 
