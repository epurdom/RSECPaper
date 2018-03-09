FIGURES= FinalFigures/*
BIBFILES= geneexpression.bib
all: bioarxiv_manuscript.pdf

OEAnalysis.Rout: OEAnalysis.R 
	Rscript $< > $@

bioarxiv_manuscript.pdf: bioarxiv_manuscript.tex OEAnalysisV2.Rout $(FIGURES) $(BIBFILES) 
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
