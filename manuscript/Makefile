
TEX=sfbay-microplastic-transport-jgr.tex
PDF=sfbay-microplastic-transport-jgr.pdf
AUX=sfbay-microplastic-transport-jgr.aux

all: $(PDF)

$(PDF): $(TEX) references.bib
	pdflatex $(TEX)
	bibtex $(AUX)
	pdflatex $(TEX)
	pdflatex $(TEX)
