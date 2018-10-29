

%.pdf: %.tex
	pdflatex $<
	bibtex $($<: .tex=.aux)
	pdflatex $<
	pdflatex $<
