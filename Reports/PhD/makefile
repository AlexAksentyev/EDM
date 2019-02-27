OBJECTS = edm_experiments systematic_error_analysis fs_and_qfs_lattices
TARGETS = $(addsuffix .pdf, $(OBJECTS))
SOURCES = $(addsuffix .tex, $(OBJECTS))


all: $(TARGETS)

%.pdf: %.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $<

clean:
	latexmk -c

.PHONY: clean 
