.SUFFIXES: .plt .tex .pdf
 
%.tex:	%.plt
	gnuplot -e "set output '$@' \
	  " $<
 
%.pdf: %.tex
	pdflatex $<
 
ratio: dens_rel_ratio_Fe3.pdf

all: dens_rel_all_Fe3.pdf
