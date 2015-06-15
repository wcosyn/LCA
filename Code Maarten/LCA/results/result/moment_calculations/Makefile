.SUFFIXES: .plt .tex .pdf
 
%.tex:	%.plt
	gnuplot -e "set output '$@' \
	  " $<
 
%.pdf: %.tex
	pdflatex $<
 
C: moments_all.pdf 
