.SUFFIXES: .plt .tex .pdf
 
%.tex:	dens_com_C.gp
	bash $< C
 
%.pdf: %.tex
	pdflatex $<
 
C: dens_com_C.pdf 
