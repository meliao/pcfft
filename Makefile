.ONESHELL:

NOTES_DIR=notes
spread-notes: 
	cd ${NOTES_DIR}
	pdflatex main.tex
	bibtex main.aux
	pdflatex main.tex
	pdflatex main.tex
	open main.pdf