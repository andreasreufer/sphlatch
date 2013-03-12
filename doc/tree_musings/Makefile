all:
	touch _tmp.aux
	rm *.aux 2>&1 >/dev/null
	pdflatex tree_musings
	bibtex tree_musings
	pdflatex tree_musings
	pdflatex tree_musings
