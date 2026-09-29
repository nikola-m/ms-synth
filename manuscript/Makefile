# Build the manuscript and Additional file 1.  `make submission` moves
# figures and tables to the end of the document, as GigaScience requires.
TEX = pdflatex -interaction=nonstopmode
.PHONY: all figs manuscript supplement submission clean
all: figs manuscript supplement
figs:
	mkdir -p figs && cp ../results/figures/*.pdf figs/
manuscript: figs
	$(TEX) ms_synth_manuscript && bibtex ms_synth_manuscript && $(TEX) ms_synth_manuscript && $(TEX) ms_synth_manuscript
supplement:
	$(TEX) ms_synth_additional_file1 && $(TEX) ms_synth_additional_file1
submission: figs
	sed 's/\\endfloatsfalse/\\endfloatstrue/' ms_synth_manuscript.tex > ms_synth_manuscript_submission.tex
	$(TEX) ms_synth_manuscript_submission && bibtex ms_synth_manuscript_submission && $(TEX) ms_synth_manuscript_submission && $(TEX) ms_synth_manuscript_submission
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.fff *.ttt *.lof *.lot
