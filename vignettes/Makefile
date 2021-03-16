DOCUMENT=wbacon
R=R.exe
PDFLATEX=pdflatex.exe

all: $(DOCUMENT).tex
	$(PDFLATEX) $(DOCUMENT).tex
	$(PDFLATEX) $(DOCUMENT).tex

$(DOCUMENT).tex:
	$(R) CMD Sweave $(DOCUMENT).Rnw

.PHONY: clean
clean:
	rm *.aux *.log *.out *.tex
	ls -A | grep "pdf" | grep -v "$(DOCUMENT).pdf" | xargs rm
