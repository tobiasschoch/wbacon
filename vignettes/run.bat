Rscript -e "Sweave('wbacon.Rnw')" 
pdflatex wbacon.tex 
bibtex wbacon
pdflatex wbacon.tex 
rm wbacon.aux wbacon.blg wbacon.log wbacon.out wbacon.tex 
