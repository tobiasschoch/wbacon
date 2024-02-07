# Dockerfile used for R CMD check; based on r-base (official Docker container)
# This container is used to check wbacon
FROM r-base
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		libv8-dev \
        libcurl4-openssl-dev \
        r-cran-v8 \
		r-cran-rmarkdown \
		r-cran-knitr \
		r-cran-devtools \
		r-cran-robustbase \
        r-cran-hexbin \
		texinfo \
        texlive-fonts-extra \
        texlive-fonts-recommended \
        texlive-latex-extra \
        texlive-latex-recommended \
        texlive-luatex \
        texlive-plain-generic \
        texlive-science \
        texlive-xetex \
        lmodern \
        qpdf \
        ghostscript \
        tidy \
	&& Rscript -e 'install.packages("robustX")' \
	&& Rscript -e 'install.packages("modi")'
