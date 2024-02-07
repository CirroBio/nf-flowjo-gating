FROM bioconductor/bioconductor_docker:RELEASE_3_17

ADD build.R /usr/local/
RUN Rscript /usr/local/build.R
