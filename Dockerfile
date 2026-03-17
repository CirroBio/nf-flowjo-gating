FROM continuumio/miniconda3

RUN conda install -y -c bioconda -c conda-forge \
    bioconductor-cytoml \
    bioconductor-cytolib \
    bioconductor-flowworkspace \
    r-tidyr \
    r-devtools && \
    conda clean -afy

RUN pip3 install anndata pandas plotly scipy
