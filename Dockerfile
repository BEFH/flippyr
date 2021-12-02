FROM mambaorg/micromamba:0.17.0
RUN micromamba install -y -n base -c conda-forge -c bioconda -c bfh \
       python=3.9.7 \
       biopython=1.79 \
       flippyr=0.5.2 \
       plink=1.90b6.21 && \
    micromamba clean --all --yes
