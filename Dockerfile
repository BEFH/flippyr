FROM mambaorg/micromamba:2.0
RUN micromamba install -y -n base -c conda-forge -c bioconda -c bfh \
       python=3.12.7 \
       pyfaidx=0.8.1.3 \
       biopython=1.84 \
       flippyr=0.6.0 \
       plink=1.90b6.21 && \
    micromamba clean --all --yes
