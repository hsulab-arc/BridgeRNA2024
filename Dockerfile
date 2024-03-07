FROM continuumio/miniconda3

# Install dependencies
RUN apt-get -y update && apt-get -y install gcc make libbz2-dev  \
    zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev g++ python2

# Install conda and dependencies
RUN conda update -n base -c defaults conda && \
    conda install -y -c conda-forge python=3.8 mamba &&  \
    mamba install -y -c conda-forge -c bioconda snakemake samtools htslib blast muscle mmseqs2 && \
    mamba install -y -c conda-forge -c r r r-essentials r-tidyverse r-ggseqlogo r-patchwork

# Install mafft-qinsi
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.490-with-extensions-src.tgz && \
    tar -xvf mafft-7.490-with-extensions-src.tgz && \
    cd mafft-7.490-with-extensions/core && \
    make clean && \
    make install && \
    cd ../extensions && \
    make clean && \
    make && \
    make install && \
    cd ../../ && \
    rm mafft-7.490-with-extensions-src.tgz && \
    rm -rf mafft-7.490-with-extensions && \
    ln -s /usr/local/bin/mafft-qinsi /usr/bin/mafft-qinsi && \
    ln -s /usr/local/bin/mafft-linsi /usr/bin/mafft-linsi && \
    ln -s /usr/local/bin/mafft-einsi /usr/bin/mafft-einsi && \
    ln -s /usr/local/bin/mafft-ginsi /usr/bin/mafft-ginsi

# Install bridgerna2024
RUN echo "CACHEBREAKER0" && \
    pip install numpy Cython && \
    pip install git+https://github.com/hsulab-arc/BridgeRNA2024.git



