################## BASE IMAGE #####################
FROM nfcore/base


################## METADATA #######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="rnaseq-transcript-nf"
LABEL software.version="2.0"
LABEL about.summary="Container image containing all requirements for rnaseq-transcript-nf"
LABEL about.home="http://github.com/IARCbioinfo/RNAseq-transcript-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/RNAseq-transcript-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/RNAseq-transcript-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN conda env create -n rnaseq-transcript-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rnaseq-transcript-nf/bin:$PATH
RUN wget https://github.com/gpertea/stringtie/releases/download/v2.0.6/stringtie-2.0.6.Linux_x86_64.tar.gz && \
    tar -xzf stringtie-2.0.6.Linux_x86_64.tar.gz && \
    rm -rf /opt/conda/envs/rnaseq-transcript-nf/bin/stringtie && \
    rm -rf /opt/conda/envs/rnaseq-transcript-nf/bin/prepDE.py && \
    mv stringtie-2.0.6.Linux_x86_64/stringtie /usr/bin/. && \
    mv stringtie-2.0.6.Linux_x86_64/prepDE.py /usr/bin/. && \
    rm -rf stringtie-2.0.6.Linux_x86_64.tar.gz && \
    rm -rf stringtie-2.0.6.Linux_x86_64
