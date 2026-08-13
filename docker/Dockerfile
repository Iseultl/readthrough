##############################################
# Stage 1 - Build geneid
##############################################

FROM debian:bookworm AS builder

ARG GENEID_VER=1.5.2

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        curl \
        ca-certificates

RUN curl -L \
    https://github.com/guigolab/geneid/archive/refs/tags/v1.5.2.tar.gz \
    -o geneid.tar.gz && \
    tar xzf geneid.tar.gz && \
    cd geneid-1.5.2 && \
    make

##############################################
# Stage 2 - Final runtime image
##############################################

FROM ileahy/biopythonimage:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    openjdk-17-jre \
    curl \
    wget \
    git \
    unzip \
    gzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

##################################################
# Install Nextflow
##################################################

RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/


##################################################
# Install annocli
##################################################

RUN pip install --upgrade pip

RUN git clone --branch fix_aliasMatch --depth 1 \
    https://github.com/apollo994/annocli.git /tmp/annocli

RUN pip install /tmp/annocli

##################################################
# Install geneid on command line
##################################################

COPY --from=builder geneid-1.5.2/bin/geneid /usr/local/bin/

RUN chmod +x /usr/local/bin/geneid

##################################################

WORKDIR /workspace

CMD ["/bin/bash"]