##############################################
# Stage 1 - Build geneid-rs
##############################################

FROM debian:stretch AS builder

ARG GENEID_VER=1.5.2

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        curl \
        ca-certificates \
        build-essential

RUN curl -O https://github.com/guigolab/geneid/archive/refs/tags/v1.5.2.tar.gz && \
    tar xz && \
    cd geneid-${GENEID_VER} && \
    make BIN=/build

FROM debian:stretch
COPY --from=builder /build/geneid /usr/local/bin/

ENTRYPOINT [ "/usr/local/bin/geneid" ]

##############################################
# Stage 2 - Final runtime image
##############################################

FROM ileahy/biopythonimage:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    openjdk-17-jre \
    curl \
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

RUN git clone --branch v1.0.0 --depth 1 \
    https://github.com/apollo994/annocli.git /tmp/annocli

RUN pip install /tmp/annocli

##################################################
# Install geneid-rs
##################################################

COPY --from=builder \
    /build/geneid /usr/local/bin/

RUN chmod +x /usr/local/bin/geneid

##################################################

WORKDIR /workspace

CMD ["/bin/bash"]