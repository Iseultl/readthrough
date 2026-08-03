##############################################
# Stage 1 - Build geneid-rs
##############################################

FROM rust:1.89-bookworm AS geneid-builder

RUN apt-get update && apt-get install -y \
    git \
    pkg-config \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN git clone --depth 1 \
    https://gitlab.com/bioinformatics6495415/bioinf-rust/genid-rs.git

WORKDIR /build/genid-rs

RUN cargo build --release


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
    python3 \
    python3-pip \
    python-is-python3 \
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

RUN python -m pip install --no-cache-dir \
    biopython \
    matplotlib \
    numpy \
    pandas

RUN git clone --branch v1.0.0 --depth 1 \
    https://github.com/apollo994/annocli.git /tmp/annocli

RUN pip install /tmp/annocli

##################################################
# Install geneid-rs
##################################################

COPY --from=geneid-builder \
    /build/genid-rs/target/release/geneid-rs \
    /usr/local/bin/

RUN chmod +x /usr/local/bin/geneid-rs

##################################################

WORKDIR /workspace

CMD ["/bin/bash"]