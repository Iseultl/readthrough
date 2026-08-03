##############################################
# Stage 1 - Build geneid-rs
##############################################

FROM rust:1-slim-bookworm AS builder
WORKDIR /usr/src/app
COPY Cargo.toml Cargo.lock ./

COPY src ./src
COPY benches ./benches
RUN cargo build --release

FROM debian:bookworm-slim

ARG INPUT_DIR=/data/input
ARG OUTPUT_DIR=/data/output

ENV INPUT_DIR=${INPUT_DIR}
ENV OUTPUT_DIR=${OUTPUT_DIR}

RUN mkdir -p ${INPUT_DIR} ${OUTPUT_DIR} \
    && chmod 777 ${INPUT_DIR} ${OUTPUT_DIR}

COPY --from=builder /usr/src/app/target/release/geneid-rs /usr/local/bin/geneid-rs
ENTRYPOINT ["/usr/local/bin/geneid-rs"]


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