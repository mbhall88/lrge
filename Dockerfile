FROM rust:slim AS builder

COPY . /lrge

WORKDIR /lrge

# https://stackoverflow.com/a/76092743/5299417
RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99custom

RUN apt-get -y update \
    && apt-get upgrade -y \
    && apt install -y pkg-config build-essential zlib1g-dev  \
    && cargo build --release \
    && strip target/release/lrge

FROM ubuntu:jammy

COPY --from=builder /lrge/target/release/lrge /bin/

RUN lrge --version

ENTRYPOINT [ "lrge" ]
