# hash:sha256:154c32122bd0a192d712de240b6a304865600f18b736eaa483e4357ab5b0ad23
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake=3.10.2-1ubuntu2.18.04.2 \
        libnlopt-dev=2.4.2+dfsg-4 \
        pkg-config=0.29.1-0ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("Matrix", "1.4-0")' \
    && Rscript -e 'remotes::install_version("arules", "1.7-3")' \
    && Rscript -e 'remotes::install_version("ca", "0.71.1")' \
    && Rscript -e 'remotes::install_version("factoextra", "1.0.7")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.3.5")'
