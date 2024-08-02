# hash:sha256:320569791615ad731d2bfea8a4b7cfeb382dab3f16fc2b9dd021ec86b00f6a04
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/jupyterlab:3.6.1-miniconda4.12.0-python3.9-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.8ubuntu1.1 \
        libgit2-dev=0.28.4+dfsg.1-2ubuntu0.1 \
        libssl-dev=1.1.1f-1ubuntu2.23 \
        pkg-config=0.29.1-0ubuntu4 \
        r-base=3.6.3-2 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        matplotlib==3.4.2 \
        pandas==2.2.2 \
        r==4.3.1 \
        r-irkernel==1.3.2 \
    && conda clean -ya

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && echo 'options(Ncpus = parallel::detectCores())' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && Rscript -e 'options(warn=2); install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("spatstat")'

COPY postInstall /
RUN /postInstall
