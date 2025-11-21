# syntax=docker/dockerfile:1

########################
# Alvo 1: base (Nextflow)
########################
FROM rocker/r-ver:4.5.1 AS base

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# deps mínimas p/ compilar pacotes R comuns
RUN apt-get update && apt-get install -y --no-install-recommends \
    nano curl git ca-certificates locales \
    build-essential libssl-dev libcurl4-openssl-dev libxml2-dev \
    # (opcionais p/ gráficos: libfontconfig1 libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev) \
 && sed -i 's/# en_US.UTF-8/en_US.UTF-8/' /etc/locale.gen && locale-gen \
 && rm -rf /var/lib/apt/lists/*

# pak (instalação rápida) + seus pacotes
RUN R -q -e "install.packages('pak', repos='https://cloud.r-project.org')" && \
    R -q -e "pak::pkg_install(c('tidyverse','data.table','remotes','future.apply','ggpubr','VIM','forecast','ape','WCluster','ggExtra','ggtree','patchwork','gghalves','randomForest','ggrepel','ggtext','extrafont','RColorBrewer'))"

# usuário não-root (boa prática p/ Nextflow)
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} nxf && useradd -m -u ${UID} -g ${GID} nxf
USER nxf
WORKDIR /workspace
CMD ["R","--version"]

########################
# Alvo 2: rstudio (temporário p/ debug)
########################
FROM rocker/rstudio:4.5.1 AS rstudio

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# deps para compilar pacotes (inclui curl headers) + locale
RUN apt-get update && apt-get install -y --no-install-recommends \
    nano curl git ca-certificates locales \
    build-essential pkg-config \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
    # (opcionais, úteis para gráficos/ggplot/ragg)
    libfontconfig1 libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
 && sed -i 's/# en_US.UTF-8/en_US.UTF-8/' /etc/locale.gen && locale-gen \
 && rm -rf /var/lib/apt/lists/*

# pak + mesmos pacotes do alvo base (mantém consistência)
RUN R -q -e "install.packages('pak', repos='https://cloud.r-project.org')" && \
    R -q -e 'pak::pkg_install(c("tidyverse","data.table","remotes","future.apply","ggpubr","VIM","forecast","ape","WCluster","ggExtra","ggtree","patchwork","gghalves","randomForest","ggrepel","ggtext","extrafont","RColorBrewer"))'

WORKDIR /home/rstudio/project
EXPOSE 8787
