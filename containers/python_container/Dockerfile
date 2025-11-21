FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 

# Dependências mínimas (adicione libs do sistema se seus pacotes Python precisarem)
RUN apt-get update && apt-get install -y --no-install-recommends build-essential gcc && rm -rf /var/lib/apt/lists/*

# venv dedicado
RUN pip install --upgrade pip

# Instala em build-time (cacheável)
RUN pip install numpy==2.3.0 pandas==2.3.3 scikit-learn==1.7.2 pyreadr==0.5.3 joblib

# Onde o VS Code monta o projeto
WORKDIR /workspaces
