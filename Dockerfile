# ============================================================
# PCR Assay Design Workflow
#
# Docker environment for:
#   01_varvamp_assay_design.py
#   02_in_silico_validation.py
# ============================================================

FROM anaconda/miniconda:26.7.1

LABEL org.opencontainers.image.title="PCR Assay Design Workflow"
LABEL org.opencontainers.image.description="VarVAMP assay design and in silico validation workflow"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

ARG CONDA_ENV=varvamp-qpcr-workflow


# ------------------------------------------------------------
# 1. System dependencies
# ------------------------------------------------------------

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    gzip \
    unzip \
    ca-certificates \
    sudo \
    && rm -rf /var/lib/apt/lists/*


# ------------------------------------------------------------
# 2. Application directory
# ------------------------------------------------------------

WORKDIR /app


# ------------------------------------------------------------
# 3. Conda environment
# ------------------------------------------------------------

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml \
    && conda clean -afy


# ------------------------------------------------------------
# 4. MFEprimer 3.1.0
#
# Workflow 02 currently expects the MFEprimer 3.x
# text-report structure.
# ------------------------------------------------------------

RUN wget -q \
    https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.1.0/mfeprimer-3.1.0-linux-amd64.gz \
    -O /tmp/mfeprimer.gz \
    && gunzip /tmp/mfeprimer.gz \
    && chmod +x /tmp/mfeprimer \
    && mv /tmp/mfeprimer /usr/local/bin/mfeprimer \
    && mfeprimer -h >/dev/null


# ------------------------------------------------------------
# 5. MARS
#
# Circular-sequence refinement used by Workflow 01.
# ------------------------------------------------------------

RUN git clone --depth 1 https://github.com/lorrainea/MARS.git /tmp/MARS \
    && cd /tmp/MARS \
    && chmod +x pre-install.sh \
    && ./pre-install.sh \
    && make -f Makefile \
    && chmod +x mars \
    && cp mars /usr/local/bin/mars \
    && rm -rf /tmp/MARS


# ------------------------------------------------------------
# 6. Workflow source files
# ------------------------------------------------------------

COPY 01_varvamp_assay_design.py /app/01_varvamp_assay_design.py
COPY 02_in_silico_validation.py /app/02_in_silico_validation.py


# ------------------------------------------------------------
# 7. Runtime directories
# ------------------------------------------------------------

RUN mkdir -p \
    /app/data/design \
    /app/data/validation \
    /app/work \
    /app/results


# ------------------------------------------------------------
# 8. Build-time checks
# ------------------------------------------------------------

RUN conda run -n ${CONDA_ENV} python --version \
    && conda run -n ${CONDA_ENV} mafft --version \
    && conda run -n ${CONDA_ENV} seqkit version \
    && conda run -n ${CONDA_ENV} blastn -version \
    && conda run -n ${CONDA_ENV} which cd-hit-est \
    && conda run -n ${CONDA_ENV} which trimal \
    && conda run -n ${CONDA_ENV} which varvamp \
    && command -v mars \
    && command -v mfeprimer


# ------------------------------------------------------------
# 9. Default command
# ------------------------------------------------------------

CMD ["bash"]