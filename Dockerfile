# Miniconda is part of the final image. Digest pins the linux/amd64 base.
FROM anaconda/miniconda:26.7.1@sha256:e73e49c19eb003fa7268e565950119e696953343c5ac20b626f81e1b04b06652 AS base
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN test "$(uname -m)" = x86_64

# Build MARS separately so compilers and source trees stay out of the final image.
FROM base AS tools-builder
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git curl unzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /tmp/mars
RUN curl -fsSL --retry 3 \
    https://codeload.github.com/lorrainea/MARS/tar.gz/cbf8f594e96db8f757f1d84edf773406b4c72701 \
    -o /tmp/mars.tar.gz \
    && echo 'a9c2fbc5b1e9042ba04e3c40729c7513c0cb039f284aaa8de1e2710636274e19  /tmp/mars.tar.gz' | sha256sum -c - \
    && tar -xzf /tmp/mars.tar.gz --strip-components=1 \
    && bash -e pre-install.sh \
    && make -j2 \
    && install -Dm755 mars /out/bin/mars \
    && mkdir -p /out/lib \
    && cp -a libsdsl/lib/. /out/lib/
RUN curl -fsSL --retry 3 \
    https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.1.0/mfeprimer-3.1.0-linux-amd64.gz \
    -o /tmp/mfeprimer.gz \
    && echo 'ee3962e38993465647cea615d6ea5d24682eb9bc8d0e36706693351231d650ec  /tmp/mfeprimer.gz' | sha256sum -c - \
    && gzip -dc /tmp/mfeprimer.gz > /out/bin/mfeprimer \
    && chmod 755 /out/bin/mfeprimer

FROM base AS runtime
LABEL org.opencontainers.image.title="PCR Assay Design Workflow" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"
# MARS uses the C++ and OpenMP runtimes. MAFFT's shell/Perl helpers are
# supplied by its Conda package dependencies and the Linux base.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libstdc++6 ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY --from=tools-builder /out/bin/ /usr/local/bin/
COPY --from=tools-builder /out/lib/ /usr/local/lib/
RUN ldconfig
COPY environment.yml /tmp/environment.yml
# Flexible priority preserves the workflow's MAFFT 7.505 bioconda pin when
# conda-forge publishes only newer MAFFT versions. Defaults remain excluded.
RUN conda config --system --remove channels defaults \
    && conda config --system --add channels bioconda \
    && conda config --system --add channels conda-forge \
    && conda config --system --set channel_priority flexible \
    && conda env create --file /tmp/environment.yml \
    && conda clean --all --yes \
    && rm -rf /root/.cache /tmp/environment.yml
ENV PYTHONUNBUFFERED=1 MPLBACKEND=Agg
WORKDIR /app
COPY 01_varvamp_assay_design.py 02_in_silico_validation.py /app/
RUN mkdir -p data/design data/validation work results /opt/workflow-manifests

# Every check runs inside the actual runtime environment; any failure aborts.
SHELL ["conda", "run", "--no-capture-output", "-n", "pcr-assay-design-workflow", "/bin/bash", "-o", "pipefail", "-c"]
RUN set -eu; \
    python --version; \
    mafft --version; \
    seqkit version; \
    blastn -version; \
    makeblastdb -version; \
    for tool in cd-hit-est trimal varvamp mars mfeprimer; do command -v "$tool"; done; \
    trimal --version; \
    varvamp --help >/dev/null; \
    mfeprimer -h >/dev/null; \
    if ldd /usr/local/bin/mars | grep -q 'not found'; then exit 1; fi; \
    python -c 'import sys, Bio, pandas, matplotlib.pyplot, PIL.Image, fitz, numpy, primer3, seqfold, varvamp; assert sys.version_info[:2] == (3, 11)'; \
    python -m pip check; \
    python3 01_varvamp_assay_design.py --help; \
    python3 02_in_silico_validation.py --help; \
    conda list --explicit > /opt/workflow-manifests/conda-explicit.txt; \
    python -m pip freeze > /opt/workflow-manifests/pip-freeze.txt; \
    dpkg-query -W > /opt/workflow-manifests/dpkg-packages.txt
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "pcr-assay-design-workflow"]
CMD ["bash"]
