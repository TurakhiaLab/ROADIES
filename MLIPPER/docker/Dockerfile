ARG CUDA_VERSION=11.8.0
ARG UBUNTU_VERSION=22.04

FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu${UBUNTU_VERSION}

ARG LIBPLL_VERSION=0.4.0

SHELL ["/bin/bash", "-lc"]

ENV DEBIAN_FRONTEND=noninteractive
ENV CUDA_HOME=/usr/local/cuda
ENV CONDA_DIR=/opt/conda
ENV CONDA_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:/usr/local/cuda/bin:${PATH}
ENV PLL_INC_DIR=/usr/local/include
ENV PLL_LIB_DIR=/usr/local/lib
ENV EPA_NG_ENV=base
ENV RAXML_NG_ENV=base
# Add CUDA and /usr/local shared libraries to the runtime search path.
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/lib:${LD_LIBRARY_PATH}

WORKDIR /workspace/MLIPPER

RUN success=0; \
    for attempt in 1 2 3 4 5; do \
      if apt-get update && \
      apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    bison \
    flex \
    git \
    libblas-dev \
    libboost-program-options-dev \
    liblapack-dev \
    libtbb-dev \
    libtool \
    pkg-config \
    python3 \
    python3-pip \
    wget \
    zlib1g-dev; then \
      rm -rf /var/lib/apt/lists/*; \
      success=1; \
      break; \
      fi; \
      echo "apt install failed on attempt ${attempt}, retrying..." >&2; \
      rm -rf /var/lib/apt/lists/*; \
      sleep 5; \
    done; \
    test "$success" = 1

RUN git clone --quiet --branch "${LIBPLL_VERSION}" --depth 1 \
    https://github.com/xflouris/libpll-2.git /tmp/libpll-2 \
    && cd /tmp/libpll-2 \
    && ./autogen.sh \
    && ./configure --prefix=/usr/local \
    && make -j4 \
    && make install \
    && ldconfig \
    && rm -rf /tmp/libpll-2

RUN wget -qO /tmp/miniconda.sh \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash /tmp/miniconda.sh -b -p "${CONDA_DIR}" \
    && rm -f /tmp/miniconda.sh \
    && conda config --system --set auto_update_conda false \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
    && conda install -y -n base -c bioconda epa-ng raxml-ng aster \
    && conda run -n base epa-ng --help >/dev/null \
    && conda run -n base raxml-ng --help >/dev/null \
    && conda run -n base astral-pro3 --help >/dev/null \
    && conda clean -afy

COPY . /workspace/MLIPPER

RUN make USE_DOUBLE=1 MLIPPER

CMD ["/bin/bash"]
