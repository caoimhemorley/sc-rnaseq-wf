ARG CUDA_VERSION
ARG UBUNTU_VERSION
FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.2.0 as scripts

FROM nvcr.io/nvidia/cuda:${CUDA_VERSION}-base-ubuntu${UBUNTU_VERSION}
ENV CUDA_VERSION "${CUDA_VERSION}"
ENV UBUNTU_VERSION "${UBUNTU_VERSION}"

LABEL MAINTAINER="Karen Fang <karen@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG JAX_VERSION
ENV JAX_VERSION "${JAX_VERSION}"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update \
	&& apt-get -qq install \
		build-essential \
		wget \
		time \
		xxd \
		curl \
		zlib1g-dev \
		libncursesw5-dev \
		libssl-dev \
		libsqlite3-dev \
		tk-dev \
		libgdbm-dev \
		libc6-dev \
		libbz2-dev \
		libffi-dev \
		liblzma-dev \
		git \
	&& rm -rf /var/lib/apt/lists/*

# gcloud sdk; needed to upload output files
ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
RUN wget "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" \
	&& tar -zxvf "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" --directory /opt \
	&& rm "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz"

ENV PATH "${PATH}:/opt/google-cloud-sdk/bin"

ARG PYTHON3_VERSION
ENV PYTHON3_VERSION "${PYTHON3_VERSION}"
RUN curl -O https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tar.xz \
	&& tar -xvf Python-${PYTHON3_VERSION}.tar.xz --directory /opt/ \
	&& rm Python-${PYTHON3_VERSION}.tar.xz
RUN cd /opt/Python-${PYTHON3_VERSION} \
	&& ./configure \
	&& make \
	&& make altinstall

RUN update-alternatives --install /usr/bin/python3 python3 /opt/Python-${PYTHON3_VERSION}/python 1 \
	&& update-alternatives --set python3 /opt/Python-${PYTHON3_VERSION}/python

# Cell Type Mapper
ARG CELL_TYPE_MAPPER_VERSION
ENV CELL_TYPE_MAPPER_VERSION "${CELL_TYPE_MAPPER_VERSION}"
RUN git clone --branch v${CELL_TYPE_MAPPER_VERSION} https://github.com/AllenInstitute/cell_type_mapper.git \
	&& cd cell_type_mapper \
	&& python3 -m pip install -e .

COPY ./requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt

# CUDA-enabled jaxlib is needed
RUN python3 -m pip install --upgrade "jax[cuda12_pip]==${JAX_VERSION}" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Torch version 2.6.0 needs libcudnn.so.9
RUN python3 -m pip install nvidia-cudnn-cu12==9.1.0.70

COPY --from=scripts /opt/scripts /opt/scripts
COPY scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts:/opt/scripts/main"

COPY resources /opt/resources
