FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.2.0 as scripts

FROM ubuntu:jammy

LABEL MAINTAINER="Heather Ward <heather@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG CONDA_ENVIRONMENT "${CONDA_ENVIRONMENT}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		time \
		xxd \
		python3-dev

ARG CELLRANGER_VERSION
# N.B. this must be passed as a build arg (using the utility script: -b CELLRANGER_DOWNLOAD_QUERY_STRING='xxyy')
#   Enter information here: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
#   The value of CELLRANGER_DOWNLOAD_QUERY_STRING = everything following the ? in the cURL request URL (no quotes)
ARG CELLRANGER_DOWNLOAD_QUERY_STRING
RUN wget -O cellranger-${CELLRANGER_VERSION}.tar.gz \
	"https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VERSION}.tar.gz?${CELLRANGER_DOWNLOAD_QUERY_STRING}" \
	&& tar -zxvf cellranger-${CELLRANGER_VERSION}.tar.gz --directory /opt \
	&& rm cellranger-${CELLRANGER_VERSION}.tar.gz

RUN ln -s /opt/cellranger-${CELLRANGER_VERSION}/cellranger /usr/local/bin/

COPY scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
