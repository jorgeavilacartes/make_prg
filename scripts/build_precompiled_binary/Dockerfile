# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg_precompiled_binary_builder:0.0.1
# This container has the python environment to build a static binary for make_prg
# To build: sudo docker build . -t leandroishilima/make_prg_precompiled_binary_builder:0.0.1
ARG ARCH=""
ARG ALPINE_VERSION="3.7"

FROM ${ARCH}python:${ALPINE_VERSION}-alpine

ARG PYINSTALLER_TAG
ENV PYINSTALLER_TAG ${PYINSTALLER_TAG:-"v3.6"}

# Official Python base image is needed or some applications will segfault.
# PyInstaller needs zlib-dev, gcc, libc-dev, and musl-dev
RUN apk --update --no-cache add \
    zlib-dev \
    musl-dev \
    libc-dev \
    libffi-dev \
    gcc \
    g++ \
    git \
    pwgen \
    && pip install --upgrade pip

# Install pycrypto so --key can be used with PyInstaller
RUN pip install \
    pycrypto

# Build bootloader for alpine with -no-pie
RUN git clone --depth 1 --single-branch --branch ${PYINSTALLER_TAG} https://github.com/pyinstaller/pyinstaller.git /tmp/pyinstaller \
    && cd /tmp/pyinstaller/bootloader \
    && CFLAGS="-Wno-stringop-overflow -Wno-stringop-truncation -no-pie" python ./waf configure --no-lsb all \
    && pip install .. \
    && rm -Rf /tmp/pyinstaller

VOLUME /src
WORKDIR /src

ADD ./bin /pyinstaller
RUN chmod a+x /pyinstaller/*

# install staticx
RUN apk add gfortran openblas-dev lapack-dev scons patchelf tk
RUN pip install staticx==0.10.0
# RUN sed -i -e 's/raise LibExistsError(basename(path))/pass/' /usr/local/lib/python3.7/site-packages/staticx/archive.py



ENTRYPOINT ["/pyinstaller/pyinstaller.sh"]
