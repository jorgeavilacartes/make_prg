# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg_precompiled_binary_builder:1.0.0
# This container has the python environment to build a static binary for make_prg
# To build: docker build . -t leandroishilima/make_prg_precompiled_binary_builder:1.0.0
# TODO: need to downgrade packages to use older versions of linux
# TODO: this is not good at all, use https://github.com/six8/pyinstaller-alpine instead
# TODO: or don't use precompiled binary at all
FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y python3 python3-dev python3-pip graphviz libgraphviz-dev zlib1g-dev libjpeg-dev
COPY ./requirements.txt /requirements.txt
RUN python3 -m pip install -r requirements.txt
RUN cp -r /usr/bin/* /usr/sbin/
