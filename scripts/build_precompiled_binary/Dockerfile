# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg_precompiled_binary_builder:1.0.0
# This container has the python environment to build a static binary for make_prg
# To build: docker build . -t leandroishilima/make_prg_precompiled_binary_builder:1.0.0
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y python3 python3-dev python3-pip graphviz libgraphviz-dev
COPY ./requirements.txt /requirements.txt
RUN python3 -m pip install -r requirements.txt
RUN cp -r /usr/bin/* /usr/sbin/
