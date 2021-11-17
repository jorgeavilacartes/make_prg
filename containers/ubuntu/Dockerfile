# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg:1.0.0
# To build: docker build . -t leandroishilima/make_prg:1.0.0
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y git python3 python3-dev python3-pip graphviz libgraphviz-dev mafft
RUN pip install git+https://github.com/leoisl/make_prg@update
