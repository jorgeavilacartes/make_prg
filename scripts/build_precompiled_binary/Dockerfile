# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg_precompiled_binary_builder:0.3.0
# This container has the python environment to build a static binary for make_prg
# To build: docker build . -t leandroishilima/make_prg_precompiled_binary_builder:0.3.0
FROM centos:7
RUN yum -y update && yum install -y python3 python3-devel pip3
COPY ./requirements.txt /requirements.txt
RUN pip3 install -r requirements.txt
