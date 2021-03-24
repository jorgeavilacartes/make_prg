# to build: sudo docker build . -t make_prg_precompiled_binary_builder:0.0.1
FROM centos:7
RUN yum -y update && yum install -y python3 python3-devel pip3
COPY ./requirements.txt /requirements.txt
RUN pip3 install -r requirements.txt