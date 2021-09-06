FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update \
 && apt-get -y install build-essential cmake

COPY . /fmmi/

RUN mkdir -p /fmmi/build \
 && cd /fmmi/build \
 && cmake .. \
 && make -j4

CMD /fmmi/build/test/doctest
