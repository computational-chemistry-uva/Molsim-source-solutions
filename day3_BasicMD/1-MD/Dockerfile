FROM ubuntu:22.04

RUN \
  apt-get update && \
  apt-get install -y --no-install-recommends llvm clang cmake make && \
  apt-get install -y --no-install-recommends python3 pybind11-dev python3-pybind11 python3-numpy python3-matplotlib && \
  apt-get clean

WORKDIR /app

COPY . /app

RUN \
  cmake . && \
  cmake --build . && \
  python3 test.py
