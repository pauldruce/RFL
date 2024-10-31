# Stage 1: Base Image
FROM ubuntu:24.04 AS base

RUN apt-get update && apt-get upgrade -y

# Stage 2: CMake Image
FROM base AS cmake
RUN apt-get install -y \
    cmake \
    clang \
    make \
    gcc \
    g++ \
    libc-dev \
    libarmadillo-dev \
    libgsl-dev          \
    git \
    doxygen             \
    graphviz            \
    libhdf5-serial-dev

# Stage 3: Builder
FROM cmake AS rfl_builder

COPY ./  /RFL/

WORKDIR /RFL

RUN cmake -B ./build . \
    && cmake --build ./build --target all -j 4

RUN ( cd build && ctest -j 2 )\
    && cmake --install ./build

# Stage 4: Final
FROM cmake AS full

# Create a user for the container
RUN useradd -m -s /bin/bash fuzzyuser
USER fuzzyuser
WORKDIR /home/fuzzyuser

COPY --from=rfl_builder --chown=fuzzyuser:fuzzyuser /RFL/lib/new_RFL/include /usr/include/RFL/
COPY --from=rfl_builder --chown=fuzzyuser:fuzzyuser /RFL/lib/new_RFL/bin/ /usr/lib/RFL

COPY --from=rfl_builder --chown=fuzzyuser:fuzzyuser /RFL/lib/RFL/include /usr/include/old_RFL/
COPY --from=rfl_builder --chown=fuzzyuser:fuzzyuser /RFL/lib/RFL/bin/ /usr/lib/old_RFL

COPY --chown=fuzzyuser:fuzzyuser . /home/fuzzyuser/RFL

# Set up environment for C++ development
ENV CXX=g++
ENV CC=gcc
ENV PATH="/usr/lib/RFL:/usr/lib/old_RFL:${PATH}"
