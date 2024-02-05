FROM alpine:3 AS build

RUN apk update && \
    apk add --no-cache \
    build-base \
    clang \
    git \
    python3 \
    gtest-dev \
    gmock \
    gcc \
    cmake                 \
    armadillo-dev \
    wget                  \
    gsl-dev            \
    doxygen               \
    graphviz              \
    hdf5-dev


# Build RFL library
WORKDIR /RFL
COPY . /RFL/
ARG BUILD_TYPE=Release
RUN cmake -B ./build -DCMAKE_BUILD_TYPE=${BUILD_TYPE} && \
    cmake --build ./build --config ${BUILD_TYPE} --target all -j 8

##############################################################
FROM alpine:3
RUN apk update && \
    apk add --no-cache \
    git \
    gcc \
    gsl \
    gtest-dev \
    gmock \
    clang \
    libstdc++ \
    armadillo \
    hdf5


RUN addgroup -S shs && adduser -S shs -G shs
USER shs

COPY --chown=shs:shs --from=build \
    /RFL/build \
    /RFL/build


ENTRYPOINT ["/RFL/build/RFL_source/new_source/tests/new_RFLTests"]
