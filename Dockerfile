FROM RFL/docker/cmake AS builder

RUN apt-get install -y          \
    libarmadillo-dev \
    libgsl-dev          \
    git \
    doxygen             \
    graphviz            \
    libhdf5-serial-dev

COPY ./  /RFL/

WORKDIR /RFL

RUN cmake -B ./build . \
    && cmake --build ./build --target all -j 4

RUN ( cd build && ctest -j 2 )\
    && cmake --install ./build

FROM RFL/docker/cmake

COPY --from=builder /RFL/lib/new_RFL/include /usr/include/RFL/
COPY --from=builder /RFL/lib/new_RFL/bin/ /usr/lib/RFL

COPY --from=builder /RFL/lib/RFL/include /usr/include/old_RFL/
COPY --from=builder /RFL/lib/RFL/bin/ /usr/lib/old_RFL
