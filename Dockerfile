FROM ubuntu:24.04 AS build
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git cmake make gcc-12 g++-12 pkg-config \
    libyaml-cpp-dev libomp-dev librange-v3-dev libseqan3-dev \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/RY-SlopeSearch
COPY . .

RUN git submodule update --init --recursive || true

RUN mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_C_COMPILER=/usr/bin/gcc-12 \
          -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
          .. && \
    make -j

FROM ubuntu:24.04 AS runtime
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libyaml-cpp0.8 libgomp1 ca-certificates \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/RY-SlopeSearch /opt/RY-SlopeSearch
RUN cp /opt/RY-SlopeSearch/build/RY-SlopeSearch /usr/local/bin/RY-SlopeSearch

RUN mkdir -p /out
WORKDIR /work
ENTRYPOINT ["RY-SlopeSearch"]
CMD ["--help"]
