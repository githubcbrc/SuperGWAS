# =============================================================================
# Stage 1: Build
# =============================================================================
FROM debian:stable-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget unzip \
    build-essential g++ \
    zlib1g-dev libopenblas-openmp-dev \
    && rm -rf /var/lib/apt/lists/*

# Eigen (header-only)
RUN wget -q https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip && \
    unzip -q eigen-3.4.0.zip && \
    mv eigen-3.4.0 /usr/include/eigen3 && \
    rm -f eigen-3.4.0.zip

WORKDIR /build
COPY src/ ./

# Common flags
ARG CXX_FLAGS="-O3 -march=x86-64 -mtune=generic -std=c++17"
ARG EIGEN_INC="-I/usr/include/eigen3"

# Compile all 14 binaries
RUN set -ex && \
    g++ ${CXX_FLAGS} sample_kmers.cpp -o sample_kmers -lz && \
    g++ ${CXX_FLAGS} filter_kmers.cpp -o filter_kmers -lz && \
    g++ ${CXX_FLAGS} ${EIGEN_INC} compute_structure.cpp -o compute_structure && \
    g++ ${CXX_FLAGS} convert_matrix.cpp -o convert_matrix -lz -fopenmp && \
    g++ ${CXX_FLAGS} -fopenmp -DEIGEN_USE_BLAS ${EIGEN_INC} kgwas.cpp -o kgwas -lopenblas -lz && \
    g++ ${CXX_FLAGS} -fopenmp build_index.cpp -o build_index && \
    g++ ${CXX_FLAGS} -fopenmp map_to_ref.cpp -o map_to_ref && \
    g++ ${CXX_FLAGS} -fopenmp map_merge_filter.cpp -o map_merge_filter && \
    g++ ${CXX_FLAGS} ${EIGEN_INC} compute_structure_all.cpp -o compute_structure_all && \
    g++ ${CXX_FLAGS} -fopenmp -DEIGEN_USE_BLAS ${EIGEN_INC} kgwas_2d.cpp -o kgwas_2d -lopenblas -lz && \
    g++ ${CXX_FLAGS} -fopenmp map_2d_pairs.cpp -o map_2d_pairs && \
    g++ ${CXX_FLAGS} ${EIGEN_INC} compute_structure_2w.cpp -o compute_structure_2w && \
    g++ ${CXX_FLAGS} ${EIGEN_INC} compute_projections_2w.cpp -o compute_projections_2w && \
    g++ ${CXX_FLAGS} -fopenmp ${EIGEN_INC} kgwas_2w.cpp -o kgwas_2w

# =============================================================================
# Stage 2: Runtime (no source code, no build tools)
# =============================================================================
FROM debian:stable-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g libgomp1 libopenblas0-openmp \
    python3 python3-numpy python3-matplotlib \
    && rm -rf /var/lib/apt/lists/*

# Copy all binaries from builder
COPY --from=builder \
    /build/sample_kmers \
    /build/filter_kmers \
    /build/compute_structure \
    /build/convert_matrix \
    /build/kgwas \
    /build/build_index \
    /build/map_to_ref \
    /build/map_merge_filter \
    /build/compute_structure_all \
    /build/kgwas_2d \
    /build/map_2d_pairs \
    /build/compute_structure_2w \
    /build/compute_projections_2w \
    /build/kgwas_2w \
    /usr/local/bin/

ENV OMP_PROC_BIND=spread \
    OMP_PLACES=cores \
    OPENBLAS_NUM_THREADS=40 \
    OMP_NUM_THREADS=40

CMD ["bash", "-c", "echo 'Available: sample_kmers, filter_kmers, compute_structure, convert_matrix, kgwas, build_index, map_to_ref, map_merge_filter, compute_structure_all, kgwas_2d, map_2d_pairs, compute_structure_2w, compute_projections_2w, kgwas_2w'"]
