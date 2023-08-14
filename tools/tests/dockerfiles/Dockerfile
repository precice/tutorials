FROM ubuntu:22.04 as base_image
USER root
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${PATH}:/home/precice/.local/bin:/opt/precice/installation/bin"
ENV LD_LIBRARY_PATH="/opt/precice/installation/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/opt/precice/installation/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/opt/precice/installation/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/opt/precice/installation:$CMAKE_PREFIX_PATH"

RUN useradd -ms /bin/bash precice && \
    mkdir -p /opt/precice/installation && \
    chown -R precice:precice /opt/precice/

FROM base_image as precice_dependecies
USER root
# Installing necessary dependecies for preCICE
RUN apt-get -qq update && \
    apt-get -qq -y install \
        build-essential \
        software-properties-common \
        ccache \
        cmake \
        curl \
        g++ \
        gfortran \
        git \
        lcov \
        libbenchmark-dev \
        libboost-all-dev \
        libeigen3-dev \
        libxml2-dev \
        lintian \
        lsb-release \
        petsc-dev \
        python3-dev \
        python3-numpy \
        python3-pip \
        python3-venv \
        pkg-config \
        wget
USER precice
RUN python3 -m pip install --user --upgrade pip


FROM precice_dependecies as precice
# Install & build precice into /opt/precice/installation
ARG PRECICE_REF=main
USER precice
WORKDIR /opt/precice
RUN git clone https://github.com/precice/precice.git -b ${PRECICE_REF} precice && \
    cd precice && \
    mkdir build && cd build &&\
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/precice/installation -DPRECICE_PETScMapping=OFF -DBUILD_TESTING=OFF && \
    make all install -j $(nproc)

FROM precice_dependecies as openfoam_adapter
COPY --from=precice /opt/precice/installation /opt/precice/installation
ENV PATH="${PATH}:/home/precice/.local/bin:/opt/precice/installation/bin"
ENV LD_LIBRARY_PATH="/opt/precice/installation/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/opt/precice/installation/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/opt/precice/installation/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/opt/precice/installation:$CMAKE_PREFIX_PATH"
ARG OPENFOAM_EXECUTABLE=openfoam2112
USER root
RUN apt-get update &&\
    wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash &&\
    apt-get -qq install ${OPENFOAM_EXECUTABLE}-dev &&\
    ln -s $(which ${OPENFOAM_EXECUTABLE} ) /usr/bin/openfoam
ARG OPENFOAM_ADAPTER_REF=master
# Build the OpenFOAM adapter 
USER precice
WORKDIR /opt/precice
RUN ls /opt/precice/installation
RUN git clone --depth 1 --branch ${OPENFOAM_ADAPTER_REF} https://github.com/precice/openfoam-adapter.git &&\
    cd openfoam-adapter && /usr/bin/${OPENFOAM_EXECUTABLE} ./Allwmake -j $(nproc)


# Will probably not work
FROM precice_dependecies as python_bindings
COPY --from=precice /opt/precice/installation /opt/precice/installation
ENV PATH="${PATH}:/home/precice/.local/bin:/opt/precice/installation/bin"
ENV LD_LIBRARY_PATH="/opt/precice/installation/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/opt/precice/installation/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/opt/precice/installation/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/opt/precice/installation:$CMAKE_PREFIX_PATH"
ENV PYTHONPATH="/opt/precice/installation/python_bindings:$PYTHONPATH"
USER precice
ARG PYTHON_BINDINGS_REF=master
WORKDIR /opt/precice
# Builds the precice python bindings for python3
RUN pip3 install --target=/opt/precice/installation/python_bindings git+https://github.com/precice/python-bindings.git@${PYTHON_BINDINGS_REF}


# Will probably not work
FROM precice_dependecies as fenics_adapter
COPY --from=precice /opt/precice/installation /opt/precice/installation
COPY --from=python_bindings /opt/precice/installation/python_bindings /opt/precice/installation/python_bindings
ENV PATH="${PATH}:/home/precice/.local/bin:/opt/precice/installation/bin"
ENV LD_LIBRARY_PATH="/opt/precice/installation/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/opt/precice/installation/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/opt/precice/installation/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/opt/precice/installation:$CMAKE_PREFIX_PATH"
ENV PYTHONPATH="/opt/precice/installation/python_bindings:/opt/precice/installation/fencis_adapter:$PYTHONPATH"
USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenics
USER precice
ARG FENICS_ADAPTER_REF=master
# Building fenics-adapter
RUN pip3 install --target=/opt/precice/installation/fencis_adapter git+https://github.com/precice/fenics-adapter.git@$FENICS_ADAPTER_REF
