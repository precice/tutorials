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

RUN useradd -ms /bin/bash precice 
USER precice

FROM base_image as precice_dependecies
USER root
# Installing necessary dependecies for preCICE
RUN apt-get -qq update && \
    apt-get -qq -y install \
        build-essential \
        software-properties-common \
        cmake \
        curl \
        g++ \
        gfortran \
        git \
        libbenchmark-dev \
        libboost-all-dev \
        libeigen3-dev \
        libxml2-dev \
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
# Install & build precice into /home/precice/precice
ARG PRECICE_REF
USER precice
WORKDIR /home/precice
RUN git clone --depth 1 https://github.com/precice/precice.git -b ${PRECICE_REF} precice && \
    cd precice && \
    mkdir build && cd build &&\
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/precice/.local/ -DPRECICE_PETScMapping=OFF -DBUILD_TESTING=OFF && \
    make all install -j $(nproc)

FROM precice_dependecies as openfoam_adapter
ARG OPENFOAM_EXECUTABLE
USER root
RUN apt-get update &&\
    wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash &&\
    apt-get -qq install ${OPENFOAM_EXECUTABLE}-dev &&\
    ln -s $(which ${OPENFOAM_EXECUTABLE} ) /usr/bin/openfoam

USER precice
COPY --from=precice /home/precice/.local/ /home/precice/.local/
ENV PATH="${PATH}:/home/precice/.local/bin"
ENV LD_LIBRARY_PATH="/home/precice/.local/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/home/precice/.local/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/home/precice/.local:$CMAKE_PREFIX_PATH"

ARG OPENFOAM_ADAPTER_REF
# Build the OpenFOAM adapter 
USER precice
WORKDIR /home/precice
RUN git clone --depth 1 --branch ${OPENFOAM_ADAPTER_REF} https://github.com/precice/openfoam-adapter.git &&\
    cd openfoam-adapter && /usr/bin/${OPENFOAM_EXECUTABLE} ./Allwmake -j $(nproc)



FROM precice_dependecies as python_bindings
COPY --from=precice /home/precice/.local /home/precice/.local
ENV PATH="${PATH}:/home/precice/.local/bin"
ENV LD_LIBRARY_PATH="/home/precice/.local/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/home/precice/.local/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/home/precice/.local:$CMAKE_PREFIX_PATH"
USER precice
ARG PYTHON_BINDINGS_REF
WORKDIR /home/precice
# Builds the precice python bindings for python3
RUN pip3 install --user git+https://github.com/precice/python-bindings.git@${PYTHON_BINDINGS_REF}


FROM precice_dependecies as fenics_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenics
USER precice
RUN pip3 install fenics-ufl
ARG FENICS_ADAPTER_REF
# Building fenics-adapter
RUN pip3 install --user git+https://github.com/precice/fenics-adapter.git@${FENICS_ADAPTER_REF}


FROM precice_dependecies as nutils_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
USER precice
# Building fenics-adapter
RUN pip3 install --user nutils



##  STILL WIP
FROM precice_dependecies as calculix_adapter
COPY --from=precice /home/precice/.local /home/precice/.local
ENV PATH="${PATH}:/home/precice/.local/bin"
ENV LD_LIBRARY_PATH="/home/precice/.local/lib:${LD_LIBRARY_PATH}"
ENV CPATH="/home/precice/.local/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/home/precice/.local:$CMAKE_PREFIX_PATH"
USER root
RUN apt-get -qq update && \
    apt-get -qq install libarpack2-dev libspooles-dev libyaml-cpp-dev
USER precice
ARG CALULIX_VERSION="2.20"
# Building fenics-adapter
WORKDIR /home/precice
RUN wget http://www.dhondt.de/ccx_${CALULIX_VERSION}.src.tar.bz2 && \
    tar xvjf ccx_${CALULIX_VERSION}.src.tar.bz2 && 

ARG CALULIX_ADAPTER_REF=v2.20.0
WORKDIR /opt/precice
RUN git clone --depth 1 --branch ${CALULIX_ADAPTER_REF} https://github.com/precice/calculix-adapter.git && \
    cd calculix-adapter && \
    make CCX="/home/precice//CalculiX/ccx_${CALULIX_VERSION}/src" CCX_VERSION="${CALULIX_VERSION}" && \
    ln -s /opt/precice/calculix_source/current/src/ccx_preCICE /home/precice/.local/bin/ccx_preCICE
