FROM ubuntu:22.04 as base_image
USER root
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
# We set a sensical value, but still have the possibilty to influence this via the build time arguments. 
# When the dockerfile is built using the systemtests.py we set the PRECICE_UID and PRECICE_GID to the user executing the systemtests.
# This ensures no file ownership problems down the line and is the most easy fix, as we normally built the containers locally
# If not built via the systemtests.py its either possible to specify manually but 1000 would be the default anyway.
ARG PRECICE_UID=1000
ARG PRECICE_GID=1000
RUN groupadd -g ${PRECICE_GID} precice && useradd -u ${PRECICE_UID} -g ${PRECICE_GID} -ms /bin/bash precice
ENV PATH="${PATH}:/home/precice/.local/bin" 
ENV LD_LIBRARY_PATH="/home/precice/.local/lib:${LD_LIBRARY_PATH}" 
ENV CPATH="/home/precice/.local/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/home/precice/.local:$CMAKE_PREFIX_PATH"
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
ARG PRECICE_PRESET
USER precice
WORKDIR /home/precice
RUN git clone https://github.com/precice/precice.git precice && \
    cd precice && \
    git checkout ${PRECICE_REF} && \
    mkdir build && cd build &&\
    cmake .. --preset=${PRECICE_PRESET} -DCMAKE_INSTALL_PREFIX=/home/precice/.local/ -DPRECICE_PETScMapping=OFF -DBUILD_TESTING=OFF && \
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
ARG OPENFOAM_ADAPTER_REF
# Build the OpenFOAM adapter 
USER precice
WORKDIR /home/precice
RUN git clone https://github.com/precice/openfoam-adapter.git &&\
    cd openfoam-adapter && \
    git checkout ${OPENFOAM_ADAPTER_REF} && \
    /usr/bin/${OPENFOAM_EXECUTABLE} ./Allwmake -j $(nproc)


FROM precice_dependecies as python_bindings
COPY --from=precice /home/precice/.local/ /home/precice/.local/
ARG PYTHON_BINDINGS_REF
USER precice
WORKDIR /home/precice
# Builds the precice python bindings for python3
# Installs also matplotlib as its needed for the elastic-tube 1d fluid-python participant.
RUN pip3 install --user git+https://github.com/precice/python-bindings.git@${PYTHON_BINDINGS_REF} && \
    pip3 install --user matplotlib 

FROM precice_dependecies as fenics_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenics
USER precice
ARG FENICS_ADAPTER_REF
# Building fenics-adapter
RUN pip3 install --user git+https://github.com/precice/fenics-adapter.git@${FENICS_ADAPTER_REF}


FROM precice_dependecies as nutils_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
USER precice
# Installing nutils - There is no adapter
RUN pip3 install --user nutils


FROM precice_dependecies as calculix_adapter
COPY --from=precice /home/precice/.local /home/precice/.local
USER root
RUN apt-get -qq update && \
    apt-get -qq install libarpack2-dev libspooles-dev libyaml-cpp-dev
ARG CALCULIX_VERSION
USER precice
#Download Calculix
WORKDIR /home/precice
RUN wget http://www.dhondt.de/ccx_${CALCULIX_VERSION}.src.tar.bz2 && \
    tar xvjf ccx_${CALCULIX_VERSION}.src.tar.bz2 && \
    rm -fv ccx_${CALCULIX_VERSION}.src.tar.bz2

ARG CALCULIX_ADAPTER_REF
WORKDIR /home/precice
RUN git clone https://github.com/precice/calculix-adapter.git && \
    cd calculix-adapter && \
    git checkout ${CALCULIX_ADAPTER_REF} &&\
    make CXX_VERSION=${CALCULIX_VERSION} ADDITIONAL_FFLAGS="-fallow-argument-mismatch" -j $(nproc) && \
    ln -s /home/precice/calculix-adapter/bin/ccx_preCICE /home/precice/.local/bin/ccx_preCICE

FROM python_bindings as su2_adapter
COPY --from=precice /home/precice/.local /home/precice/.local
USER root
RUN apt-get -qq update && \
    apt-get -qq install swig
ARG SU2_VERSION
USER precice

# Download and build SU2 (We could also use pre-built binaries from the SU2 releases)
WORKDIR /home/precice
RUN wget https://github.com/su2code/SU2/archive/refs/tags/v${SU2_VERSION}.tar.gz && \
    tar xvzf v${SU2_VERSION}.tar.gz && \
    rm -fv v${SU2_VERSION}.tar.gz
RUN pip3 install --user mpi4py
ARG SU2_ADAPTER_REF
WORKDIR /home/precice
ENV SU2_RUN="/home/precice/SU2_RUN" 
ENV SU2_HOME="/home/precice/SU2-${SU2_VERSION}"
ENV PATH="/home/precice/su2-adapter/run:$SU2_RUN:$PATH"
ENV PYTHONPATH="$SU2_RUN:$PYTHONPATH"
RUN git clone https://github.com/precice/su2-adapter.git && \
    cd su2-adapter &&\
    git checkout ${SU2_ADAPTER_REF} &&\
    ./su2AdapterInstall
RUN cd "${SU2_HOME}" &&\
    ./meson.py build -Denable-pywrapper=true --prefix=$SU2_RUN &&\
    ./ninja -C build install
