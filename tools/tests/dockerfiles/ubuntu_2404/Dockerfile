FROM ubuntu:24.04 AS base_image
USER root
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
# We set a sensical value, but still have the possibilty to influence this via the build time arguments. 
# When the dockerfile is built using the systemtests.py we set the PRECICE_UID and PRECICE_GID to the user executing the systemtests.
# This ensures no file ownership problems down the line and is the most easy fix, as we normally built the containers locally
# If not built via the systemtests.py its either possible to specify manually but 1000 would be the default anyway.
ARG PRECICE_UID=1000
ARG PRECICE_GID=1000
RUN set -eux; \
    if getent group "${PRECICE_GID}" >/dev/null; then \
        existing_group="$(getent group "${PRECICE_GID}" | cut -d: -f1)"; \
        if [ "${existing_group}" != "precice" ]; then \
            groupmod -n precice "${existing_group}"; \
        fi; \
    else \
        groupadd -g "${PRECICE_GID}" precice; \
    fi; \
    if getent passwd "${PRECICE_UID}" >/dev/null; then \
        existing_user="$(getent passwd "${PRECICE_UID}" | cut -d: -f1)"; \
        if [ "${existing_user}" != "precice" ]; then \
            usermod -l precice -d /home/precice -m "${existing_user}"; \
        fi; \
        usermod -g "${PRECICE_GID}" precice; \
    else \
        useradd -u "${PRECICE_UID}" -g "${PRECICE_GID}" -ms /bin/bash precice; \
    fi; \
    mkdir -p /home/precice; \
    chown -R "${PRECICE_UID}:${PRECICE_GID}" /home/precice
ENV PATH="${PATH}:/home/precice/.local/bin" 
ENV LD_LIBRARY_PATH="/home/precice/.local/lib:${LD_LIBRARY_PATH}" 
ENV CPATH="/home/precice/.local/include:$CPATH"
# Enable detection with pkg-config and CMake
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig:$PKG_CONFIG_PATH"
ENV CMAKE_PREFIX_PATH="/home/precice/.local:$CMAKE_PREFIX_PATH"
USER precice

FROM base_image AS precice_dependecies
USER root
# Dependencies for preCICE and common dependencies for the tutorials/tests
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
        libopenmpi-dev \
        libxml2-dev \
        lsb-release \
        openmpi-bin \
        petsc-dev \
        python3-dev \
        python3-numpy \
        python3-pip \
        python3-venv \
        pkg-config \
        wget \
        inotify-tools
# The following are dependencies of gmsh, needed by some tutorials
RUN apt-get -qq update && \
    apt-get -qq install \
    libglu1-mesa \
    libxrender1 \
    libxcursor1 \
    libxft2 \
    libxinerama1
USER precice

FROM precice_dependecies AS precice
# Install & build precice into /home/precice/precice
ARG PRECICE_PR
ARG PRECICE_REF
ARG PRECICE_PRESET
USER precice
WORKDIR /home/precice
RUN git clone https://github.com/precice/precice.git precice && \
    cd precice && \
    if [ -n "${PRECICE_PR}" ]; then git fetch origin pull/${PRECICE_PR}/head; fi && \
    git checkout ${PRECICE_REF} && \
    mkdir build && cd build &&\
    cmake .. --preset=${PRECICE_PRESET} -DCMAKE_INSTALL_PREFIX=/home/precice/.local/ -DBUILD_TESTING=OFF && \
    make all install -j $(nproc)

FROM precice_dependecies AS openfoam_adapter
ARG OPENFOAM_EXECUTABLE
USER root
RUN apt-get update &&\
    wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash &&\
    apt-get -qq install ${OPENFOAM_EXECUTABLE}-dev &&\
    ln -s $(which ${OPENFOAM_EXECUTABLE} ) /usr/bin/openfoam
USER precice
COPY --from=precice /home/precice/.local/ /home/precice/.local/
ARG OPENFOAM_ADAPTER_PR
ARG OPENFOAM_ADAPTER_REF
# Build the OpenFOAM adapter 
USER precice
WORKDIR /home/precice
RUN git clone https://github.com/precice/openfoam-adapter.git &&\
    cd openfoam-adapter && \
    if [ -n "${OPENFOAM_ADAPTER_PR}" ]; then git fetch origin pull/${OPENFOAM_ADAPTER_PR}/head; fi && \
    git checkout ${OPENFOAM_ADAPTER_REF} && \
    /usr/bin/${OPENFOAM_EXECUTABLE} ./Allwmake -j $(nproc)


FROM precice_dependecies AS python_bindings
COPY --from=precice /home/precice/.local/ /home/precice/.local/
ARG PYTHON_BINDINGS_REF
USER precice
WORKDIR /home/precice
# Builds the precice python bindings for python3
# Installs also matplotlib as its needed for the elastic-tube 1d fluid-python participant.
RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip3 install git+https://github.com/precice/python-bindings.git@${PYTHON_BINDINGS_REF} && \
    pip3 install matplotlib 

FROM precice_dependecies AS fenics_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenics
USER precice
ARG FENICS_ADAPTER_REF
# Building fenics-adapter
RUN python3 -m venv --system-site-packages /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip3 install git+https://github.com/precice/fenics-adapter.git@${FENICS_ADAPTER_REF}

FROM precice_dependecies AS fenicsx_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenicsx
USER precice
ARG FENICSX_ADAPTER_REF
# Building fenicsx-adapter
RUN python3 -m venv --system-site-packages /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip3 install git+https://github.com/precice/fenicsx-adapter.git@${FENICSX_ADAPTER_REF}


FROM precice_dependecies AS calculix_adapter
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

ARG CALCULIX_ADAPTER_PR
ARG CALCULIX_ADAPTER_REF
WORKDIR /home/precice
RUN git clone https://github.com/precice/calculix-adapter.git && \
    cd calculix-adapter && \
    if [ -n "${CALCULIX_ADAPTER_PR}" ]; then git fetch origin pull/${CALCULIX_ADAPTER_PR}/head; fi && \
    git checkout ${CALCULIX_ADAPTER_REF} &&\
    make CXX_VERSION=${CALCULIX_VERSION} ADDITIONAL_FFLAGS="-fallow-argument-mismatch" -j $(nproc) && \
    ln -s /home/precice/calculix-adapter/bin/ccx_preCICE /home/precice/.local/bin/ccx_preCICE

FROM python_bindings AS su2_adapter
COPY --from=precice /home/precice/.local /home/precice/.local
USER root
RUN apt-get -qq update && \
    apt-get -qq install swig python3-mpi4py python3-setuptools
ARG SU2_VERSION
USER precice

# Download and build SU2 (We could also use pre-built binaries from the SU2 releases)
# The sed command applies a patch needed for Ubuntu 24.04.
WORKDIR /home/precice
RUN wget https://github.com/su2code/SU2/archive/refs/tags/v${SU2_VERSION}.tar.gz && \
    tar xvzf v${SU2_VERSION}.tar.gz && \
    rm -fv v${SU2_VERSION}.tar.gz
RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate
ARG SU2_ADAPTER_PR
ARG SU2_ADAPTER_REF
WORKDIR /home/precice
ENV SU2_RUN="/home/precice/SU2_RUN" 
ENV SU2_HOME="/home/precice/SU2-${SU2_VERSION}"
ENV PATH="/home/precice/su2-adapter/run:$SU2_RUN:$PATH"
ENV PYTHONPATH="$SU2_RUN:$PYTHONPATH"
RUN git clone https://github.com/precice/su2-adapter.git && \
    cd su2-adapter &&\
    if [ -n "${SU2_ADAPTER_PR}" ]; then git fetch origin pull/${SU2_ADAPTER_PR}/head; fi && \
    git checkout ${SU2_ADAPTER_REF} &&\
    ./su2AdapterInstall
RUN cd "${SU2_HOME}" &&\
    sed -i '1s/^/#include <cstdint>\n/' SU2_CFD/src/output/filewriter/CParaviewXMLFileWriter.cpp &&\
    ./meson.py build -Denable-pywrapper=true --prefix=$SU2_RUN &&\
    ./ninja -C build install

FROM precice_dependecies AS dealii_adapter
USER root
RUN apt-get update &&\
    apt-get -qq install libdeal.ii-dev libdeal.ii-doc cmake make g++
USER precice
COPY --from=precice /home/precice/.local/ /home/precice/.local/
ARG DEALII_ADAPTER_PR
ARG DEALII_ADAPTER_REF
# Build the deal.II adapter 
USER precice
WORKDIR /home/precice
ENV PATH="/home/precice/dealii-adapter/:$PATH"
RUN git clone https://github.com/precice/dealii-adapter.git &&\
    cd dealii-adapter && \
    if [ -n "${DEALII_ADAPTER_PR}" ]; then git fetch origin pull/${DEALII_ADAPTER_PR}/head; fi && \
    git checkout ${DEALII_ADAPTER_REF} && \
    cmake . && \
    make -j $(nproc)

FROM precice_dependecies AS mercurydpm_adapter
USER precice
COPY --from=precice /home/precice/.local /home/precice/.local
WORKDIR /home/precice
ENV PATH="/home/precice/mercurydpm/build/Drivers/PreCICE:$PATH"
# The channel-transport-particles tutorial is hosted in a frozen fork of the MercuryDPM repository.
# That's why there are no arguments for version, ref, or pr.
RUN git clone https://bitbucket.org/davidscn/mercurydpm.git &&\
    cd mercurydpm && \
    git checkout channel-transport-tutorial && \
    mkdir -p build && cd build && \
    cmake -D MercuryDPM_PreCICE_COUPLING="ON" .. && \
    make -j $(nproc) ChannelTransport

FROM precice_dependecies AS micro_manager
USER precice
WORKDIR /home/precice
COPY --from=precice /home/precice/.local/ /home/precice/.local/

ARG MICRO_MANAGER_REF
RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    git clone https://github.com/precice/micro-manager.git && \
    cd micro-manager && \
    if [ -n "${MICRO_MANAGER_PR}" ]; then git fetch origin pull/${MICRO_MANAGER_PR}/head; fi && \
    git checkout ${MICRO_MANAGER_REF} && \
    pip3 install . && \
    pip3 install pybind11
WORKDIR /home/precice

FROM precice_dependecies AS dune
USER root
COPY --from=precice /home/precice/.local/ /home/precice/.local/
ARG DUNE_VERSION
ARG DUMUX_VERSION
USER precice
WORKDIR /home/precice/dune
ENV PATH="/home/precice/dune/dune-common/bin:${PATH}"
# This cascade of git clones often led to HTTP 502. Adding some delay and retries in between as a workaround.
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-common.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-common; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-geometry.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-geometry; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-grid.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-grid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-istl.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-istl; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-subgrid.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-subgrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git -b releases/${DUMUX_VERSION}; then break; fi; \
        rm -rf dumux; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone -b cell_problems https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield.git; then break; fi; \
        rm -rf dumux-phasefield; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-SPGrid.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-SPGrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-localfunctions.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-localfunctions; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
# Up to here, all the dependencies are intended for dumux.
# DUNE modules needed for the dune adapter:
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-foamgrid.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-foamgrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-foamgrid.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-foamgrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/staging/dune-functions.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-functions; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/staging/dune-typetree.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-typetree; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/staging/dune-uggrid.git -b releases/${DUNE_VERSION}; then break; fi; \
        rm -rf dune-uggrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
# Modules specific to the dune adapter:
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://github.com/maxfirmbach/dune-elastodynamics.git -b master; then break; fi; \
        rm -rf dune-elastodynamics; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://github.com/precice/dune-adapter.git -b main; then break; fi; \
        rm -rf dune-adapter; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
# Make dune-adapter/dune-precice-howto/build-cmake/examples/dune-perpendicular-flap discoverable
ENV PATH="/home/precice/dune/dune-adapter/dune-precice-howto/build-cmake/examples/:${PATH}"
RUN ./dune-common/bin/dunecontrol all

FROM dune AS dumux_adapter
USER root
COPY --from=micro_manager /home/precice/venv /home/precice/venv
RUN . /home/precice/venv/bin/activate
ARG DUMUX_VERSION
USER precice
WORKDIR /home/precice/dune
# build core DUNE, DuMuX and the adapter
ARG DUMUX_ADAPTER_PR
ARG DUMUX_ADAPTER_REF
RUN git clone https://github.com/precice/dumux-adapter.git &&\
    cd dumux-adapter && \
    if [ -n "${DUMUX_ADAPTER_PR}" ]; then git fetch origin pull/${DUMUX_ADAPTER_PR}/head; fi && \
    git checkout ${DUMUX_ADAPTER_REF} && \
    cd .. &&\
    ./dune-common/bin/dunecontrol --only=dumux-precice
