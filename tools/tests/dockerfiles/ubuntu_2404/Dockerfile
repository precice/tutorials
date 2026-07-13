FROM ubuntu:24.04 AS base_image
# We set a senseful value, but still have the possibility to influence this via the build time arguments. 
# When the dockerfile is built using the systemtests.py we set the PRECICE_UID and PRECICE_GID to the user executing the system tests.
# This ensures no file ownership problems down the line and is the easiest fix, as we normally build the containers locally
# If not built via the systemtests.py, it is possible to specify it manually, but 1000 would be the default anyway.
ARG PRECICE_UID=1000
ARG PRECICE_GID=1000

USER root
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive

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
USER precice
### end of base_image stage ###


FROM base_image AS precice_dependencies
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
### end of precice_dependencies stage ###


FROM precice_dependencies AS precice
# Build & install precice into /home/precice/precice
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
    make all install -j $(nproc) && \
    cd ../.. && rm -rf precice

ENV PATH="${PATH}:/home/precice/.local/bin" 
ENV LD_LIBRARY_PATH="/home/precice/.local/lib" 
ENV CPATH="/home/precice/.local/include"
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig"
ENV CMAKE_PREFIX_PATH="/home/precice/.local"
### end of precice stage ###


FROM precice AS python_bindings
ARG PYTHON_BINDINGS_REF

USER precice
WORKDIR /home/precice
# Builds the precice python bindings for python3
# Installs also matplotlib as it is needed for the elastic-tube 1d fluid-python participant.
RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip3 install git+https://github.com/precice/python-bindings.git@${PYTHON_BINDINGS_REF} && \
    pip3 install matplotlib 
### end of python_bindings stage ###


FROM precice AS micro_manager
ARG MICRO_MANAGER_REF

USER precice
WORKDIR /home/precice
RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    git clone https://github.com/precice/micro-manager.git && \
    cd micro-manager && \
    if [ -n "${MICRO_MANAGER_PR}" ]; then git fetch origin pull/${MICRO_MANAGER_PR}/head; fi && \
    git checkout ${MICRO_MANAGER_REF} && \
    pip3 install . && \
    cd .. && \
    pip3 install pybind11
### end of micro_manager stage ###


#####################################################################
#   Non-dependent adapter build stages - Order alphabetically       #
#####################################################################

FROM precice AS calculix_adapter
ARG CALCULIX_VERSION
ARG CALCULIX_ADAPTER_PR
ARG CALCULIX_ADAPTER_REF

USER root
RUN apt-get -qq update && \
    apt-get -qq install libarpack2-dev libspooles-dev libyaml-cpp-dev
USER precice

WORKDIR /home/precice
RUN wget -nv http://www.dhondt.de/ccx_${CALCULIX_VERSION}.src.tar.bz2 && \
    tar xvjf ccx_${CALCULIX_VERSION}.src.tar.bz2 && \
    rm -fv ccx_${CALCULIX_VERSION}.src.tar.bz2

WORKDIR /home/precice
RUN git clone https://github.com/precice/calculix-adapter.git && \
    cd calculix-adapter && \
    if [ -n "${CALCULIX_ADAPTER_PR}" ]; then git fetch origin pull/${CALCULIX_ADAPTER_PR}/head; fi && \
    git checkout ${CALCULIX_ADAPTER_REF} &&\
    make CXX_VERSION=${CALCULIX_VERSION} ADDITIONAL_FFLAGS="-fallow-argument-mismatch" -j $(nproc) && \
    ln -s /home/precice/calculix-adapter/bin/ccx_preCICE /home/precice/.local/bin/ccx_preCICE
### end of calculix_adapter stage ###


FROM precice AS dealii_adapter
ARG DEALII_ADAPTER_PR
ARG DEALII_ADAPTER_REF

USER root
RUN apt-get update &&\
    apt-get -qq install libdeal.ii-dev libdeal.ii-doc cmake make g++

USER precice
WORKDIR /home/precice
RUN git clone https://github.com/precice/dealii-adapter.git &&\
    cd dealii-adapter && \
    if [ -n "${DEALII_ADAPTER_PR}" ]; then git fetch origin pull/${DEALII_ADAPTER_PR}/head; fi && \
    git checkout ${DEALII_ADAPTER_REF} && \
    cmake . && \
    make -j $(nproc)

ENV PATH="/home/precice/dealii-adapter/:$PATH"
### end of dealii_adapter stage ###


FROM micro_manager AS dumux_adapter
# The micro-manager is needed for some dumux-adapter cases
ARG DUNE_VERSION_DUMUX
ARG DUMUX_VERSION
ARG DUMUX_ADAPTER_PR
ARG DUMUX_ADAPTER_REF

USER precice
WORKDIR /home/precice/dumux

# This cascade of git clones often led to HTTP 502. Adding some delay and retries in between as a workaround.
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-common.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
        rm -rf dune-common; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-geometry.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
        rm -rf dune-geometry; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-grid.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
        rm -rf dune-grid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-istl.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
        rm -rf dune-istl; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-subgrid.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
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
        if git clone --depth 1 https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield.git -b cell_problems; then break; fi; \
        rm -rf dumux-phasefield; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-SPGrid.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
        rm -rf dune-SPGrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-localfunctions.git -b releases/${DUNE_VERSION_DUMUX}; then break; fi; \
        rm -rf dune-localfunctions; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN ./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts cmake &&\
    ./dune-common/bin/dunecontrol make -j $(nproc)

RUN git clone https://github.com/precice/dumux-adapter.git &&\
    cd dumux-adapter && \
    if [ -n "${DUMUX_ADAPTER_PR}" ]; then git fetch origin pull/${DUMUX_ADAPTER_PR}/head; fi && \
    git checkout ${DUMUX_ADAPTER_REF} && \
    cd .. &&\
    ./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts cmake &&\
    ./dune-common/bin/dunecontrol make -j $(nproc)

ENV PATH="/home/precice/dumux/dune-common/bin:${PATH}"
### end of dumux_adapter stage ###


FROM precice AS dune_adapter
ARG DUNE_VERSION_DUNE
ARG DUNE_ADAPTER_PR
ARG DUNE_ADAPTER_REF

USER precice
WORKDIR /home/precice/dune
ENV PATH="/home/precice/dune/dune-common/bin:${PATH}"

# This cascade of git clones often led to HTTP 502. Adding some delay and retries in between as a workaround.
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-common.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-common; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-geometry.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-geometry; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-grid.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-grid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-istl.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-istl; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/core/dune-localfunctions.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-localfunctions; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-foamgrid.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-foamgrid; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/staging/dune-functions.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-functions; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/staging/dune-typetree.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
        rm -rf dune-typetree; \
        if [[ $i -eq 3 ]]; then exit 1; fi; \
        sleep 10; \
    done
RUN for i in $(seq 1 3); do \
        if git clone --depth 1 https://gitlab.dune-project.org/staging/dune-uggrid.git -b releases/${DUNE_VERSION_DUNE}; then break; fi; \
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
RUN git clone https://github.com/precice/dune-adapter.git&&\
    cd dune-adapter && \
    if [ -n "${DUNE_ADAPTER_PR}" ]; then git fetch origin pull/${DUNE_ADAPTER_PR}/head; fi && \
    git checkout ${DUNE_ADAPTER_REF} && \
    cd .. &&\
    ./dune-common/bin/dunecontrol all

# Make dune-adapter/dune-precice-howto/build-cmake/examples/dune-perpendicular-flap discoverable
ENV PATH="/home/precice/dune/dune-adapter/dune-precice-howto/build-cmake/examples/:${PATH}"
### end of dune_adapter stage ###


FROM precice AS fenics_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
ARG FENICS_ADAPTER_REF

USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenics

USER precice
RUN python3 -m venv --system-site-packages /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip3 install git+https://github.com/precice/fenics-adapter.git@${FENICS_ADAPTER_REF}
### end of fenics_adapter stage ###


FROM precice AS fenicsx_adapter
COPY --from=python_bindings /home/precice/.local /home/precice/.local
ARG FENICSX_ADAPTER_REF

USER root
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qq install --no-install-recommends fenicsx

USER precice
RUN python3 -m venv --system-site-packages /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip3 install git+https://github.com/precice/fenicsx-adapter.git@${FENICSX_ADAPTER_REF}
### end of fenicsx_adapter stage ###


FROM precice AS mercurydpm_adapter

USER precice
WORKDIR /home/precice
# The channel-transport-particles tutorial is hosted in a frozen fork of the MercuryDPM repository.
# That's why there are no arguments for version, ref, or pr.
RUN git clone https://bitbucket.org/davidscn/mercurydpm.git &&\
    cd mercurydpm && \
    git checkout channel-transport-tutorial && \
    mkdir -p build && cd build && \
    cmake -D MercuryDPM_PreCICE_COUPLING="ON" .. && \
    make -j $(nproc) ChannelTransport

ENV PATH="/home/precice/mercurydpm/build/Drivers/PreCICE:$PATH"
### end of mercurydpm_adapter stage ###


FROM precice AS openfoam_adapter
ARG OPENFOAM_EXECUTABLE
ARG OPENFOAM_ADAPTER_PR
ARG OPENFOAM_ADAPTER_REF

USER root
RUN apt-get update &&\
    wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash &&\
    apt-get -qq install ${OPENFOAM_EXECUTABLE}-dev &&\
    ln -s $(which ${OPENFOAM_EXECUTABLE} ) /usr/bin/openfoam

USER precice
WORKDIR /home/precice
RUN git clone https://github.com/precice/openfoam-adapter.git &&\
    cd openfoam-adapter && \
    if [ -n "${OPENFOAM_ADAPTER_PR}" ]; then git fetch origin pull/${OPENFOAM_ADAPTER_PR}/head; fi && \
    git checkout ${OPENFOAM_ADAPTER_REF} && \
    /usr/bin/${OPENFOAM_EXECUTABLE} ./Allwmake -j $(nproc)
### end of openfoam_adapter stage ###


FROM python_bindings AS su2_adapter
ARG SU2_VERSION
ARG SU2_ADAPTER_PR
ARG SU2_ADAPTER_REF

USER root
RUN apt-get -qq update && \
    apt-get -qq install swig python3-mpi4py python3-setuptools

USER precice
WORKDIR /home/precice
# Download and build SU2 (We could also use pre-built binaries from the SU2 releases)
# The sed command applies a patch needed for Ubuntu 24.04.
RUN wget -nv https://github.com/su2code/SU2/archive/refs/tags/v${SU2_VERSION}.tar.gz && \
    tar xvzf v${SU2_VERSION}.tar.gz && \
    rm -fv v${SU2_VERSION}.tar.gz
RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate

# Some of these variables are used by the scripts that follow
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
### end of su2_adapter stage ###


#####################################################################
#            Build stages using third-party base images             #
#####################################################################


FROM ubuntu:26.04 AS aste

ARG PRECICE_UID=1000
ARG PRECICE_GID=1000
ARG PRECICE_PR
ARG PRECICE_REF
ARG PRECICE_PRESET
ARG ASTE_PR
ARG ASTE_REF

USER root
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive

# Set the same user as the rest of the system tests components, for consistency.
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

# Additional dependencies for ASTE
RUN apt-get -qq update && \
    apt-get -qq -y install \
        libmetis-dev \
        libvtk9-dev

USER precice
WORKDIR /home/precice

# Build & install precice into /home/precice/precice
RUN git clone https://github.com/precice/precice.git precice && \
    cd precice && \
    if [ -n "${PRECICE_PR}" ]; then git fetch origin pull/${PRECICE_PR}/head; fi && \
    git checkout ${PRECICE_REF} && \
    mkdir build && cd build &&\
    cmake .. --preset=${PRECICE_PRESET} -DCMAKE_INSTALL_PREFIX=/home/precice/.local/ -DBUILD_TESTING=OFF && \
    make all install -j $(nproc) && \
    cd ../.. && rm -rf precice

ENV PATH="${PATH}:/home/precice/.local/bin" 
ENV LD_LIBRARY_PATH="/home/precice/.local/lib"
ENV CPATH="/home/precice/.local/include"
ENV PKG_CONFIG_PATH="/home/precice/.local/lib/pkgconfig"
ENV CMAKE_PREFIX_PATH="/home/precice/.local"

RUN git clone https://github.com/precice/aste.git aste && \
    cd aste && \
    if [ -n "${ASTE_PR}" ]; then git fetch origin pull/${ASTE_PR}/head; fi && \
    git checkout ${ASTE_REF}

RUN python3 -m venv /home/precice/venv && \
    . /home/precice/venv/bin/activate && \
    pip install -r /home/precice/aste/requirements.txt && \
    pip install vtk

RUN cd /home/precice/aste && \
    mkdir build && cd build &&\
    cmake .. -DCMAKE_INSTALL_PREFIX=/home/precice/.local/ && \
    make all install -j $(nproc) && \
    cd ../.. && rm -rf aste

WORKDIR /home/precice