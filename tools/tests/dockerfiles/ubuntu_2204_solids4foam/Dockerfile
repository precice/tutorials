# A Dockerfile specifically for testing solids4foam,
# based on the ubuntu_2404/Dockerfile.

FROM solids4foam/solids4foam:v2.3-openfoam-v2412-latest as solids4foam

ARG PRECICE_UID=1000
ARG PRECICE_GID=1000
ARG OPENFOAM_ADAPTER_PR
ARG OPENFOAM_ADAPTER_REF

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

# The base image provides all binaries under /root/, but we prefer to build as a regular user.
RUN mv /root/OpenFOAM/ /home/precice/OpenFOAM/ && \
    chown -R "${PRECICE_UID}:${PRECICE_GID}" /home/precice/OpenFOAM

# Set a general openfoam executable, to make the component template version-independent
RUN ln -sf /usr/bin/openfoam2412 /usr/bin/openfoam

RUN apt-get -qq update && \
    apt-get -qq -y install wget git

# preCICE v3.4 was the last release to support Ubuntu 22.04 (used as a base image),
# so there is no reason to accept the preCICE version as an argument and build it from source.
# Rework this part once a newer base image is used.
RUN  wget https://github.com/precice/precice/releases/download/v3.4.1/libprecice3_3.4.1_jammy.deb && \
     apt-get -qq -y install ./libprecice3_3.4.1_jammy.deb

USER precice
WORKDIR /home/precice

RUN git clone https://github.com/precice/openfoam-adapter.git &&\
    cd openfoam-adapter && \
    if [ -n "${OPENFOAM_ADAPTER_PR}" ]; then git fetch origin pull/${OPENFOAM_ADAPTER_PR}/head; fi && \
    git checkout ${OPENFOAM_ADAPTER_REF} && \
    /usr/bin/openfoam2412 ./Allwmake -j $(nproc)
