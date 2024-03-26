# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov


FROM ubuntu:22.04

# prevent any user interaction during install
ENV DEBIAN_FRONTEND=noninteractive

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN apt-get  -y update
RUN apt-get -y install --no-install-recommends mc git curl wget unzip
RUN apt-get -y install --no-install-recommends cmake ninja-build
RUN apt-get -y install --no-install-recommends build-essential clang
RUN apt-get -y install --no-install-recommends zlib1g-dev libssl-dev libbz2-dev libreadline-dev libsqlite3-dev liblzma-dev libffi-dev

RUN apt-get -y install --no-install-recommends python3 python3-dev python3-numpy python3-setuptools python3-distutils python3-openssl python3-venv

RUN apt-get -y install --no-install-recommends ca-certificates

RUN apt-get -y install --no-install-recommends python-dev-is-python3


# Download and extract TensorFlow
RUN wget https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-1.15.0.tar.gz
RUN tar -xzf libtensorflow-cpu-linux-x86_64-1.15.0.tar.gz -C /usr/local && rm libtensorflow-cpu-linux-x86_64-1.15.0.tar.gz

# Set environment variables for TensorFlow
ENV LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

# Download and extract PyTorch
RUN cd /root && wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -O libtorch.zip
RUN cd /root && unzip libtorch.zip && rm libtorch.zip
RUN cp -r /root/libtorch/* /usr/local && rm -rf /root/libtorch

# setting up python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

#RUN apt-get clean && rm -rf /var/lib/apt/lists/*
