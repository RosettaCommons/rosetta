# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

FROM ubuntu:18.04


# prevent any user interaction during install
ENV DEBIAN_FRONTEND=noninteractive

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN apt-get -y update
RUN apt-get -y install --no-install-recommends mc git curl
RUN apt-get -y install --no-install-recommends cmake ninja-build
RUN apt-get -y install --no-install-recommends build-essential clang
RUN apt-get -y install --no-install-recommends zlib1g-dev libssl-dev libbz2-dev libreadline-dev libsqlite3-dev liblzma-dev libffi-dev

RUN apt-get -y install --no-install-recommends python2.7 python-dev
RUN apt-get -y install --no-install-recommends python3 python3-dev python3-numpy python3-setuptools python3-distutils python3-openssl python3-venv

RUN apt-get -y install --no-install-recommends clang-tidy
RUN apt-get -y install --no-install-recommends ca-certificates

RUN apt-get -y install --no-install-recommends

#RUN apt-get clean && rm -rf /var/lib/apt/lists/*
