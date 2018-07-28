FROM ubuntu:xenial

# Set an encoding to make things work smoothly.
ENV LANG en_US.UTF-8

# Install basic requirements.
RUN apt-get -qq update && \
    apt-get install -y \
                   file \
                   vim-tiny \
                   bzip2 \
                   make \
                   patch \
                   sudo \
                   tar \
                   openssl \
                   openssh-client \
                   curl \
                   rsync \
                   util-linux \
                   wget && \
    apt-get -qq -y autoremove && \
    apt-get autoclean

# Install system compilers for rosetta builds
RUN apt-get -qq update \
    && apt-get install -y \
                   gcc \
                   g++ \
    && apt-get -qq -y autoremove \
    && apt-get autoclean

# Run common commands
COPY scripts/run_commands /opt/docker/bin/run_commands
RUN /opt/docker/bin/run_commands

# Add a file for users to source to activate the `conda`
# environment `root` and the devtoolset compiler. Also
# add a file that wraps that for use with the `ENTRYPOINT`.
COPY linux-anvil/entrypoint_source /opt/docker/bin/entrypoint_source
COPY scripts/entrypoint /opt/docker/bin/entrypoint

# Ensure that all containers start with tini and the user selected process.
# Activate the `conda` environment `root` and the devtoolset compiler.
# Provide a default command (`bash`), which will start if the user doesn't specify one.
ENTRYPOINT [ "/opt/conda/bin/tini", "--", "/opt/docker/bin/entrypoint" ]
CMD [ "/bin/bash" ]
