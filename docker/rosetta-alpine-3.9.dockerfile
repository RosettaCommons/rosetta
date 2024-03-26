# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

FROM alpine:3.9

RUN apk update
RUN apk add mc python git bash build-base zlib-dev
RUN apk add python3 cmake ninja
RUN apk add libexecinfo-dev

RUN apk add clang clang-dev

ENTRYPOINT cd ~ && /bin/bash
