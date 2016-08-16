#!/bin/bash 

## (c) Copyright Rosetta Commons Member Institutions.
## (c) This file is part of the Rosetta software suite and is made available under license.
## (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
## (c) For more information, see http://www.rosettacommons.org. Questions about this can be
## (c) addressed to University of Washington CoMotion, email: license@uw.edu.
##
# Author:  Mike Tyka

if [[ $# -eq 0 ]]; then
echo "Usage:  cat Myfile | column 1 4 5 "
echo "   will filter out columns 1 4 and 5 from the stream produced by 'cat Myfile' "
fi
if [[ $# -eq 1 ]]; then cat /dev/stdin | awk  "{ print \$$1 }"; fi
if [[ $# -eq 2 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2 }"; fi
if [[ $# -eq 3 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3 }"; fi
if [[ $# -eq 4 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3,\$$4 }"; fi
if [[ $# -eq 5 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3,\$$4,\$$5 }"; fi
if [[ $# -eq 6 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3,\$$4,\$$5,\$$6 }"; fi
if [[ $# -eq 7 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3,\$$4,\$$5,\$$6,\$$7 }"; fi
if [[ $# -eq 8 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3,\$$4,\$$5,\$$6,\$$7,\$$8 }"; fi
if [[ $# -eq 9 ]]; then cat /dev/stdin | awk  "{ print \$$1,\$$2,\$$3,\$$4,\$$5,\$$6,\$$7,\$$8,\$$9 }"; fi
