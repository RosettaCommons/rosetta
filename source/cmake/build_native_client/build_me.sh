#!/bin/bash

export CC=/Users/mtyka/src/nacl_sdk/pepper_23/toolchain/mac_x86_newlib/bin/x86_64-nacl-gcc
export CXX=/Users/mtyka/src/nacl_sdk/pepper_23/toolchain/mac_x86_newlib/bin/x86_64-nacl-g++
cmake ./
make -j 4


