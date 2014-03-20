#!/bin/bash

export CC=/Users/christoffernorn/nacl_sdk/pepper_33/toolchain/mac_x86_newlib/bin/x86_64-nacl-gcc
export CXX=/Users/christoffernorn/nacl_sdk/pepper_33/toolchain/mac_x86_newlib/bin/x86_64-nacl-g++
cmake ./
make VERBOSE=1


