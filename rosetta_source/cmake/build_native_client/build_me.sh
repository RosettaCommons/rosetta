#!/bin/bash

export CC=../../external/native_client_sdk_0_2_803_0/toolchain/mac_x86/bin/nacl-gcc
export CXX=../../external/native_client_sdk_0_2_803_0/toolchain/mac_x86/bin/nacl-g++
cmake ./
make -j 4


