#!/bin/bash

# First setup nacl_sdk
# Find instructions for that here: https://developer.chrome.com/native-client/sdk/download

make clean

cd ../
python make_project.py all
cd build_native_client 

cd ../../../source/external/zlib-1.2.8/
~/nacl_sdk/pepper_33/toolchain/mac_x86_newlib/bin/x86_64-nacl-gcc -c *.c
~/nacl_sdk/pepper_33/toolchain/mac_x86_newlib/bin/x86_64-nacl-ar -rsc libz-static.a *.o
mv libz-static.a ../../cmake/build_native_client/
cd ../../cmake/build_native_client

export CC=/Users/christoffernorn/nacl_sdk/pepper_33/toolchain/mac_x86_newlib/bin/x86_64-nacl-gcc
export CXX=/Users/christoffernorn/nacl_sdk/pepper_33/toolchain/mac_x86_newlib/bin/x86_64-nacl-g++
cmake ./
make VERBOSE=1


