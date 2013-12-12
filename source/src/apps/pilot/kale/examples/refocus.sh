#!/usr/bin/env sh

if [ $# -ge 1 ]; then
    ln -sf $1 current_example.cc
    touch current_example.cc
fi

ls -l current_example.cc
