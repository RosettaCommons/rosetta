#!/usr/bin/env python

from string import strip

def load(file):
    file = open(file, 'r').readlines()
    loaded_file = []
    for line in file:
        loaded_file.append(strip(line))

    return loaded_file

