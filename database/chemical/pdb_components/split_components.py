#!/usr/bin/env python

# components.cif is 0-9, A-Z, so we don't need to store in a dict temporarily

line5 = "0"
outf = open("components.0.cif", "w")

with open("components.cif") as f:
    for line in f:
        if line[0:5] == "data_":
            if line[5] != line5:
                outf.close()
                line5 = line[5]
                outf = open("components.{}.cif".format(line5), "w")
        outf.write(line)
