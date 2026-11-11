#!/usr/bin/env python

stem = None
outf = None

def get_stem(line):
    ccd_id = line.split("_")[1].strip()
    if len(ccd_id) <= 3:
        return ccd_id[0]
    elif len(ccd_id) <= 5 :
        return ccd_id[:-2]
    else:
        raise ValueError(f"ISSUE WITH CCD IDs -- they shouldn't be more than five characters yet. (Got {ccd_id}) Likely formatting issue in file!")

with open("components.cif") as f:
    for line in f:
        if line[0:5] == "data_":
            if get_stem(line) != stem:
                if outf is not None:
                    outf.close()
                stem = get_stem(line)
                outf = open("components.{}.cif".format(stem), "a")
        outf.write(line)
