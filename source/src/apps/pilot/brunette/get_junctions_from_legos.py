from argparse import ArgumentParser
import os
from os import system
from os.path import basename, exists, splitext
from sys import exit, stderr, stdout
from itertools import repeat
import multiprocessing as mp
import glob
import copy
import shutil
import math

def collect_junctions(location):
    data = {}
    results = glob.glob("{}/*.jobs".format(location))
    for fn in results:
        fl = open(fn)
        for line in fl:
            descriptions= line.split("|")
            lego =descriptions[0].split(" ")[0]
            junctions =[]
            for description in descriptions[1:-1]:
                pot_junc = description.split(",")[1]
                path = pot_junc.split("/")
                if(path[-2]!="DHRs"):
                    print(path[-2])
                    item.append(path[-1])

                #print(description)


data = collect_junctions("job_files")

