#!/usr/bin/python

import re

letters = "ABCFGHJKXYZ"
for line in open("raw_clusters") :
    data = re.split(" +",line.strip())
    print "%4d %4d %4d :" % (int(data[0]),int(data[1]),int(data[2])),
    members = sorted(data[4:])
    for i in members :
	print letters[int(i)-1],
    print


