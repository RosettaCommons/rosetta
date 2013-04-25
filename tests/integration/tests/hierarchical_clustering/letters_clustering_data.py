#!/usr/bin/python

letters = "ABCFGHJKXYZ"

for i in range(1,len(letters)) :
    for j in range(0,i+1) :
	print i+1,j+1,abs( ord(letters[i]) - ord(letters[j]) )



