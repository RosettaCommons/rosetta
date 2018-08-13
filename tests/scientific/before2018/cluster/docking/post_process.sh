#!/bin/bash

echo filename > pdblist.txt
find . -name "*fasc" >>pdblist.txt
R --vanilla < PostProcess.R
