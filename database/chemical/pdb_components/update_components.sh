#!/bin/bash

curl -o ./components.cif.gz https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz 
gunzip components.cif.gz
rm components.*.cif
python split_components.py
rm components.cif
