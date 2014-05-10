#!/bin/bash
# Edited by Emily Koo to create difference plot if native PDB is included in dir

echo "Generating contact map"
ls *.pdb > PDBsInDir
DecoyCnt=$( wc -l PDBsInDir | cut -f1 -d" " )

Current_structure=1
while [ $Current_structure -le $DecoyCnt ]
  do
  Current_structure_name=$( head -$Current_structure PDBsInDir | tail -1 )
  grep ATOM $Current_structure_name > protein.pdb
  contacts.py protein.pdb > $Current_structure_name.con
  rm protein.pdb
  Current_structure=$[ $Current_structure + 1 ] 
done

cat AdsState*.con > Ads.con

# Should work if there's either one or two types of files
ContactMap.py Ads.con native.pdb.con 

rm *con PDBsInDir
    
echo "Contact map DONE"
