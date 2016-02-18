
# coding: utf-8

# In[1]:

import os
import sys
from collections import defaultdict


# In[12]:

codes = defaultdict()
PDB_NAMES = open("../../input_output/3-letter_codes/pdb_sugar.codes", 'r')
for line in PDB_NAMES:
    if line[0] == '#':continue
    line = line.strip()
    if not line: continue
    lineSP = line.split()
    if len(lineSP) < 3: continue
    
    code = lineSP[0]
    short_iupac = lineSP[2]
    codes[code] = short_iupac
PDB_NAMES.close()
codes


# In[14]:

ORIGINAL = open("original_glycan_relax_conformers_unformatted.txt", 'r')
NEW = open("original_glycan_relax_conformers_formatted.txt", 'w')

for line in ORIGINAL:
    if line[0] == '#': continue
    line = line.strip()
    if not line: continue
        
    print line
    lineSP = line.split()
    
    
    res2_name = lineSP[1]
    res2_link = lineSP[2]
    
    res1_name = lineSP[3]
    res1_link = lineSP[4]
    
    if not codes.has_key(res1_name):
        codes[res1_name] = "->?)-Unknown-"+res1_name
    if not codes.has_key(res2_name):
        codes[res2_name] = "->?)-Unknown-"+res2_name
        
    columns = [codes[res1_name].replace('?', res1_link), codes[res2_name].replace('?', res2_link)]
    columns.extend(lineSP[5:])
    
    new_line = "\t".join(columns)
    print new_line
    NEW.write(new_line+"\n")
ORIGINAL.close()
NEW.close()


# In[3]:

x = "A?B"


# In[4]:

x.replace('?', "Happy")


# In[ ]:



