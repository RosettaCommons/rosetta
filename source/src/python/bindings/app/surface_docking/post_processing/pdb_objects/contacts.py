#!/usr/bin/python
import sys
import string
import calccontacts

Nargs = len(sys.argv)-1

if Nargs < 1:
    print 'Usage: contacts.py <pdb> parse=<parserule> mode=<mode> output=<output> <file>'
    print "Any order of arguments is acceptable and most arguments are optional\n"
    print 'pdb is mandatory\n'
    print "parserule has several options (optional; default = 'all')"
    print "'all':  whole protein; all chains and all heteroatoms"
    print "'A', 'AB', etc. - any chain(s) from the protein including heteroatoms of that chainID"
    print "'A-B', 'AB-C':  interface between chains separated by '-' "
    print "'A-ligand', 'all-ligand' - a subset or the whole protein against ligand (binding sites)\n"
    print "mode has several options (optional; default = 'fullatom')"
    print "fullatom - all atom-atom contacts"
    print "calpha - all CA-CA pairs within 7.0A"
    print "sccentroid - all sidechain centroid-centroid pairs within 6.0A\n"
    print "output can be either 'list' or 'profile' (optional; default = 'list')"
    print "'list' is a list of individual contacts"
    print "'profile' is table of total number of contacts for each residue"
    print "'both' will output both forms; you should use this only if you turn on 'file'\n"
    print "'file', if added to the arglist, will output a file for the type of output you specified"
    print "'pdb.contactmap' for 'list' or 'pdb.contactprofile' for 'profile'"
    sys.exit(0)

#default arguments
parserule = 'all'
output = 'list'
mode = 'fullatom'
dest = 'print'

arglist = sys.argv[1:Nargs+1]
for arg in arglist:
    if arg[-3:] == 'pdb':
        pdb = arg
        continue
    if arg[0:5] == 'parse':
        parserule = string.split(arg, '=')[1]
        continue
    if arg[0:4] == 'mode':
        mode = string.split(arg, '=')[1]
        continue
    if arg[0:6] == 'output':
        output = string.split(arg, '=')[1]
        continue
    if arg == 'file':
        dest = 'file'
        
contactObject = calccontacts.contactProtein(mode)
contactObject.initializePDB(pdb)
contactObject.find_contacts(parserule)

if output in ['list', 'both']:
    contactObject.outputMap(dest, pdb)
if output in ['profile', 'both']:
    contactObject.initContactProfile()
    contactObject.outputProfile(dest, pdb)
