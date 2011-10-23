#! /usr/bin/env python
# Adapted from fitting.py (Copyright (c) 2005 Robert L. Campbell)

from pymol import cmd

def align_features(obj1,select1,obj2,select2, debug=False):
  """
DESCRIPTION

  "align_features" aligns obj1 and obj2 by aligning select1 and select1

  Be careful when creating your selection strings. Within the
  selections, do not include the object name because the chain, residue
  name, residue number etc. of selection1 of object1 are converted to
  match those in selection2.  If the object names are included in the
  selections, no atoms will be selected since an atom cannot exist in
  both object1 and object2 at the same time.

USAGE 

  align_features object1, selection1, object2, selection2

  DO NOT include object names in selections!

EXAMPLES
  ####### test_align_features.pml ##########
  # copy into test_align_features.pml in a directory with 1r6jH.pdb and align_features.py
  #
  # run by: 'pymol align_features.py'
  #

  run align_features.py

  load 1r6jH.pdb, 1r6jH.pdb_332
  cmd.select("inst_332", "/1r6jH.pdb_332//A/222/H or /1r6jH.pdb_332//A/222/N or /1r6jH.pdb_332//A/222/CA")
  orient inst_332

  load 1r6jH.pdb, 1r6jH.pdb_445
  cmd.select("inst_445", "/1r6jH.pdb_445//A/220/H or /1r6jH.pdb_445//A/220/N or /1r6jH.pdb_445//A/220/CA")
  align_features 1r6jH.pdb_332, inst_332, 1r6jH.pdb_445, inst_445

  """

  m=cmd.get_model("%s & %s" % (obj1,select1))
  n=cmd.get_model("%s & %s" % (obj2,select2))

  if len(m.atom) != len(n.atom):
    print "Unable to preform fit because the selections have unequal number of atoms:"
    print "  select1=%s & %s, has %s atoms." % (obj1, select1, len(m.atom))
    print "  select2=%s & %s, has %s atoms." % (obj2, select2, len(n.atom))
    return
  else:
    total = len(m.atom)

  # for the atoms to be used in fit:
  list_m, list_n = [], []
  for at in m.atom:
    list_m.append((at.id,at.chain,at.resn,at.resi,at.name,at.segi, at.alt))
  for at in n.atom:
    list_n.append((at.id,at.chain,at.resn,at.resi,at.name,at.segi, at.alt))


  # the fit command is busted: it requires the things being aligned to
  # have the same identification, ie same chain, resn, resi, segi and
  # alt.

  # set a new segi for the atoms to be used in fit command and to
  # allow resetting later
  seg_fit="1fit"

  # change the chain,resn,resi,segi and alt of select1 to match
  # select2
  for i in range(total):
    cmd.do("alter %s & id %s, chain='%s'" % (obj1,list_m[i][0],list_n[i][1]))
    cmd.do("alter %s & id %s, resn='%s'" % (obj1,list_m[i][0],list_n[i][2]))
    cmd.do("alter %s & id %s, resi=%s" % (obj1,list_m[i][0],list_n[i][3]))
    cmd.do("alter %s & id %s, segi='%s'" % (obj1,list_m[i][0],seg_fit))
    cmd.do("alter %s & id %s, alt='%s'" % (obj1,list_m[i][0],list_n[i][6]))
    # change the segid for obj2 and select2
    cmd.do("alter %s & id %s, segi='%s'" % (obj2,list_n[i][0],seg_fit))

  if debug:
    print "Fitting %s and %s\n     to %s and %s" % (obj1,select1,obj2,select2)
  
    print "Altered to:"
    print "%s & %s & segi %s\n" % (obj1,select1,seg_fit),
    print "%s & %s & segi %s\n" % (obj2,select2,seg_fit),
    print "--------------------------------------------\n"

  rms = cmd.align(
    "%s & %s & segi %s" % (obj1,select1,seg_fit),
    "%s & %s & segi %s" % (obj2,select2,seg_fit),
    quiet=0)

  if debug:
    cmd.delete("%s_fitting" % obj1)
    cmd.delete("%s_fitting" % obj2)
    # create new objects to show the fit atoms
    cmd.create("%s_fitting" % obj1, "%s & %s & segi %s" % (obj1,select1,seg_fit))
    cmd.create("%s_fitting" % obj2, "%s & %s & segi %s" % (obj2,select2,seg_fit))

  # reset chain,resn,resi,segi & alt of obj1 & select1 from stored list
  for atoms_m in list_m:
    cmd.do("alter %s & id %s, chain='%s'" % (obj1,atoms_m[0],atoms_m[1]))
    cmd.do("alter %s & id %s, resn='%s'" % (obj1,atoms_m[0],atoms_m[2]))
    cmd.do("alter %s & id %s, resi=%s" % (obj1,atoms_m[0],atoms_m[3]))
    cmd.do("alter %s & id %s, segi='%s'" % (obj1,atoms_m[0],atoms_m[5]))
    cmd.do("alter %s & id %s, alt='%s'" % (obj1,atoms_m[0],atoms_m[6]))
  # reset segi of obj2 & select2 from stored list
  for atoms_n in list_n:
    cmd.do("alter %s & id %s, segi='%s'" % (obj2,atoms_n[0],atoms_n[5]))

  if debug:
    print "RMSD for fitting selection %s of %s onto \n                 selection %s of %s = %6.3f" % (select1, obj1, select2, obj2, rms)

cmd.extend("align_features", align_features)
