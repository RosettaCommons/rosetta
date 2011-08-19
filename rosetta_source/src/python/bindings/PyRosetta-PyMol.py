# Eric Kim
# script to communicate between PyRosetta and PyMol
# run rosetta libraries and this script inside PyMol command-line window
# PyRosetta and PyMol must be built with matching Python versions
# March 2009
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

with_pymol = True
try:
  import pymol
except:
  with_pymol = False

pymol_state = 1
def pymol_show(pose):
  global pymol_state

  if not with_pymol: return
  coord = {}
  for resi in range(1, pose.n_residue()+1):
    res = pose.residue(resi)
    resn = pose.pdb_info().number(resi)
    for atomi in range(1, res.natoms()+1):
      name = res.atom_name(atomi).strip()
      coord[(resn, name)] = res.atom(atomi).xyz()

  from pymol import cmd

  model = cmd.get_model("decoy")
  for a in model.atom:
    key = (int(a.resi), a.name.strip())
    if coord.has_key(key):
      a.coord=[coord[key].x,coord[key].y,coord[key].z]
    else:
      print key,"is missing"

  cmd.load_model(model,"decoy",pymol_state)
  cmd.forward()
  cmd.refresh()
#	cmd.orient("decoy",pymol_state)
  pymol_state += 1
  time.sleep(.1)


def pymol_load(f):
  if not with_pymol: return
  from pymol import cmd
  cmd.load(f, "decoy")
  cmd.show("cartoon", "decoy")


from rosetta import *
rosetta.init()

pose = Pose("hw4/ala.pdb")
pymol_load("hw4/ala.pdb")
pymol_show(pose)

scorefxn = core.scoring.ScoreFunction()
scorefxn.set_weight(core.scoring.fa_atr, 1.0)
scorefxn.set_weight(core.scoring.fa_rep, 1.0)
scorefxn.set_weight(core.scoring.hbond_sr_bb, 1.0)
scorefxn.set_weight(core.scoring.hbond_lr_bb, 1.0)
scorefxn.set_weight(core.scoring.hbond_bb_sc, 1.0)
scorefxn.set_weight(core.scoring.hbond_sc, 1.0)

#switch = SwitchResidueTypeSetMover("centroid")
#switch.apply(pose)
#scorefxn = create_score_function("score3")




import random, math, time

def perturb(new_pose):
  res = random.randrange(1,11)
  if random.randrange(0,2)==0:
    new_pose.set_phi(res, pose.phi(res)+random.gauss(0, 25))
  else:
    new_pose.set_psi(res, pose.psi(res)+random.gauss(0, 25))

from rosetta.core.fragment import *
fragset = core.fragment.ConstantLengthFragSet(3)
fragset.read_fragment_file("hw4/A10_3mers.frag")
movemap = core.kinematics.MoveMap()
movemap.set_bb(True)
mover_3mer = ClassicFragmentMover(fragset,movemap)
log = open("hw4/log", 'w')

def fold(pose):
  kt = 1

  low_pose = Pose()
  low_score = scorefxn(pose)
  new_pose = Pose()
  maxit = 1000
  start_time = time.time()
  stat = {"E":0,"A":0,"R":0}
  for it in range(0, maxit):
    if it==maxit/2: pose.assign(low_pose)

    score = scorefxn(pose)
    new_pose.assign(pose)

    if it<maxit/2: mover_3mer.apply(new_pose)
    else: perturb(new_pose)

    new_score = scorefxn(new_pose)
    delE =  new_score-score

    if delE<0:
      score = new_score
      pose.assign(new_pose)
      action = "E"
    else:
      if random.random()<math.exp(-delE/kt):
        score = new_score
        pose.assign(new_pose)
        action = "A"
      else: action = "R"

    if score<low_score:
      low_score = score
      low_pose.assign(pose)

    print it,"\t",score,"\t",low_score,"\t",action
    #log.write(str(score)+"\t"+str(low_score)+"\t"+str(action)+"\n")
    stat[action] += 1
    if action<>"R": pymol_show(pose)

  duration = time.time()-start_time
  print "took ",duration,"s\t",duration/float(maxit),"s/it"
  #print "E:",stat["E"]/float(maxit),"\t","A:",stat["A"]/float(maxit),"\t","R:",stat["R"]/float(maxit)
  print low_score
  pymol_show(low_pose)
  return low_pose

def angle(a):
  return math.fmod(a,180)

"""
for i in range(0, 1):
  pose2 = Pose()
  pose2.assign(pose)
  low_pose = fold(pose2)
  print i,
  for i in range(1,11):
    log.write(str(angle(low_pose.phi(i)))+"\t"+str(angle(low_pose.psi(i)))+"\n")
"""

import threading

def run():
  thread = threading.Thread(target=fold, args=(pose,))
  thread.setDaemon(1)
  thread.start()
