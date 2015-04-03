#!/usr/bin/python
#David Masica


from pymol import cmd, preset
from os import popen, getcwd, chdir
from commands import getoutput
from string import strip
from time import sleep

me = getoutput('whoami')

"""
Master python script for utilities running rosetta SurfaceMode.
Each function is preceded by a header describing its usage

"""
print "Pymol add-ons for use with Rosetta SurfaceMode."
print "See instructions in ~/home/.pymolSurfaceMode.py for usage."

cmd.hide("everything", "all")

"""
Build Rosetta's periodic table
"""

cmd.alter("elem C",    "vdw=2.0000")
cmd.alter("elem N",    "vdw=1.7500")
cmd.alter("elem O",    "vdw=1.5500")
cmd.alter("elem S",    "vdw=1.9000")
cmd.alter("elem H",    "vdw=1.0000")
cmd.alter("resn LEU+VAL+GLY+ILE+PRO+MET+CYS+ALA+PHE and elem H and not name H", "vdw=1.2000")
cmd.alter("name Al",   "vdw=2.3972")
cmd.alter("name K",    "vdw=1.8712")
cmd.alter("name Omic", "vdw=1.7766")
cmd.alter("name Hydx", "vdw=1.7766")
cmd.alter("name Sica", "vdw=1.8532")

cmd.alter("name Pha",  "vdw=2.1500")
cmd.alter("name OPha", "vdw=1.7000")
cmd.alter("name OHha", "vdw=1.7700")
cmd.alter("name Ca2p", "vdw=1.3700")
cmd.alter("name Hha",  "vdw=0.9000")

cmd.alter("name CNH2",  "vdw=2.0000")
cmd.alter("name COO ",  "vdw=2.0000")
cmd.alter("name CH1 ",  "vdw=2.0000")
cmd.alter("name CH2 ",  "vdw=2.0000")
cmd.alter("name CH3 ",  "vdw=2.0000")
cmd.alter("name aroC",  "vdw=2.0000")
cmd.alter("name Ntrp",  "vdw=1.7500")
cmd.alter("name Nhis",  "vdw=1.7500")
cmd.alter("name NH20",  "vdw=1.7500")
cmd.alter("name Nlys",  "vdw=1.7500")
cmd.alter("name Narg",  "vdw=1.7500")
cmd.alter("name Npro",  "vdw=1.7500")
cmd.alter("name OH  ",  "vdw=1.5500")
cmd.alter("name ONH2",  "vdw=1.5500")
cmd.alter("name OOC ",  "vdw=1.5500")
cmd.alter("name S   ",  "vdw=1.9000")
cmd.alter("name Nbb ",  "vdw=1.7500")
cmd.alter("name CAbb",  "vdw=2.0000")
cmd.alter("name CObb",  "vdw=2.0000")
cmd.alter("name OCbb",  "vdw=1.5500")
cmd.alter("name Phos",  "vdw=2.1500")
cmd.alter("name Hpol",  "vdw=1.0000")
cmd.alter("name Hapo",  "vdw=1.2000")
cmd.alter("name Haro",  "vdw=1.2000")
cmd.alter("name HNbb",  "vdw=1.0000")

cmd.alter("name Hice",  "vdw=0.8000")
cmd.alter("name Oice",  "vdw=1.6000")
cmd.alter("name Ccod",  "vdw=1.7000")
cmd.alter("name Ocod",  "vdw=1.6000")


cmd.select("protein", "all  and not hetatm" )
cmd.select("surface", "hetatm and not resn doc")
cmd.disable("surface")
cmd.disable("protein")
cmd.show("cartoon")
cmd.set("cartoon_fancy_helices",1)
cmd.spectrum("count",selection="(protein)&e. c")
#preset.publication("protein")
cmd.show("spheres","surface")
cmd.clip("atoms",5,"all")
cmd.select("prot_int", "protein w. 5 of surface")
cmd.disable('prot_int')
cmd.select("surf_int", "surface w. 5 of protein")
cmd.disable('surf_int')
cmd.hide("everything","surface")
cmd.show("spheres"   ,"surf_int")
cmd.show("sticks","((byres (prot_int))&(!(n;c,o,h|(n. n&!r. pro))))")
cmd.space('cmyk')
cmd.bg_color('black')
cmd.select("BackBone", "name c+ca+o+n")
cmd.select("SC_Carbon", "elem c and not BackBone")
cmd.color("grey60","SC_Carbon")
cmd.disable('SC_Carbon')
############## OSTEOCALCIN STUFF ###############
cmd.color("red","name OE11")
cmd.color("red","name OE12")
cmd.color("red","name OE21")
cmd.color("red","name OE22")
cmd.remove("resn cys and name hg")
cmd.show("sticks","resn cys and not BackBone")
cmd.show("sticks","resn cys and name CA+HA")
###########################################
cmd.create("ProtObj","protein")
cmd.hide("everything","ProtObj")
cmd.show("surface", "ProtObj")
cmd.set('transparency',0.65)
#cmd.set('cartoon_transparency',0.2)
#cmd.color("grey50","protein")
cmd.color("phosphorus", "name P")
cmd.color("firebrick", "name o3p")
cmd.color("firebrick", "name o2p")
cmd.color("firebrick", "name o1p")
cmd.color("firebrick", "name og and resn sep")
cmd.color("phosphorus", "name Pha")
cmd.color("firebrick", "name OPha")
cmd.color("oxygen", "name OHha")
cmd.select("carbon", "elem C and surface")
cmd.color("gray", "carbon")
cmd.disable('carbon')
cmd.color("grey70","ProtObj")
cmd.set("spec_power",200)
cmd.set("spec_reflect",1.5)
cmd.set("depth_cue",0)

cmd.color("gray",      "name CNH2")
cmd.color("gray",      "name COO ")
cmd.color("gray",      "name CH1 ")
cmd.color("gray",      "name CH2 ")
cmd.color("gray",      "name CH3 ")
cmd.color("gray",      "name aroC")
cmd.color("blue",      "name Ntrp")
cmd.color("blue",      "name Nhis")
cmd.color("blue",      "name NH20")
cmd.color("blue",      "name Nlys")
cmd.color("blue",      "name Narg")
cmd.color("blue",      "name Npro")
cmd.color("red",        "name OH ")
cmd.color("red",       "name ONH2")
cmd.color("red",       "name OOC ")
cmd.color("orange",    "name S   ")
cmd.color("blue",      "name Nbb ")
cmd.color("gray",      "name CAbb")
cmd.color("gray",      "name CObb")
cmd.color("red",       "name OCbb")
cmd.color("red",       "name Oice")
cmd.color("phosphorus","name Phos")
cmd.color("white",     "name Hpol")
cmd.color("white",     "name Hapo")
cmd.color("white",     "name Haro")
cmd.color("white",     "name HNbb")
cmd.color("white",     "name Hice")
cmd.color("gray",      "name Ccod")
cmd.color("red",       "name Ocod")


util.performance(0)
def fcn(mk_file):
    file = open(mk_file, 'r').readlines()
    new_file = []
    for i in range(len(file)):
        new_file.append(strip(file[i]))

    return new_file

"""
Allows easy translation and rotation of protein with respect to surface.
This is useful when aligning protein and surface previous the start of a
Rosetta SurfaceMode run. When finished type save name.pdb and protein and
surface will be saved in single pdb. Note: It is easiest to use when protein
and surface are loaded as seperate molecules.
"""
# USAGE
# F1 = translate +1 angstrom
# F2 = translate -1 angstrom
# F3 = rotate +1 degree
# F4 = rotate -1 degree
# Note: Holding down any of the above keys results in a smooth continuous
#       motion
# Specify axis by simply typing x, y, or z
# Note: You can always reset the clipping plane and depth cue by typing z

cmd.alias("x","x()")
cmd.alias("y","y()")
cmd.alias("z","z()")

def x():
    print "translate or rotate about the x axis"
    print "F1 = translate positive"
    print "F2 = translate negative"
    print "F3 = rotate positive"
    print "F4 = rotate negative"
    cmd.set_key ('F1', cmd.translate,("[1,0,0]","protein"))
    cmd.set_key ('F2', cmd.translate,("[-1,0,0]","protein"))
    cmd.set_key ('F3', cmd.rotate,('x',1,"protein"))
    cmd.set_key ('F4', cmd.rotate,('x',-1,"protein"))

def y():
    print "translate or rotate about the y axis"
    print "F1 = translate positive"
    print "F2 = translate negative"
    print "F3 = rotate positive"
    print "F4 = rotate negative"
    cmd.set_key ('F1', cmd.translate,("[0,1,0]","protein"))
    cmd.set_key ('F2', cmd.translate,("[0,-1,0]","protein"))
    cmd.set_key ('F3', cmd.rotate,('y',1,"protein"))
    cmd.set_key ('F4', cmd.rotate,('y',-1,"protein"))


def z():
    print "translate or rotate about the z axis"
    print "F1 = translate positive"
    print "F2 = translate negative"
    print "F3 = rotate positive"
    print "F4 = rotate negative"
    cmd.set_key ('F1', cmd.translate,("[0,0,1]","protein"))
    cmd.set_key ('F2', cmd.translate,("[0,0,-1]","protein"))
    cmd.set_key ('F3', cmd.rotate,('z',1,"protein"))
    cmd.set_key ('F4', cmd.rotate,('z',-1,"protein"))
    cmd.clip("atoms",5,"all")
    cmd.set("depth_cue",0)


"""
Usefull once you have saved a pdb that contains both your protein and surface of interest
in the orientation you want. This function will then scp Rosetta's surface module to the
cluster with a copy of said pdb inside.
"""
# USAGE
# Simply type: remote_set_up directory, pdb
# directory = the name you would like to call your jazz:/home/directory
# Where pdb is the name of your pdb. Do not include the .pdb extension

def remote_set_up(directory, pdb):
    command1 = 'scp -r /home/cabalo1/SurfaceProject/Algorithms/scripts/surface_module jazz:' + directory
    command2 = 'scp ' + pdb + '.pdb jazz:' + directory
    popen(command1)
    popen(command2)

cmd.extend('remote_set_up', remote_set_up)

"""
Launchs Rosetta SurfaceMode jobs from pymol GUI
"""
# USAGE
# Simply type: launch directory, pdb
# directory = the name of the jazz:/home/directory/ that your pdb is in
# pdb = is the name of your pdb. Do not include the .pdb extension

def launch(directory, pdb, mode):
    command1 = 'xhost +'
    command2 = 'xterm -e ssh -X jazz \" export PATH=$PATH:/opt/condor-6.6.6/bin; export CONDOR_CONFIG=/opt/condor-6.6.6/etc/condor_config; cd ' + directory +  '; ./accessories/remote.launch ' + pdb + ' ' + mode + '\"'
    popen(command1)
    popen(command2)

cmd.extend('launch', launch)

"""
Gives user an xterm window inside the directory of interest. Usefull for checking progress of a run remotely
"""
# USAGE
# Simply type: status directory
# directory = the name of the jazz:/home/directory/ that your pdb is in

def status(directory):
    command1 = 'ssh jazz \"cd ' + directory + '; xterm\"'
    popen(command1)

cmd.extend('status', status)

"""
Goes to jazz, retrieves dock_surface structures, and loads into pymol. Designed for picking starting
structures from dock_surface mode
"""
# USAGE
# Simply type: get_start_structs source, destination
# source = a path to the directory containing dock_surfaceZZ directory you  are interested in
# desination = name of directory you would like to keep dock_surfaceZZ in on your home machine

def get_start_structs(source, destination):
    chdir( '/home/' + me )
    os.makedirs(destination,mode=0511)
    os.chmod( '/home/' + me + '/' + destination, 466 )
    command1 = 'ssh jazz \"scp -r ' + source + ' coltrane:' + destination + '/\"'
    popen(command1)
    #chdir( '/home/' + me )
    #chdir( destination + '/dock_surfaceZZ'  )

cmd.extend('get_start_structs', get_start_structs)

"""
Loads range of structures specified user, very useful for examining many structures quickly
"""

# USAGE
# Simply type: load begin, end
# begin = first structure in range of interest
# end = last structure in range of interest

def load_range(begin, end):
    begin = int(begin)
    begin = begin - 1
    end = int(end)
    cmd.reinitialize()
    cmd.bg_color('grey50')
    getoutput( 'ls | grep \".pdb\" > pdb_list' )
    pdb_list = open('pdb_list', 'r')
    pdb_list = fcn('pdb_list')
    getoutput( 'ls | grep \".pdb\" > pdb_list' )
    pdb_list = open('pdb_list', 'r')
    pdb_list = fcn('pdb_list')
    for i in range(begin, end):
        cmd.load(pdb_list[i])
        cmd.disable(pdb_list[i][:15])
    preset.publication('all and not hetatm')
    cmd.show('spheres','hetatm')
    cmd.clip("atoms",5,"all")
    cmd.set("depth_cue", 0)

cmd.extend('load_range', load_range)


