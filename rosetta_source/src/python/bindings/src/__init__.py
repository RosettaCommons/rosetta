# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.  # (c) This file is part of the Rosetta software suite and is made available under license.  # (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.  # (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os, sys, platform, os.path

import warnings
warnings.filterwarnings("ignore", "to-Python converter for .+ already registered; second conversion method ignored.", RuntimeWarning, "^rosetta\\.")

# Create global 'Platform' that will hold info of current system
if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == 'win32'  : Platform = 'windows'
else: Platform = "_unknown_"
PlatformBits = platform.architecture()[0][:2]

if Platform == "cygwin": print 'Importing Rosetta libs, this could take more then a few minutes on old hardware...'

# double-checked right order...
import utility
import utility.excn

import numeric

import basic
import basic.datacache

import core
import core.graph
import core.chemical.orbitals
import core.scoring
import core.scoring.methods
import core.scoring.constraints
import core.scoring.etable
import core.kinematics

import core.io.silent
import core.pose


import protocols
import protocols.moves
import protocols.canonical_sampling
import protocols.simple_moves
import protocols.jumping
import protocols.jd2.archive
import protocols.abinitio

import protocols.filters
import protocols.docking
import protocols.init

import rosetta.utility.file
import rosetta.core.chemical
import rosetta.core.conformation
import rosetta.core.id
import rosetta.core.io
import rosetta.core.io.pdb
import rosetta.core.fragment
import rosetta.core.kinematics
import rosetta.core.pack
import rosetta.core.pack.task
import rosetta.core.pose
import rosetta.core.scoring.hbonds
import rosetta.protocols.loops
import rosetta.protocols.wum
import rosetta.protocols.relax
import rosetta.core.pose.signals

import rosetta.protocols.simple_moves

from rosetta.core.id import *
from rosetta.core.conformation import *
from rosetta.core.chemical import *
from rosetta.core.pose import Pose

from rosetta.core.import_pose import pose_from_pdb
from rosetta.core.io.pdb import dump_pdb
from rosetta.core.pose import make_pose_from_sequence

from rosetta.core.scoring import *
from rosetta.core.kinematics import *
from rosetta.core.fragment import *
from rosetta.core.pack.task import *
from rosetta.core.pack.task.operation import *

from rosetta.protocols.moves import *
from rosetta.protocols.simple_moves import *
from rosetta.protocols.abinitio import *
from rosetta.protocols.docking import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import *

from rosetta.protocols.rigid import *
from rosetta.protocols.simple_moves import *


class CD:
    ''' Class to represent named tuples '''
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = '|'
        for i in dir(self):
            if not i.startswith('__'): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


# add iter property to Pose
def _Pose_residue_iterator(obj):
    def __pose_iter():
        for i in range(obj.total_residue()): yield obj.residue(i+1)
    return __pose_iter()

Pose.__iter__ = _Pose_residue_iterator

#    def __iter__(self):
#        def __pose_iter():
#            for i in range(self.total_residue()): yield self.residue(i+1)
#        return __pose_iter()



def add_extend(vectype):
  def extendfunc(vec,othervec):
    for i in othervec: vec.append(i)

  vectype.extend = extendfunc


for k,v in utility.__dict__.items():
  if k.startswith("vector1_"):
    add_extend(v)


def new_vector1_init(self,arg1=None,arg2=False):
  self.__old_init()
  if hasattr(arg1,"__iter__"): self.extend(arg1)
  elif type(arg1) is type(1):
    for i in xrange(arg1):
      self.append(arg2)


def replace_init(cls,init):
  cls.__old_init = cls.__init__
  cls.__init__ = init


def Vector1(l):
    ''' Create Vector1 object deducing type from the given list
    '''
    if   all( map(lambda x: type(x) == int, l) ): t = utility.vector1_int
    elif all( map(lambda x: type(x) == float or type(x) == int, l) ): t = utility.vector1_double
    elif all( map(lambda x: type(x) == str, l) ): t = utility.vector1_string
    elif all( map(lambda x: type(x) == bool, l) ): t = utility.vector1_bool
    elif all( map(lambda x: type(x) == core.id.AtomID, l) ): t = utility.vector1_AtomID
    else: raise Exception('Vector1: attemting to create vector of unknow type or mixed type vector init_list=' + str(l) )

    v = t()
    for i in l: v.append(i)
    return v


class PyJobDistributor:
	def __init__( self, pdb_name, nstruct, scorefxn ):
		self.pdb_name = pdb_name
		self.nstruct = nstruct
		self.current_num = 0		#current decoy number
		self.current_name = " "		#current decoy name
		self.job_complete = False	#job status
		self.scorefxn = scorefxn	#used for final score calculation
		self.native_pose = 0		#used for rmsd calculation
		self.additional_decoy_info = ' '  #used for any additional decoy information you want stored
		self.start_decoy()		#initializes the job distributor

	def start_decoy( self ):
		if self.job_complete == True:
			return
		i = 1
		file_exists = True
		while (file_exists == True and i <= self.nstruct):
			current_name = self.pdb_name + "_" + str(i) + ".pdb"
			if os.path.exists(current_name) == False:
				current_name_temp = current_name + ".in_progress"
				if os.path.exists(current_name_temp) == False:
					file_exists = False	#if such a file is not found, i is the current decoy #
					f = open(current_name_temp, 'w')
					f.write("this decoy is in progress")
					f.close()
					self.current_name = current_name
			i = i + 1
		self.current_num = i - 1
		if (file_exists == True):
			self.job_complete = True

	def output_decoy( self, pose):
		current_name_temp = self.current_name + ".in_progress"
		if os.path.exists(current_name_temp) == False:
			return

		dump_pdb(pose, self.current_name) #outputs pdb file
		os.remove(current_name_temp)

		score_tag = ".fasc"
		if (pose.is_fullatom() == False):
			score_tag = ".sc"

		scorefile = self.pdb_name + score_tag
		if os.path.exists(scorefile) == False:
			f = open(scorefile, 'w')
			f.write("pdb name: " + self.pdb_name + "     nstruct: " + str(self.nstruct) + '\n')
			f.close

		score = self.scorefxn(pose)	#calculates total score
		score_line = pose.energies().total_energies().weighted_string_of( self.scorefxn.weights())
		output_line = "filename: " + self.current_name + " total_score: " + str(round(score,2))
		if (self.native_pose != 0 ):	#calculates an rmsd if a native pose is defined
			rmsd = CA_rmsd(self.native_pose, pose )
			output_line = output_line + " rmsd: " + str(round(rmsd,2))
		f = open(scorefile, 'a')
		f.write(output_line + ' ' + score_line + self.additional_decoy_info + '\n') #outputs scorefile
		f.close

		self.start_decoy()

def generate_nonstandard_residue_set( params_list ):
	res_set = ChemicalManager.get_instance().nonconst_residue_type_set("fa_standard")
	res_set.read_files(params_list, ChemicalManager.get_instance().atom_type_set("fa_standard"),
                                    ChemicalManager.get_instance().element_set('fa_standard'),
	                                ChemicalManager.get_instance().mm_atom_type_set("fa_standard"),
                                    ChemicalManager.get_instance().orbital_type_set("fa_standard"),
	                   )
    	return res_set

def standard_task_factory():
	tf = TaskFactory()
	tf.push_back(InitializeFromCommandline())
	#tf.push_back(IncludeCurrent())
	tf.push_back(NoRepackDisulfides())
	return tf

def standard_packer_task(pose):
	tf = standard_task_factory()
	task = tf.create_task_and_apply_taskoperations(pose)
	return task

def add_extra_options():
    rosetta.protocols.abinitio.AbrelaxApplication.register_options()
    rosetta.protocols.abinitio.IterativeAbrelax.register_options()
    rosetta.protocols.abinitio.register_options_broker()

################################################################################
# TOOLBOX
# these methods are useful in PyRosetta and intended to demonstrate
# proper syntax for various common activities
# for those interested, Rosetta hag a Surface Area calculator and a Radius of
# Gyration calculator, they are EnergyMethods sa and rg respectively,
# create an empty ScoreFunction and use ScoreFunction.set_weight to use
# these calculations
# other common geometry calculators are CA_rmsd and all_atom_rmsd

# EVAN CHECK THIS
def generate_resfile_from_pose( pose , resfilename , input_sc = True ):
    """
    Writes a resfile for  <pose>  named  <resfilename> , optionally allowing
    input side chains to be used in packing

    example:
        generate_resfile_from_pose(pose,'1YY8.resfile')
    See also:
        Pose
        PackRotamersMover
        TaskFactory
    """
    f = open(resfilename, 'w')
    id = "NATRO"
    start = ''
    if input_sc:
        start = 'USE_INPUT_SC\n'
    f.write(start+'start\n')
    info = pose.pdb_info()
    # pose_from_sequence returns empty PDBInfo, Pose() makes NULL
    if info and info.nres():
        for i in range (1,pose.total_residue()+1):
            num = pose.pdb_info().number(i)
            chain = pose.pdb_info().chain(i)
            f.write(str(num).rjust(4) + str(chain).rjust(3) + str(id).rjust(7) + '  \n')
    else:
        for i in range (1,pose.total_residue()+1):
            num = i
            chain = ' '
            f.write(str(num).rjust(4) + str(chain).rjust(3) + str(id).rjust(7) + '  \n')
    f.close()

def generate_resfile_from_pdb( pdbfilename , resfilename , input_sc = True ):
	"""
	Writes a resfile for PDB file  <pdbfilename>  named  <resfilename> ,
	optionally allowing input side chains to be used in packing

	example:
	    generate_resfile_from_pdb('1YY8.pdb','1YY8.resfile')
	See also:
	    Pose
	    PackRotamersMover
	    TaskFactory
	"""
	p = pose_from_pdb(pdbfilename)
	generate_resfile_from_pose(p,resfilename,input_sc)

# replaces the residue at  <resid>  in  <pose>  with  <new_res>
def mutate_residue(pose, resid, new_res):
	"""
	Replaces the residue at  <resid>  in  <pose>  with  <new_res>
	note: <new_res>  is the single letter name for the desired ResidueType

	example:
	    mutate_residue(pose,30,A)
	See also:
	    Pose
	    PackRotamersMover
	"""
	if (pose.is_fullatom() == False):
		IOError("mutate_residue only works with fullatom poses")

	scorefxn = create_score_function('standard')
	pack_task = TaskFactory.create_packer_task(pose)
	pack_task.initialize_from_command_line()

	v1 = rosetta.utility.vector1_bool()
	mut_res = aa_from_oneletter_code(new_res)

	for i in range(1,21):
		if (i == mut_res):
			v1.append(True)
		else:
			v1.append(False)

	for i in range(1,pose.total_residue()+1):
		if (i != resid):
			pack_task.nonconst_residue_task(i).prevent_repacking()

	pack_task.nonconst_residue_task(resid).restrict_absent_canonical_aas( v1 )

	packer = protocols.simple_moves.PackRotamersMover(scorefxn, pack_task)
	packer.apply(pose)
	return pose

import urllib

# removes non ATOM lines from  <pdb_file>  and writes to  <pdb_file>.clean.pdb
def cleanATOM( pdb_file , edit = -4 ):
    """
    Writes a PDB file from  <pdb_file> with all lines not beginning with
    ATOM removed tp  <pdb_file>.clean.pdb
    note: the second argument, <edit>, if for PDB files not ending in .pdb

    example:
        cleanATOM('1YY9.pdb')
    See also:
        Pose
        Pose.dump_pdb
        pose_from_pdb
        pose_from_rcsb
        cleanCRYS
    """
    if not edit:
        edit = 255
    # why is it pdb_file[:-4] the whole way?
    if os.path.exists( os.getcwd() + '/' + pdb_file ):
        print 'if the file',pdb_file[:edit]+'.clean.pdb already exists, it will be overwritten'
#        os.system("grep \"ATOM\" %s.pdb > %s.clean.pdb"%(pdb_file[:edit],pdb_file[:edit]))
        fid = open(pdb_file,'r')
        data = fid.readlines()
        fid.close()
        good = []
        for i in data:
            if i[:5] == 'ATOM ':
                # maybe rules for ligands and DNA and/or water...?
                good.append(i)
        fid = open(pdb_file[:edit]+'.clean.pdb','w')
        fid.writelines(good)
        fid.close()
        print 'PDB',pdb_file,'successfully cleaned, non-ATOM lines removed\nclean data written to',pdb_file[:edit]+'.clean.pdb'
    else:
        raise IOError('No such file or directory named '+pdb_file)

# removes redundant crystal contacts, isolate monomer
def cleanCRYS( pdb_file , olig = 2 ):
    """
    Writes a PDB file for a monomer of  <pdb_file>  if it is a  <olig>-mer
    to  <pdb_file>.mono
    note: this is a simple sequence comparison

    example:
        cleanCRYS('1YY8.pdb',2)
    See also:
        Pose
        Pose.dump_pdb
        pose_from_pdb
        pose_from_rcsb
        cleanATOM
    """
    if os.path.exists( os.getcwd() + '/' + pdb_file ):
        print 'if the file',pdb_file[:-4]+'.mono.pdb already exists, it will be overwritten'
        pose = pose_from_pdb(pdb_file)
        tot = pose.total_residue()
        seq = pose.sequence()
        frags = ['']*olig
        match = [False]*(olig-1)
        olig = float(olig)
        frac = int(round(tot/olig))
        for f in range(int(olig)):
            frags[f] = seq[:frac]
            seq = seq[frac:]
        for f in range(int(olig-1)):
            match[f] = (frags[0]==frags[f+1])
        if sum(match)==(olig-1):
           for i in range(frac*int(olig-1)):
               pose.delete_polymer_residue(frac+1)
           pose.dump_pdb(pdb_file[:-4]+'.mono.pdb')
           print 'PDB',pdb_file,'successfully cleaned, redundant monomers removed\nmonomer data written to',pdb_file[:-4]+'.mono.pdb'
        else:
            print pdb_file,'is not a '+str(olig)+'-mer'
    else:
        raise IOError('No such file or directory named '+pdb_file)

# returns a pose made from seq, all phi/psi/omega set to 180
def pose_from_sequence( seq , res_type = 'fa_standard' ):
    """
    Returns a pose generated from amino acid single letters in  <seq>  using
    the  <res_type>  ResidueType

    example:
        pose=pose_from_sequence('THANKSEVAN')
    See also:
        Pose
        make_pose_from_sequence
        pose_from_pdb
        pose_from_rcsb
    """
    pose=Pose()
    make_pose_from_sequence(pose,seq,res_type)
    for i in range(0,pose.total_residue()):
        pose.set_phi(i+1,180)
        pose.set_psi(i+1,180)
        pose.set_omega(i+1,180)
    pose.pdb_info(rosetta.core.pose.PDBInfo())
    pose.pdb_info().name(seq[:4])
#	print pose
    return pose

# retreives pdbData from rcsb for  <pdb_code>
# ADD NAMING OPTION
def load_from_rcsb( pdb_code ):
    """
    Writes PDB data for RCSB data for  <pdb_code>  into the file  <pdb_code>.pdb

    example:
        load_from_rcsb('1YY8')
    See also:
        Pose
        pose_from_pdb
        pose_from_rcsb
        pose_from_sequence
        cleanATOM
        cleanCRYS
    """
    if pdb_code:    # if something input...
        pdb_code = pdb_code.upper()
        try:
            filename = urllib.urlretrieve('http://www.rcsb.org/pdb/files/' + pdb_code + '.pdb')[0]
        except:
            raise IOError('Cannot access the PDB database, please check your Internet access')
        else:
            if (os.path.getsize(filename) > 1500):    # arbitrary 1500...then pdb_code was invalid
                pdb_file = open(filename)
                pdb_data = pdb_file.readlines()
                pdb_file.close()

                pdb_code = pdb_code + '.pdb'
                if os.path.exists( os.getcwd() + '/' + pdb_code ):
                    print 'the file',pdb_code,'already exists, this file will be overwritten'
                #if input('Do you want to overwrite ' + pdbCode + '.pdb')
                pdb_file = open(pdb_code,'w')
                pdb_file.writelines(pdb_data)
                pdb_file.close()

                print 'PDB',pdb_code[:-4],'successfully loaded from rcsb into',pdb_code
#                if auto_clean:
#                    cleanATOM(pdb_code)
            else:
                raise IOError('Invalid PDB code')
        os.remove(filename)    # remove temp file

# packaged my method to fit naming
def pose_from_rcsb( pdb_code , ATOM = True , CRYS = True ):
    """
    Returns a pose for RCSB PDB  <pdb_code> , also writes this data to
    <pdb_code>.pdb, optionally calls cleanATOM and cleanCYRS

    example:
        pose=pose_from_rcsb('1YY8')
    See also:
        Pose
        pose_from_pdb
        pose_from_sequence
        load_from_rcsb
        cleanATOM
        cleanCRYS
    """
    load_from_rcsb(pdb_code)
    if ATOM:
        cleanATOM(pdb_code+'.pdb')
        pdb_code = pdb_code+'.clean'
    if CRYS:
        cleanCRYS(pdb_code+'.pdb')
        pdb_code = pdb_code+'.mono'
    pose = pose_from_pdb(pdb_code+'.pdb')
    return pose

# fills secstruct into pose...some wierd stuff, DSSP
def get_secstruct( pose , space = 8 , page = 80 ):
    """
    Predicts the secondary structure of  <pose> , loading this data into
    the pose's secstruct information and printing the prediction to screen

    example:
        get_secstruct(pose)
    See also:
        Pose
        Pose.secstruct
        Pose.sequence
        pose_from_pdb
        pose_from_sequence
        pose_from_rcsb
    """
    dssp = rosetta.protocols.jumping.Dssp(pose)
    dssp.insert_ss_into_pose(pose)
    seq = pose.sequence()
    sec = pose.secstruct()
    count = 0
    while len(seq):
        count += 1
        num = ''
        for i in range(page*(count-1)+1,page*count,space):
            num += str(i).ljust(space)
        print num+'\n'+seq[:page]+'\n'+sec[:page]+'\n'
        seq = seq[page:]
        sec = sec[page:]

# returns an HBondSet filled with the hbonding info
def get_hbonds( pose ):
    """
    Returns a HBondSet of the hydrogen bonding in  <pose>
    note: more info in rosetta.core.scoring.hbonds

    example:
        hbset=get_hbonds(pose)
        hbset.hbond(1).acc_res()
    See also:
        Pose
        Energies
        HBondSet
        ScoreFunction
        create_score_function
    """
    if not pose.energies().energies_updated():
        raise IOError('Energy is not updated, please score the pose first!')
    hbset = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set(pose,False,hbset)
    return hbset


# end of TOOLBOX
################################################################################


# PyMOL link code ----------------------------------------------------------------------------------
import uuid, socket, bz2
from array import array

#if Platform == "cygwin":
#        import gzip
#else:
#        import bz2

class PySocketClient:
    def __init__(self, udp_port=65000, udp_ip = '127.0.0.1'):
        self.udp_ip, self.udp_port = udp_ip, udp_port

        # of course this will not works... but it will be readjusted automatically
        self.last_accepted_packet_size = 1024*8  # ... well maybe next time...
        self.uuid = uuid.uuid4()
        self.sentCount = array('H', [0, 0, 0])  # packet_id, N, count
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    def sendMessage(self, msg):
        count = 1
        if len(msg) > self.last_accepted_packet_size:
            count = len(msg) /  self.last_accepted_packet_size + 1

        #print 'Sending messgae in %s packets...' % count

        for i in range(count):
            self.sentCount[1] = i
            self.sentCount[2] = count
            m = msg[i*self.last_accepted_packet_size:(i+1)*self.last_accepted_packet_size]
            self.sendRAWMessgae(m)

        self.sentCount[0] += 1


    def sendRAWMessgae(self, msg):
        buf = array('c', self.uuid.bytes)
        buf.extend( self.sentCount.tostring() )
        #buf.extend( array('H', [1,2]).tostring() )
        buf.extend( msg )

        self.socket.sendto(buf, (self.udp_ip, self.udp_port) )


#PyMOL_Mover = rosetta.protocols.moves.PyMolMover  # for now we use one implementation for both Rosetta C++ and PyRosetta
#__Deprecated_do_not_use_
class PyMOL_Mover(rosetta.protocols.moves.PyMolMover):
    def __init__(self, *args): rosetta.protocols.moves.PyMolMover.__init__(self, *args)

    def send_colors(self, pose, colors, default_color=protocols.moves.XC_blue):
        cm = rosetta.utility.map_int_int()
        for c in colors: cm[c] = colors[c]

        rosetta.protocols.moves.PyMolMover.send_colors(self, pose, cm, default_color)



class WinPyMOL_Mover(rosetta.protocols.moves.PyMolMover):
    def __init__(self, keep_history=False, update_energy=False, energy_type='total_score'):
        rosetta.protocols.moves.PyMolMover.__init__(self)
        self.keep_history = keep_history
        self.update_energy = update_energy
        self.energy_type = energy_type
        self.link = PySocketClient()

    def getPoseName(self, pose):
        #if 'name' in pose.__dict__: return pose.name[:255]
        if not pose.pdb_info(): return 'pose'
        else:
            p1,p2,p3 = pose.pdb_info().name().rpartition('.pdb')
            name = p1 or p3
            return name[:255]


    def apply(self, pose):
        #print 'This PyMOL_Mover apply...', pose

        name = self.getPoseName(pose)
        #print name

        # Creating message...
        os = rosetta.utility.OStringStream();  pose.dump_pdb(os)
        #if Platform == "cygwin":
        #        message = 'PDB.gzip' + chr(len(name)) + name + gzip.compress(os.str())
        #else:
        message = 'PDB.bz2 ' + chr(self.keep_history) + chr(len(name)) + name + bz2.compress(os.str())
        self.link.sendMessage(message)

        if self.update_energy: self.send_energy(pose , self.energy_type)



    def _send_RAW_Energies_old(self, pose, energyType, enegies):
        energyType = energyType[:255]
        name = self.getPoseName(pose)
        message = 'Ener.bz2' + chr(self.keep_history) + chr(len(name)) + name \
                             + chr(len(energyType)) + energyType \
                             + bz2.compress( array('f', enegies).tostring())
        self.link.sendMessage(message)


    def _send_RAW_Energies(self, pose, energyType, energies, autoscale=True):
        energyType = energyType[:255]
        name = self.getPoseName(pose)

        if autoscale: energies = self._scale_energy(energies)

        e = ''
        info = pose.pdb_info()
        # pose_from_sequence fills PDBInfo object but nres=0, Pose makes NULL
        if info and info.nres():
            for i in range(len(energies)):
                chain = info.chain(i+1)[0]
                res = info.number(i+1)
                icode = info.icode(i+1)
                e += '%s%4d%c%02x' % (chain, res, icode, energies[i])
        else:
            for i in range(len(energies)):
                chain = ' '
                res = i+1
                e += '%s%4d %02x' % (chain, res, energies[i])

        message = 'Ener.bz2' + chr(self.keep_history) + chr(len(name)) + name \
                             + chr(len(energyType)) + energyType \
                             + bz2.compress(e)
        self.link.sendMessage(message)


    def _scale_energy(self, energies):
        ''' scale give array from 0 to 255
        '''
        #for i in range(len(energies)):
        #    energies[i] = i*1./len(energies)

        r = [0]*len(energies)
        mi = min(energies)
        ma = max(energies)
        if ma - mi < 1e-100: ma+= 1e-100
        for i in range(len(energies)):
            r[i] = int( (energies[i] - mi)*255. / (ma-mi) )
        return r


    # energy output functions
    ################################################################################
    def send_energy(self, input_pose, energy_type='total_score'):
        ''' Send cummulative energy to PyMOL
        '''
        if not input_pose.energies().energies_updated():
            raise IOError('PyMOL_Mover::send_specific_energy: Energy is not updated, please score the pose first!')

        #self.apply(input_pose)
        output = [0.]*input_pose.total_residue()
        score_type = score_type_from_name(energy_type)
        weight = 1.
        if not energy_type == 'total_score':
            weight = input_pose.energies().weights()[score_type]
        for i in range( 0 , len(output) ):
            output[i] = input_pose.energies().residue_total_energies(i+1)[score_type]*weight
        self._send_RAW_Energies(input_pose,energy_type,output)


    def send_colors(self, pose, colors, default_color='blue'):
        ''' Color protein by using color dict as map resiue_num --> color_name
        '''
        energies = [ X11Colors[default_color][0] ] * pose.total_residue()
        for r in colors: energies[r-1] = X11Colors[ colors[r] ][0]
        self._send_RAW_Energies(pose, 'X11Colors', energies, autoscale=False)


    # returns all variables in workspace matching target
    #    main_vars must be globals() or vars() when called
    #    i.e. pose=Pose(), target=pose will return matches=['pose']
    def _find_variable_name( target , main_vars ):
        # searches 'main_vars' for the object matching 'target', returns all matches
        matches = [0]
        matches = [index for index, value in main_vars.items() if value == target]
        if len(matches) > 1:
            print 'Consider reassignment, multiple variables have this object data'
        return matches

################################################################################
###############################

# Workable but unfinished luxury PyMOL_Mover features

    # a generic message sender
    #    type_name tags the data, how it is to be interpreted
    #    size sets the length of data units, must match and be single digit
    def send_any( self , ptype , pose , data , size = 6 ):
        # type_name MUST be the same on both sides, same on pymol side
        to_send = str(size)
        # size MUST be a single digit!
        # all data[i] MUST be at least size long!
        info = pose.pdb_info()
        for i in range( 1 , pose.total_residue() + 1 ):
            if info and info.nres():
                pdb = (str(info.number(i))+info.icode(i)+info.chain(i)).rjust(6)
            else:
                pdb = (str(i)+' A').rjust(6)
            dat = str(data[i-1])[0:size]
            if not len(dat)==size:
                raise IOError('Error: not all data is the same size!')

            to_send += '%s%s' % (pdb, dat)
        name = self.getPoseName(pose)
        message = ptype + chr(self.keep_history) + chr(len(name)) + name + bz2.compress(to_send)
        self.link.sendMessage(message)


    # polar identity per residue
    def send_polars( self , pose ):
        data = [0]*pose.total_residue()
        # send 0 or 1, if polar or not
        for i in range( 1 , pose.total_residue() + 1 ):
            data[i-1] = int(pose.residue(i).is_polar())
        self.send_any( 'pol.bz2 ' , pose , data , 1 )

    # movemap dof info per residue, pose limits the length viewed
    def send_movemap( self , pose , movemap ):
        data = [0]*pose.total_residue()
        # sent oddly...11 first digit for bb, second for sc
        #    1=off, 2=on
        for i in range( 1 , pose.total_residue() + 1 ):
            data[i-1] = 11 + int(movemap.get_bb(i))*10 + int(movemap.get_chi(i))
        self.send_any( 'mm1.bz2 ' , pose , data , 2 )

    #...this one is weird, colors pdb by foldtree
    def send_foldtree( self , pose , foldtree='' ):
        # if not sent, use poses fold_tree
        if not foldtree:
            foldtree = pose.fold_tree()
        data = [0]*pose.total_residue()
        # remove jump data for identifying later
        njump = foldtree.num_jump()
        starts = [0]*njump
        stops = [0]*njump
        in_jumps = []
        # list of start and stop edges
        for x in range( 0 , njump ):
            ed = foldtree.jump_edge(x+1)
            s1 = ed.start()
            s2 = ed.stop()
            if s1 < s2:
                starts[x] = s1
                stops[x] = s2
            else:
                starts[x] = s2
                stops[x] = s1
        # each residue is either: in loop, edge, cutpoint, or jump_point
        #    entry and exit from jumps is tricky
        for i in range( 1 , pose.total_residue() + 1 ):
            # keeps the iden relative to the jump, not entry
            if foldtree.is_jump_point(i):
                data[i-1] = 1
                for j in range( 0 , len(starts) ):
                    if i == starts[j]:
                        in_jumps.append(j+1)
                    elif i == stops[j]:
                        in_jumps.remove(j+1)
            elif foldtree.is_cutpoint(i):
                data[i-1] = 0
            # color varies for in jump, up to 7 jumps can easily be viewed
            #    edit later for more jumps
            elif in_jumps:
                data[i-1] = 2+max(in_jumps)
            else:
                data[i-1] = 2
        self.send_any( 'ft1.bz2 ' , pose , data , 1 )


    # send a foldtree visualization
    def view_foldtree( self , pose , foldtree = '' ):
        if not foldtree:
            foldtree = pose.fold_tree()
        njump = foldtree.num_jump()
        to_send = str(njump).rjust(2)    # 2 spaces
        for j in range( 1 , njump + 1 ):
            jump = foldtree.jump_edge(j)
            to_send = to_send + str(jump.start()).rjust(3) + str(jump.stop()).rjust(3) + str(foldtree.cutpoint_by_jump(j)).rjust(3)    # 3 spaces each, 9 per jump
        for i in range( 1 , pose.total_residue() + 1 ):
            if foldtree.is_cutpoint(i):
                to_send = to_send + '0'
            else:
                to_send = to_send + str(pose.chain(i))
            #to_send = to_send + str(int(foldtree.is_cutpoint(i)))    # 1 space per residue, map to rosetta resi
        name = self.getPoseName(pose)
        message = 'ftd.bz2 ' + chr(self.keep_history) + chr(len(name)) + name + bz2.compress(to_send)
        self.link.sendMessage(message)


    # send hbonds, labels absent/difficult, uses pymol "dist" function
    def send_hbonds( self , pose ):
        if not pose.energies().energies_updated():
            raise IOError('PyMOL_Mover::send_hbonds: Energy is not updated, please score the pose first!')
        hbset = rosetta.core.scoring.hbonds.HBondSet()
        rosetta.core.scoring.hbonds.fill_hbond_set(pose,False,hbset) # why False? what does this flag do
        # also, there are 4 additional flags! all default False but what do these do?
        to_send = str(hbset.nhbonds()).rjust(5)
        info = pose.pdb_info()
        energies = [0.]*hbset.nhbonds()
        # for min, max
        for j in range(0,hbset.nhbonds()):
            energies[j] = hbset.hbond(j+1).energy()
        maxe = max(energies)
        mine = min(energies)
        for i in range(1,hbset.nhbonds()+1):
            hb = hbset.hbond(i)
            accatm = pose.residue(hb.acc_res()).atom_name(hb.acc_atm())
            donatm = pose.residue(hb.don_res()).atom_name(hb.don_hatm())
            # each hbond 6+4+6+4+2 = 22 chars
            if info and info.nres():
                to_send = to_send + (str(info.number(hb.acc_res()))+info.icode(hb.acc_res())+info.chain(hb.acc_res())).rjust(6) + accatm + (str(info.number(hb.don_res()))+info.icode(hb.don_res())+info.chain(hb.don_res())).rjust(6) + donatm + ('%02x' % int((energies[i-1]-mine)*255./(maxe-mine)))
            else:
                to_send = to_send + str(hb.acc_res()).rjust(5) + accatm + str(hb.don_res()).rjust(5) + donatm + str(hb.energy).rjust(5)
        name = self.getPoseName(pose)
        message = 'hbd.bz2 ' + chr(self.keep_history) + chr(len(name)) + name + bz2.compress(to_send)
        self.link.sendMessage(message)


    # send a list of x and y to plot in pymol
    # name of object, connect as a color str or empty '' for no color/connect
    #    scale turns aces on or off, axis_color colors the axes,
    #    num of internal numberings on the scale
    def send_graph( self , name , connect , x_array , y_array , z_array = False , scale =True , axis_color = 'blue' , num = 0 ):
    # sends data 'name' of 'i_array' to PyMOL,
    #    connects based on 'connect' string, additional axis options
        if not z_array:
            z_array = [0]*len(x_array)
        if not len(x_array) == len(y_array) == len(z_array):
            raise IOError('FAIL: arrays are not the same length!')
        to_send = connect.rjust(7) + str(int(scale)) + axis_color.rjust(7) + str(num).rjust(6)

        for i in range( 0 , len(x_array) ):
            to_send += str(x_array[i])[:9].rjust(9)
            to_send += str(y_array[i])[:9].rjust(9)
            to_send += str(z_array[i])[:9].rjust(9)
        message = 'grp1.bz2' + chr(self.keep_history) + chr(len(name)) + name + bz2.compress(to_send)
        self.link.sendMessage(message)


    # send a point to add to a plot in pymol
    # options same as send_graph, banner is the point title
    #    rescale determines if axes are reset due to new point
    def send_point( self , name , connect , x , y , z = False , rescale = True , scale = True , axis_color = 'blue' , num = 0 , banner = '' ):
    # sends point to 'name' data in PyMOL
    #    connects based on 'connect' string, additional axis options
        if not z:
            z = 0
        to_send = str(x)[:9].rjust(9) + str(y)[:9].rjust(9) + str(z)[:9].rjust(9) + connect.rjust(7) + str(int(rescale)) + str(int(scale)) + axis_color.rjust(7) + str(num).rjust(6) + banner
        message = 'pnt.bz2 ' + chr(self.keep_history) + chr(len(name)) + name + bz2.compress(to_send)
        self.link.sendMessage(message)


    # displays sigs number of characters for each energy, labels on CA
    def label_energy( self , input_pose , energy_type='total_score' , sigs = 6 ):
        if not input_pose.energies().energies_updated():
            raise IOError('PyMOL_Mover::send_energy: Energy is not updated, please score the pose first!')

        # fix this later, see send_energy above!
        output = [0.]*input_pose.total_residue()
        if not energy_type == 'total_score':
            score_type = score_type_from_name(energy_type)
        for i in range( 0 , len(output) ):
            if energy_type == 'total_score':
                output[i] = input_pose.energies().residue_total_energy(i+1)
            else:
                output[i] = input_pose.energies().residue_total_energies(i+1)[score_type]
        self.send_any( 'lbE1.bz2' , input_pose , output , sigs )



###############################
################################################################################


'''
    def _output_energies( input_pose ):#, main_vars ):
        # could be implemented other ways, obtains a string of the pose name
        #pose_name = find_variable_name( input_pose , main_vars )[0]
        energies , weights = energies_to_list( input_pose )
        pymol = PyMOL_Mover()
        pymol.apply(input_pose)
        for i in range( 0 , len(energies) ):
            ##### KEY PORTION
            pymol.sendEnergies( input_pose , weights[i] , energies[i] )



    # unused!
    def energies_to_list( input_pose ):
        weights_obj = input_pose.energies().weights()
        weights = ['total_energy']
        for i in range( 1 , rosetta.core.scoring.end_of_score_type_enumeration + 1 ):
            score_name = name_from_score_type( ScoreType(i) )
            score_type = score_type_from_name(score_name)
            if not weights_obj[score_type] == 0.0:
                weights.append(score_name)
        output = []
        for j in range( 0 , len(weights)):
            if not weights[j] == 'total_energy':
                score_type = score_type_from_name(weights[j])
            temp = [0.]*input_pose.total_residue()
            output.append(temp)
            for k in range( 0 , input_pose.total_residue() ):
                if not weights[j] == 'total_energy':
                    output[j][k] = input_pose.energies().residue_total_energies(k+1)[score_type]
                else:
                    output[j][k] = input_pose.energies().residue_total_energy(k+1)
        return output , weights
    '''
# --------------------------------------------------------------------------------------------------
class PyMOL_Observer(core.pose.PosePyObserver):
    ''' Responds to general events (changes of geometry and energies) to pose and sends updates to
        pymol.
    '''
    def __init__(self, keep_history=False):
        rosetta.core.pose.PosePyObserver.__init__(self)
        self.pymol = rosetta.PyMOL_Mover(keep_history=keep_history)


    def generalEvent(self, event):
        #print 'PyMOL_Observer...'
        #print 'PyMOL_Observer:generalEvent', event.pose
        self.pymol.apply(event.getPose())

# --------------------------------------------------------------------------------------------------
# Energy method decorator, simplify creation of users EnergyMethods

_mem_EnergyMethods_ = []
_mem_EnergyCreators_ = []

_ScoreTypesRegestryByType_ = [
    CD(base=core.scoring.methods.ContextIndependentTwoBodyEnergy, first=PyRosettaTwoBodyContextIndepenedentEnergy_first, last=PyRosettaTwoBodyContextIndepenedentEnergy_last, methods={}),
    CD(base=core.scoring.methods.ContextDependentTwoBodyEnergy, first=PyRosettaTwoBodyContextDependentEnergy_first, last=PyRosettaTwoBodyContextDependentEnergy_last, methods={}),
    CD(base=None, first=PyRosettaEnergy_first, last=PyRosettaEnergy_last, methods={}),
]

ScoreTypesRegestry = {}

def defineEnergyMethodCreator(class_, scoreType):
    class Abstract_EnergyMethodCreator(rosetta.core.scoring.methods.EnergyMethodCreator):
        def __init__(self):
            rosetta.core.scoring.methods.EnergyMethodCreator.__init__(self)

        def create_energy_method(self, energy_method_options):
            e = self.EnergyMethodClass()
            _mem_EnergyMethods_.append(e)
            return e

        def score_types_for_method(self):
            sts = rosetta.utility.vector1_ScoreType();  sts.append( self.scoreType )
            return sts

    class_name = class_.__name__ + '_Creator'
    new_class = type(class_name, (Abstract_EnergyMethodCreator,), {'EnergyMethodClass' : class_, 'scoreType' : rosetta.core.scoring.ScoreType(scoreType)})
    #globals()[class_name ] = new_class

    return new_class


class EnergyMethod:
    def __init__(self, scoreName=None, scoreType=None, version=1):
        self.scoreName = scoreName
        self.scoreType = scoreType
        self.version = version

    def __call__(self, original_class):
        self.scoreName = self.scoreName or original_class.__name__
        if self.scoreType is None:  # trying to automatically determent first avaliable scoreType
            for s in _ScoreTypesRegestryByType_:
                if (s.base is None)  or  issubclass(original_class, s.base):
                    self.scoreType = max( s.methods.keys() or [s.first-1] ) + 1
                    if self.scoreType > s.last: raise Exception('Can not find free ScoreType to create %s! (looking in range [%s, %s])' % (self.scoreName, s.first, s.last) )
                    s.methods[self.scoreType] = self.scoreName
                    ScoreTypesRegestry[self.scoreType] = self.scoreName
                    break

        def clone(self_): return type(self_)()
        def f_version(self_): return self.version
        def indicate_required_context_graphs(self_, v): pass

        creator = defineEnergyMethodCreator(original_class, self.scoreType)

        if 'clone' not in original_class.__dict__: original_class.clone = clone
        if 'version' not in original_class.__dict__: original_class.version = f_version
        if 'indicate_required_context_graphs' not in original_class.__dict__: original_class.indicate_required_context_graphs = indicate_required_context_graphs

        original_class.creator = creator
        original_class.scoreType = rosetta.core.scoring.ScoreType(self.scoreType)

        _mem_EnergyCreators_.append( creator() )
        rosetta.core.scoring.methods.PyEnergyMethodRegistrator( _mem_EnergyCreators_[-1] )

        return original_class

# --------------------------------------------------------------------------------------------------

# By Michael Pacella
def etable_atom_pair_energies(atom1, atom2, sfxn):
    ''' Usage: lj_atr, lj_rep, solv=etable_atom_pair_energies(atom1, atom2, sfxn)
	Description: given a pair of atoms and scorefunction, use the precomputed
	'etable' to return LJ attractive, LJ repulsive, and LK solvation energies
    '''

    score_manager=core.scoring.ScoringManager.get_instance()
    etable_ptr=score_manager.etable(sfxn.energy_method_options().etable_type())
    etable=etable_ptr.get()
    etable_energy=core.scoring.etable.EtableEnergy(etable,sfxn.energy_method_options())

	#constructing AtomPairEnergy container to hold computed energies
    ape=core.scoring.etable.AtomPairEnergy()

	#setting all energies in the AtomPairEnergy to zero prior to calculation
    ape.attractive, ape.bead_bead_interaction, ape.repulsive, ape.solvation=0.0, 0.0, 0.0, 0.0

	#calculating distance squared and setting it in the AtomPairEnergy
    ape.distance_squared = atom1.xyz().distance_squared(atom2.xyz())

	#evaluate energies from pre-calculated etable, using a weight of 1.0
	#in order to match the raw energies from eval_ci_2b
    etable_energy.atom_pair_energy(atom1,atom2,1.0,ape)

	#calculating atom_atom scores
    lj_atr=ape.attractive
    lj_rep=ape.repulsive
    solv=ape.solvation

    return lj_atr, lj_rep, solv





# --------------------------------------------------------------------------------------------------

class PyRosettaException(Exception):
    #def __init__(self): pass
    def __str__(self): return 'PyRosettaException'


class PythonPyExitCallback(utility.PyExitCallback):
    def exit_callback(self):
        raise PyRosettaException()

    def __init__(self): utility.PyExitCallback.__init__(self)


__global_PythonPyExitCallback__ = None

def init(*args, **kargs):
    global __global_PythonPyExitCallback__

    __global_PythonPyExitCallback__ = PythonPyExitCallback()
    utility.PyExitCallback.set_PyExitCallBack(__global_PythonPyExitCallback__)

    #utility.set_pyexit_callback()  # make sure that all mini sys exit just generate exceptions


    #if not args: args = ["app"
    #                  ,"-database"
    #                  ,os.path.join( os.path.expanduser("~"), "rosetta_database")
    #                  ]
    #if not args: args = ["app", "-database", "rosetta_database", "-ex1", "-ex2aro"]

    # Figure out database dir...
    if os.path.isdir('rosetta_database'):
        database = os.path.abspath('rosetta_database')
        print 'Found rosetta_database at %s, using it...' % database

    elif 'PYROSETTA_DATABASE' in os.environ:
        database = os.path.abspath( os.environ['PYROSETTA_DATABASE'] )
        print 'PYROSETTA_DATABASE environment variable was set to: %s... using it...' % database

    elif os.path.isdir(os.environ['HOME'] + '/rosetta_database'):
        database = os.path.abspath(os.environ['HOME'] + '/rosetta_database')
        print 'Found rosetta_database at home folder, ie: %s, using it...' % database

    elif sys.platform == "cygwin" and os.path.isdir('/rosetta_database'):
        database = os.path.abspath('/rosetta_database')
        print 'Found rosetta_database at root folder, ie: %s, using it...' % database

    else:
        print 'Could not find rosetta_database! Check your paths or set PyRosetta environment vars. Exiting...'
        sys.exit(1)


    if not args: args = ["app", "-database", database, "-ex1", "-ex2aro"]

    args.extend( kargs.get('extra_options', '').split(' ') )

    v = utility.vector1_string()
    #v = utility.vector1_string('--database %s' % database)
    #v = Vector1('--database %s' % database)
    v.extend(args)
    #v.append('--database')
    #v.append(database)

    print version()
    protocols.init.init(v)


def version():
    return "PyRosetta 2.011 [r%s] retrieved from: %s" % (rosetta.core.minirosetta_svn_version(), rosetta.core.minirosetta_svn_url()) + \
    '\n(C) Copyright Rosetta Commons Member Institutions.\nCreated in JHU by Sergey Lyskov and PyRosetta Team.\n'


X11Colors = {
    'black': (0, 0, 0, 0),
    'AntiqueWhite': (1, 250, 235, 215),
    'BlanchedAlmond': (2, 255, 235, 205),
    'BlueViolet': (3, 138, 43, 226),
    'CadetBlue': (4, 95, 158, 160),
    'CornflowerBlue': (5, 100, 149, 237),
    'DarkBlue': (6, 0, 0, 139),
    'DarkCyan': (7, 0, 139, 139),
    'DarkGoldenrod': (8, 184, 134, 11),
    'DarkGray': (9, 169, 169, 169),
    'DarkGreen': (10, 0, 100, 0),
    'DarkGrey': (11, 169, 169, 169),
    'DarkKhaki': (12, 189, 183, 107),
    'DarkMagenta': (13, 139, 0, 139),
    'DarkOliveGreen': (14, 85, 107, 47),
    'DarkOrange': (15, 255, 140, 0),
    'DarkOrchid': (16, 153, 50, 204),
    'DarkRed': (17, 139, 0, 0),
    'DarkSalmon': (18, 233, 150, 122),
    'DarkSeaGreen': (19, 143, 188, 143),
    'DarkSlateBlue': (20, 72, 61, 139),
    'DarkSlateGray': (21, 47, 79, 79),
    'DarkSlateGrey': (22, 47, 79, 79),
    'DarkTurquoise': (23, 0, 206, 209),
    'DarkViolet': (24, 148, 0, 211),
    'DebianRed': (25, 215, 7, 81),
    'DeepPink': (26, 255, 20, 147),
    'DeepSkyBlue': (27, 0, 191, 255),
    'DimGray': (28, 105, 105, 105),
    'DimGrey': (29, 105, 105, 105),
    'DodgerBlue': (30, 30, 144, 255),
    'FloralWhite': (31, 255, 250, 240),
    'ForestGreen': (32, 34, 139, 34),
    'GhostWhite': (33, 248, 248, 255),
    'GreenYellow': (34, 173, 255, 47),
    'HotPink': (35, 255, 105, 180),
    'IndianRed': (36, 205, 92, 92),
    'LavenderBlush': (37, 255, 240, 245),
    'LawnGreen': (38, 124, 252, 0),
    'LemonChiffon': (39, 255, 250, 205),
    'LightBlue': (40, 173, 216, 230),
    'LightCoral': (41, 240, 128, 128),
    'LightCyan': (42, 224, 255, 255),
    'LightGoldenrod': (43, 238, 221, 130),
    'LightGoldenrodYellow': (44, 250, 250, 210),
    'LightGray': (45, 211, 211, 211),
    'LightGreen': (46, 144, 238, 144),
    'LightGrey': (47, 211, 211, 211),
    'LightPink': (48, 255, 182, 193),
    'LightSalmon': (49, 255, 160, 122),
    'LightSeaGreen': (50, 32, 178, 170),
    'LightSkyBlue': (51, 135, 206, 250),
    'LightSlateBlue': (52, 132, 112, 255),
    'LightSlateGray': (53, 119, 136, 153),
    'LightSlateGrey': (54, 119, 136, 153),
    'LightSteelBlue': (55, 176, 196, 222),
    'LightYellow': (56, 255, 255, 224),
    'LimeGreen': (57, 50, 205, 50),
    'MediumAquamarine': (58, 102, 205, 170),
    'MediumBlue': (59, 0, 0, 205),
    'MediumOrchid': (60, 186, 85, 211),
    'MediumPurple': (61, 147, 112, 219),
    'MediumSeaGreen': (62, 60, 179, 113),
    'MediumSlateBlue': (63, 123, 104, 238),
    'MediumSpringGreen': (64, 0, 250, 154),
    'MediumTurquoise': (65, 72, 209, 204),
    'MediumVioletRed': (66, 199, 21, 133),
    'MidnightBlue': (67, 25, 25, 112),
    'MintCream': (68, 245, 255, 250),
    'MistyRose': (69, 255, 228, 225),
    'NavajoWhite': (70, 255, 222, 173),
    'NavyBlue': (71, 0, 0, 128),
    'OldLace': (72, 253, 245, 230),
    'OliveDrab': (73, 107, 142, 35),
    'OrangeRed': (74, 255, 69, 0),
    'PaleGoldenrod': (75, 238, 232, 170),
    'PaleGreen': (76, 152, 251, 152),
    'PaleTurquoise': (77, 175, 238, 238),
    'PaleVioletRed': (78, 219, 112, 147),
    'PapayaWhip': (79, 255, 239, 213),
    'PeachPuff': (80, 255, 218, 185),
    'PowderBlue': (81, 176, 224, 230),
    'RosyBrown': (82, 188, 143, 143),
    'RoyalBlue': (83, 65, 105, 225),
    'SaddleBrown': (84, 139, 69, 19),
    'SandyBrown': (85, 244, 164, 96),
    'SeaGreen': (86, 46, 139, 87),
    'SkyBlue': (87, 135, 206, 235),
    'SlateBlue': (88, 106, 90, 205),
    'SlateGray': (89, 112, 128, 144),
    'SlateGrey': (90, 112, 128, 144),
    'SpringGreen': (91, 0, 255, 127),
    'SteelBlue': (92, 70, 130, 180),
    'VioletRed': (93, 208, 32, 144),
    'WhiteSmoke': (94, 245, 245, 245),
    'YellowGreen': (95, 154, 205, 50),
    'aquamarine': (96, 127, 255, 212),
    'azure': (97, 240, 255, 255),
    'beige': (98, 245, 245, 220),
    'bisque': (99, 255, 228, 196),
    'AliceBlue': (100, 240, 248, 255),
    'blue': (101, 0, 0, 255),
    'blue1': (102, 0, 0, 255),
    'blue2': (103, 0, 0, 238),
    'blue3': (104, 0, 0, 205),
    'blue4': (105, 0, 0, 139),
    'brown': (106, 165, 42, 42),
    'burlywood': (107, 222, 184, 135),
    'chartreuse': (108, 127, 255, 0),
    'chocolate': (109, 210, 105, 30),
    'coral': (110, 255, 127, 80),
    'cornsilk': (111, 255, 248, 220),
    'cyan': (112, 0, 255, 255),
    'firebrick': (113, 178, 34, 34),
    'gainsboro': (114, 220, 220, 220),
    'gold': (115, 255, 215, 0),
    'goldenrod': (116, 218, 165, 32),
    'gray': (117, 190, 190, 190),
    'gray0': (118, 0, 0, 0),
    'gray10': (119, 26, 26, 26),
    'gray100': (120, 255, 255, 255),
    'gray20': (121, 51, 51, 51),
    'gray30': (122, 77, 77, 77),
    'gray40': (123, 102, 102, 102),
    'gray50': (124, 127, 127, 127),
    'gray60': (125, 153, 153, 153),
    'gray70': (126, 179, 179, 179),
    'gray80': (127, 204, 204, 204),
    'gray90': (128, 229, 229, 229),
    'green': (129, 0, 255, 0),
    'green1': (130, 0, 255, 0),
    'green2': (131, 0, 238, 0),
    'green3': (132, 0, 205, 0),
    'green4': (133, 0, 139, 0),
    'honeydew': (134, 240, 255, 240),
    'ivory': (135, 255, 255, 240),
    'khaki': (136, 240, 230, 140),
    'lavender': (137, 230, 230, 250),
    'linen': (138, 250, 240, 230),
    'magenta': (139, 255, 0, 255),
    'maroon': (140, 176, 48, 96),
    'moccasin': (141, 255, 228, 181),
    'navy': (142, 0, 0, 128),
    'orange': (143, 255, 165, 0),
    'orchid': (144, 218, 112, 214),
    'peru': (145, 205, 133, 63),
    'pink': (146, 255, 192, 203),
    'plum': (147, 221, 160, 221),
    'purple': (148, 160, 32, 240),
    'red': (149, 255, 0, 0),
    'red1': (150, 255, 0, 0),
    'red2': (151, 238, 0, 0),
    'red3': (152, 205, 0, 0),
    'red4': (153, 139, 0, 0),
    'salmon': (154, 250, 128, 114),
    'seashell': (155, 255, 245, 238),
    'sienna': (156, 160, 82, 45),
    'snow': (157, 255, 250, 250),
    'snow1': (158, 255, 250, 250),
    'snow2': (159, 238, 233, 233),
    'snow3': (160, 205, 201, 201),
    'snow4': (161, 139, 137, 137),
    'tan': (162, 210, 180, 140),
    'thistle': (163, 216, 191, 216),
    'tomato': (164, 255, 99, 71),
    'turquoise': (165, 64, 224, 208),
    'violet': (166, 238, 130, 238),
    'wheat': (167, 245, 222, 179),
    'white': (168, 255, 255, 255),
    'yellow': (169, 255, 255, 0),
}
