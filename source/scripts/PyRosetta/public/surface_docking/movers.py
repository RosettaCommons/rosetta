# Mover classes for surface docking protocol
# Emily Koo

# Import python modules
import math
from random import randint
from os.path import exists
import glob
from append_data import *

# Import Rosetta modules
from rosetta import *
from rosetta.protocols.rigid import *
import rosetta.core.scoring.solid_surface
import rosetta.core.scoring.constraints

# Import custom modules
from surf_param import *
from constraints import *

#=========================== Setup movers ==============================
class Abinitio:
    def __init__(self, PDB, frag3, frag9):
        print "Creating Ab initio instance..."
        print ">> Loading fragment files..."
        print PDB
        self.fragset_9mer = core.fragment.ConstantLengthFragSet(9, frag9)
        self.fragset_3mer = core.fragment.ConstantLengthFragSet(3, frag3)
        
        print ">> Setting up movemap..."
        self.movemap = MoveMap()
        self.movemap.set_bb(True)
        self.movemap.set_chi(False)
        self.mover_9mer = ClassicFragmentMover(self.fragset_9mer, self.movemap) 
        self.mover_3mer = ClassicFragmentMover(self.fragset_3mer, self.movemap)
        
        print ">> Setting up ClassicAbinitio mover..."
        self.abinitio = ClassicAbinitio(self.fragset_3mer, self.fragset_9mer, self.movemap)
        # Set rg to 0 to get more extended conformers
        self.abinitio.set_score_weight(rg, 0.0)

    def apply(self, pose):
 
        # Score pose
        self.score_low = create_score_function("RS_centroid.wts")
        self.score_low(pose)
        self.score_low.show(pose)
        print ".............Abinitio starting................"
        self.abinitio.init(pose)
        self.abinitio.apply(pose)
        print "................Abinitio done................."
        print "\n"
        
    def info(self):
        return self.__doc__
            
    def get_name(self):
        return self.__class__.__name__
        
class CentroidRelax:
    # Set-up default parameters
    kT = 0.5
    inner_cycles = 6
    outer_cycles = 43   # will be changed
    n_moves = 6        # number of small/shear moves
    max_angle = 30

    def __init__(self):
        print "Creating centroid relax instance..."
        print ">> Loading score function..."
        self.score_low = create_score_function("RS_centroid.wts")
        self.smallmin_type = "linmin"
        self.shearmin_type = "lbfgs_armijo"
        
        # Construct movemap
        self.movemap = rosetta.MoveMap()
        self.movemap.set_bb(True)
    
        self._setupMovers()
        
    def _setupMovers(self):
        print ">> Setting up small and shear movers..."
        # Small trial
        self.small_mover = rosetta.SmallMover(self.movemap, self.kT, self.n_moves/2)
        self.small_mover.angle_max(self.max_angle)
        self.smallmin_mover = rosetta.MinMover(self.movemap, self.score_low, self.smallmin_type, 0.1, True)
        self.small_seq = rosetta.SequenceMover()
        self.small_seq.add_mover(self.small_mover)
        self.small_seq.add_mover(self.smallmin_mover)
        
        # Shear trial
        self.shear_mover = rosetta.ShearMover(self.movemap, self.kT, self.n_moves)
        self.shear_mover.angle_max(self.max_angle)
        self.shearmin_mover = rosetta.MinMover(self.movemap, self.score_low,self.shearmin_type, 0.01, True)
        self.shear_seq = rosetta.SequenceMover()
        self.shear_seq.add_mover(self.shear_mover)
        self.shear_seq.add_mover(self.shearmin_mover)

    def _setupTrial(self, pose):
        print ">> Setting up trial movers..."
        self.mc = MonteCarlo(pose, self.score_low, self.kT)
        self.small_trial = rosetta.TrialMover(self.small_seq, self.mc)
        self.shear_trial = rosetta.TrialMover(self.shear_seq, self.mc)
        
        self.small_rep = rosetta.RepeatMover(self.small_trial, self.inner_cycles)

    def set_nmoves(self, nmoves):
        self.n_moves = nmoves
        
    def set_max_angle(self, angle):
        # Changes  max angles for small and shear movers
        self.max_angle = angle
        
    def apply(self, pose):
    
        # Score pose
        self.score_low(pose)
        
        # Set up MonteCarlo and trial movers
        self._setupTrial(pose)
        
        # Set outer cycles = protein length
        self.outer_cycles = pose.total_residue()
        
        print ">> Starting centroid relaxation..."
        for outerLoop in range(self.outer_cycles):
            self.score_low.show(pose)
            self.mc.recover_low(pose)
            self.shear_trial.apply(pose)
            self.small_rep.apply(pose)
        
        self.mc.show_counters()
        self.mc.recover_low(pose)
        print ">> Finished centroid relaxation."
        
    def info(self):
        return self.__doc__
            
    def get_name(self):
        return self.__class__.__name__
        
class FullAtomRelax:
    # Set-up default parameters
    kT = 0.5
    n_moves = 6         # number of small/shear moves
    tolerance = 0.01   # min mover tolerance
    
    def __init__(self, score_high, score_pack, std_scorefxn, nosmallshear):
        print "Creating full atom relax instance..."
        print ">> Initializing variables..."
        
        # Score functions
        self.score_high = score_high
        self.score_pack = score_pack
        self.std_scorefxn = std_scorefxn
        
        self.nosmallshear = nosmallshear    
        # Initial min movers for small and shear
        self.smallmin_type = "lbfgs_armijo"
        self.shearmin_type = "lbfgs_armijo"
        
        # Initialize constraint files
        self.constraints = 'False' # constraints toggle. false by default
        self.sol_cst = None
        self.ads_cst = None
 
        # Other variables
        self.max_angle = 30
        self.outer_cycles = 5
        self.inner_cycles = 0
        self.outer_cycle = 0   # Current ramp cycle
        self.ref_cycles = 0
        self.curr_cycle = 0 # Current ref cycle
        self.protStart = 0
        self.totalRes = 0
        
        # Packer task
        self.task = None      
        
        # Movemaps
        # Both bb and chi flexible
        self.movemap = None
        # Only chi flexible
        self.chi_mm = None

        self.state = "ads" 
        self.sol_constraints = False
        self.ads_constraints = False
        self.name = "" #name of decoy

    def set_scorefxn(self, scorefxn):
        # Set scorefxn to specified scorefxn
        self.score_high = scorefxn
    
    def set_mc(self, pose):
        # Get input pose to set up mc and trial movers
        self.mc = MonteCarlo(pose, self.score_high, self.kT)

    
    def set_params(self, pose):
    
        # Set mc
        self.set_mc(pose)
        
        # Set movemap
        self.set_movemap(pose.total_residue(), pose.num_jump() + 1)                  
        
        # Set up basic movers
        self._setupMovers()  
        
        # Set up trial movers
        self._setupTrial(pose)   
        
        # Set up Dock MCM
        self.dockMCM = DockMCMProtocol(pose.num_jump(), self.score_high, self.score_pack) 
        
        # Set up surface energies to score only protein
        self._setupSurfe(pose)

    def set_nmoves(self, nmoves):
        self.n_moves = nmoves
        
    def set_max_angle(self, angle):
        # Changes max angle for shear mover
        self.max_angle = angle
        # Set max angle for small mover
        self.small_mover.angle_max('H',angle/3)
        self.small_mover.angle_max('E',angle/2)
        self.small_mover.angle_max('L',angle)

        
    def set_movemap(self, totalRes, protStart):  
        # Get user input
        self.protStart = protStart
        self.totalRes = totalRes
        
        # Update movemap
        self._setupMovemap()
        
    def loadSurf(self, file):
        print ">> Loading surface vectors..."
        pdb = load(file)
        self.SURFA0 = get_markers(pdb)[0]
        self.SURFA1 = get_markers(pdb)[1]
        self.SURFA2 = get_markers(pdb)[2]

    def loadConstraints(self, PDB):
        
        if self.sol_cst is None:
            # Atom pair and dihedral constraints
            self.sol_cst = load_constraints(PDB)[0]
        if self.ads_cst is None:
            self.ads_cst = load_constraints(PDB)[1]

    def set_constraints_weight(self, weight):
        self.score_high.set_weight(atom_pair_constraint, weight)
        self.score_high.set_weight(dihedral_constraint, weight)

    def _setupSurfe(self, pose):
        
        print ">> Setting up surface energy..."
        # Set Surface Energy to 0 - Get score only from protein
        surfe = core.scoring.solid_surface.SurfaceEnergies()
        surfe.set_total_residue(pose.total_residue())
        surfe.set_residue_range_not_surface(pose.num_jump() + 1, pose.total_residue())
        pose.set_new_energies_object(surfe)
        self.score_high(pose)

    def _setupMovemap(self):
        
        print ">> Setting up movemap..."
        # Construct movemap
        self.movemap = rosetta.MoveMap()
        self.movemap.set_bb(False)
        self.movemap.set_chi(False)

        # Only allow protein bb and chi to move
        for residue in range(self.protStart, self.totalRes + 1):
            self.movemap.set_bb(residue, True)
            self.movemap.set_chi(residue, True)
        
        # If prot only, no jumps
        if self.protStart > 1:
            # Allow jump between surface and protein to move
            self.movemap.set_jump(self.protStart - 1, True)
        
    def _setupMovers(self):
        print ">> Setting up small and shear movers..."
        # Small mover
        self.small_mover = rosetta.SmallMover(self.movemap, self.kT, self.n_moves)
        self.small_mover.angle_max(self.max_angle)
        self.smallmin_mover = rosetta.MinMover(self.movemap, self.score_high, self.smallmin_type, self.tolerance, True)
        
        # Shear mover
        self.shear_mover = rosetta.ShearMover(self.movemap, self.kT, self.n_moves/2)
        self.shear_mover.angle_max(self.max_angle)
        self.shearmin_mover = rosetta.MinMover(self.movemap, self.score_high, self.shearmin_type, self.tolerance, True)
    
        # Set up specific min movers
        self.min_dfp_mover = rosetta.MinMover(self.movemap, self.score_high,"lbfgs_armijo", self.tolerance, True)
        self.min_lin_mover = rosetta.MinMover(self.movemap, self.score_high,"linmin", self.tolerance, True)
        
    def _setupTrial(self, pose):
        print ">> Setting up trial movers..."
                
        # Check to see if movemap is set
        if self.movemap is None:
            self.set_movemap(pose.total_residue(), pose.num_jump() + 1)    
        
        # Small seq movers
        self.small_seq = rosetta.SequenceMover()
        self.small_seq.add_mover(self.small_mover)
        self.small_seq.add_mover(self.smallmin_mover)
        
        # Shear seq movers
        self.shear_seq = rosetta.SequenceMover()
        self.shear_seq.add_mover(self.shear_mover)
        self.shear_seq.add_mover(self.shearmin_mover)
        
        # Trial small/shear movers
        self.smallmin_trial = rosetta.TrialMover(self.small_seq, self.mc)
        self.shearmin_trial = rosetta.TrialMover(self.shear_seq, self.mc)

        # Trial min movers
        self.min_dfp_trial = rosetta.TrialMover(self.min_dfp_mover, self.mc)
        self.min_lin_trial = rosetta.TrialMover(self.min_lin_mover, self.mc)
        
        # Repeat mover for inner/outer cycles
        self.min_seq = rosetta.SequenceMover()
        self.min_seq.add_mover(self.smallmin_trial)
        self.min_seq.add_mover(self.shearmin_trial)

        self.inner_min_rep = rosetta.RepeatMover(self.min_seq, self.inner_cycles)
        
    def _setMinmovers(self, minmover):
        # Update min mover types
        print ">>> Updating min mover types..."
        self.smallmin_mover.min_type(minmover)
        self.shearmin_mover.min_type(minmover)
        
            
    def _slideProt(self, pose):
        print ">>> Move protein onto surface..."
        # Symmetry move to get protein in middle of surface
        surface_symm_move(pose, pose.num_jump() + 1, pose.total_residue(), self.SURFA0, self.SURFA1, self.SURFA2)
  
        # Slide protein into contact with surface
        self.slide_into_contact = FaDockingSlideIntoContact(pose.num_jump())
        self.slide_into_contact.apply(pose)

        #self.std_scorefxn.show(pose)

    def _moveProt(self, pose):
    
        print ">>> Moving protein away from surface..."
        # Moving apart protein and surface
        self.trans_mover = RigidBodyTransMover(pose, pose.num_jump())
        self.trans_mover.step_size(50)
        self.trans_mover.apply(pose)
        
    def _randMovers(self, pose):
     	print ">>> Randomizing protein position..."

        # Moving apart protein and surface
        self._moveProt(pose)
    
        # Spin about a random axis
        #self.spin_mover = RigidBodySpinMover(pose.num_jump())
        #self.spin_mover.apply(pose)
        self.random_mover = RigidBodyRandomizeMover(pose, pose.num_jump(), partner_downstream)
        self.random_mover.apply(pose) 
        #self.std_scorefxn.show(pose)
        # Slide protein to surface
        self._slideProt(pose)
        
    def applyConstraints(self, pose, state):
       
        print ">>> Reading constraints..."
        # Read in constraints        
        if state == "sol":
            if (self.sol_constraints is False):
                apply_constraints(pose, self.sol_cst)
                self.sol_constraints = True;

        
        elif state == "ads":
            if (self.ads_constraints is False):
                apply_constraints(pose, self.ads_cst)
                self.ads_constraints = True;
   
        # IMPORTANT! Or constraints will not be activated for other movers
        self.mc.reset_scorefxn(pose, self.score_high)
        #self.score_high.show(pose)  
        
    def _setupPackertask(self, pose):
        
        print ">> Setting Side-chain packer task..."
        self.task = standard_packer_task(pose)
        self.task.restrict_to_repacking()
        self.task.or_include_current(True)
        
        print ">> Setting repacking range..."
        for r in range(1, pose.total_residue()+1):
            if r < pose.num_jump() + 1:
                # Don't pack surface
                self.task.nonconst_residue_task(r).prevent_repacking()
            else:
                # Pack protein
                self.task.nonconst_residue_task(r).restrict_to_repacking()            
        
    def fullRepack(self, pose):
        # Full repack of protein only (PackRotamersMover)
        
        # Set up packer task if not done
        if self.task is None:
            self._setupPackertask(pose)
        
        #print "Before full repack", #self.score_high.show(pose)
        self.full_repack = PackRotamersMover(self.score_pack, self.task)
        self.full_repack.apply(pose)
        #print "After full repack", #self.score_high.show(pose)

        """
        # Followed by dfp min with only chi flexible
        if self.chi_mm is None:
            self._setChimap(pose)
        
        chi_minmover = rosetta.MinMover(self.chi_mm, self.score_pack, "lbfgs_armijo", self.tolerance, True)
        chi_minmover.apply(pose)
        print "After chi minmover", #self.score_high.show(pose)
        """
    def _setChimap(self, pose):
            
        self.chi_mm = MoveMap()
        self.chi_mm.set_bb(False)
        self.chi_mm.set_chi(False)
        for residue in range(pose.num_jump() + 1, pose.total_residue() + 1):
            self.chi_mm.set_chi(residue, True)

    def _rotamerTrialmin(self, pose):
        # Rotamer trial min mover
        #print "Before Rotamer Trial Min scorefxn", #self.score_high.show(pose)
        self.rt_min = RotamerTrialsMinMover(self.score_high, self.task)
        self.rt_min.apply(pose)
        #print "After Rotamer Trial Min scorefxn", #self.score_high.show(pose)
        
    def _dock(self, pose):
        # Docking protocol
        print ">>> Dock MCM protocol..."
        #print "Before:"
        #self.score_high.show(pose)
        self.dockMCM.apply(pose)   
        #print "After:"
        #self.score_high.show(pose)                    

    def _fullatomRelax(self, pose):
        
        print ">> Outer cycle = ", self.outer_cycle
        # Adsorbed cycles only
        if self.curr_cycle == 1 and self.outer_cycle == 1 and self.state == "ads":
            # Take protein and surface apart and 
            # randomly orient and dock protein during first adsorbed state cycle
            self._randMovers(pose)
            # Rigid body docking 
            self._dock(pose)
                
        if not self.nosmallshear:
            # Initial settings =================================
            self.mc.reset(pose)
            self.mc.reset_scorefxn(pose, self.score_high)
            self.score_high(pose)
        
            # Outer cycle refinements ============================
            self._setMinmovers("lbfgs_armijo")
            print ">>> Small/shear refinements outer cycle..."
            self.min_seq.apply(pose)
        
            # Inner cycle refinements ============================
            self._setMinmovers("linmin")
    
            print ">>> Small/shear refinements inner cycles..."
            self.inner_min_rep.apply(pose) # shear trial, small trial
        
            # Recover lowest scoring decoy
            self.mc.recover_low(pose)   
        
        # DockMCM 
        if self.state == "ads":
            # Rigid body docking protocol on surface
            self._dock(pose)
        
    def apply(self, pose):

        # Check to make sure movemap is right
        #if self.totalRes is not pose.total_residue() or self.protStart is not pose.num_jump() + 1:
        #    self.set_params(pose)
            
        # Do full atom relax
        self._fullatomRelax(pose)

        # Reset applied constraints toggles
        self.ads_constraints = False
        self.sol_constraints = False
        
    def info(self):
        return self.__doc__
            
    def get_name(self):
        return self.__class__.__name__
 
 
