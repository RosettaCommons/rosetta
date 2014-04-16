/*
	test_metalloproteins.cc
	A pilot app to test the automated setup of metalloproteins.

	File history:
		--Created 27 February 2014.
*/

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <devel/init.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
//#include <core/pose/metalloproteins/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <stdio.h>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/relax/FastRelax.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/kinematics/MoveMap.hh>
//#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>

#define PI 3.1415926535897932384626433832795

//OPT_KEY( Integer, peptidesize)
//OPT_KEY( Real, energycutoff)
//OPT_KEY( Integer, relaxrnds)
//OPT_KEY( IntegerVector, dpositions)

//Enzdes code from Florian:
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//utility::vector1 <core::Size> emptyintlist;

	//NEW_OPT( peptidesize, "The number of residues in the peptide, including terminal cysteines.  Default 9.", 9);
	//NEW_OPT( energycutoff, "The energy above which structures are discarded.  Only used if a value is specified by the user.", 1000.0);
	//NEW_OPT( relaxrnds, "The number of rounds of relaxation with the fastrelax mover.  Default 3.", 3);
	//NEW_OPT( dpositions, "The positions at which a D-alanine should be placed.  Default empty list.", emptyintlist);
}

/*void detect_covalent_bonds_to_metals(
	core::pose::Pose & pose
) {
	core::Size const nres = pose.n_residue(); //Number of residues in the pose.
 
	for(core::Size ir=1; ir<=nres; ir++) { //Loop through all residues, looking for metals
		if(pose.residue(ir).is_metal()) {
			printf("Residue %lu is a metal!\n", ir); fflush(stdout);

			utility::vector1 < core::id::AtomID > atomidlist = core::pose::metalloproteins::find_metalbinding_atoms( pose, ir,  1.05);
			core::pose::metalloproteins::add_covalent_linkages_to_metal( pose, ir, atomidlist );
	
		}
	}

	return;
} */

/*void set_up_metal_distance_constraints (
	core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn
) {
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	using namespace core::id;

	core::Size const nres = pose.n_residue(); //Number of residues in the pose.
	core::Real weightmultiplier = 1.0; //Constraint weight multiplier, for cases in which the atom_pair_constraint has already been turned on but is not equal to 1.0.
	if(sfxn->get_weight(atom_pair_constraint) < 1.0e-10) { //If this score term is turned off, turn it on.
		sfxn->set_weight(atom_pair_constraint, 1.0);
	} else { //Otherwise, if this score term is already turned on, set the weight multiplier appropriately.
		weightmultiplier = 1.0 / sfxn->get_weight(atom_pair_constraint); //Note -- we will take the square root of this later on.
	}

	for (core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues in the protein, looking for metals
		if(pose.residue(ir).is_metal()) { //When we find a metal, check what it's bound to and set up distance constraints accordingly.
			core::Size const n_resconn = pose.residue(ir).n_residue_connections(); //Number of residue connections

			//Make a vector of indices of virtual atoms in this residue:
			utility::vector1 < core::Size > virtindices;
			for(core::Size ia=1, ia_max=pose.residue(ir).natoms(); ia<=ia_max; ++ia) {
				if(pose.residue(ir).is_virtual(ia)) virtindices.push_back(ia);
			}

			for (core::Size jr=1; jr<=n_resconn; ++jr) { //Loop through all of the metal's residue connections
				if(jr > virtindices.size()) { //We'll use virtual atoms for the constraints.  This means that if we have more connections than virtual residues, we can't do this.
					std::string message = "Error!  The number of connections to a metal is greater than the number of virtual atoms in that metal's residuetype.  Unable to continue.";
					utility_exit_with_message(message);
				}

				core::chemical::ResConnID const & curconnection = pose.residue(ir).connect_map(jr); //The current residue connection
				core::Size const otherres = curconnection.resid(); //The other residue to which this one is connected
				core::Size const otherres_conid = curconnection.connid(); //The other residue's connection id for the bond to the metal.
				core::Size const otherres_atom = pose.residue(otherres).residue_connection(otherres_conid).atomno(); //The atom index of the atom to which this metal is bonded.

				AtomID virtID(virtindices[jr], ir);
				AtomID otherID(otherres_atom, otherres);

				pose.set_xyz(virtID, pose.xyz(otherID)); //Move this residue's virt to the other residue's metal-binding atom position.
				HarmonicFuncOP hfunc = new HarmonicFunc( 0.0, 0.1 / sqrt(weightmultiplier)); //Harmonic function for constraining the position.
				AtomPairConstraintOP pairconst = new AtomPairConstraint(virtID, otherID, hfunc); //Atom pair constraint holding the virt at the position of the metal-binding atom.
				pose.add_constraint(pairconst); //Add the constraint to the pose, and carry on.
			}
		}
	}

	return;
}*/

/*void set_up_metal_constraints (
	core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn
) {

	//set_up_metal_distance_constraints(pose, sfxn);
	core::pose::metalloproteins::auto_setup_all_metal_constraints(pose, sfxn);

	return;
}*/

/**********
  MAIN!!!
**********/
int main(int argc, char *argv[]) {
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options; 
	using namespace basic::options::OptionKeys; 
	using namespace core::scoring::constraints; 
	using namespace core::id;
	using namespace core::pose;
	using namespace core;
	using namespace chemical;
	using namespace conformation;

	//using namespace protocols::toolbox::match_enzdes_util; //Florian's stuff

	numeric::random::RandomGenerator RG( 923749 ); //Random generator and seed

	printf("Starting test_metalloproteins.cc\n");
	printf("Pilot app created 27 February 2014 by Vikram K. Mulligan, Baker Laboratory.\n");
	register_options(); //Parse the input flags and set global variables appropriately.
	devel::init(argc, argv); //Initialize Rosetta

	if(!option[in::file::s].user()) {
		printf("Error!  The user must specify an input file with the -in:file:s flag.  Crashing.\n"); fflush(stdout);
		exit(1);
	}

	//Default scorefunction:
	core::scoring::ScoreFunctionOP sfxn = core::scoring::getScoreFunction();
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.

	//Movers:
	protocols::relax::FastRelax frlx(sfxn, 1);
	//protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn);
	//core::optimization::MinimizerOptions minoptions("dfpmin_armijo_nonmonotone", 0.000000001, true, false, false);
	core::optimization::MinimizerOptions minoptions("dfpmin", 0.000000001, true, false, false);
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_bb(true); //Allow backbone motion
	mm->set_chi(true); //Allow side-chain motion
	mm->set_jump(true); //Allow "jumps" (relative motion of different molecules relative one another)
	//core::optimization::AtomTreeMinimizer minimizer; //A torsion-space minimizer
	//core::optimization::CartesianMinimizer cminimizer; //A cartesian-space minimizer

	core::pose::Pose mypose;
	core::import_pose::pose_from_pdb (mypose, option[in::file::s]()[1]);
	//const core::Size nres = mypose.n_residue();
	printf("Import complete.\n"); fflush(stdout);

	mypose.dump_scored_pdb("1_initial.pdb", (*sfxn)); printf("Wrote 1_initial.pdb.\n"); fflush(stdout);

	//sfxn->set_weight(atom_pair_constraint, 1.0);
	//sfxn->set_weight(angle_constraint, 1.0);

	//Set up covalent geometry:
	//printf("Detecting bonds to metals.\n"); fflush(stdout);
	//core::pose::metalloproteins::auto_setup_all_metal_bonds(mypose, 1.05, true);
	//detect_covalent_bonds_to_metals(mypose);

	//Set up constraints:
	//printf("Setting up metal constraints.\n"); fflush(stdout);
	//set_up_metal_constraints(mypose, sfxn);

	//mypose.dump_scored_pdb("2_working.pdb", (*sfxn) ); printf("Wrote 2_working.pdb.\n"); fflush(stdout);

	frlx.apply(mypose);

	mypose.dump_scored_pdb("2_final.pdb", (*sfxn)); printf("Wrote 2_final.pdb.\n"); fflush(stdout);

	printf("*****JOB COMPLETED.*****\n");

	return 0;
}

