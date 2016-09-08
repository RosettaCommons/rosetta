/*
	test_D_scoring.cc
	A pilot app to test the Rama scorefunction's behaviour with terminal residues
	in backbone-cyclized peptides.

	File history:
		--Created 5 February 2014.
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
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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

// This should be removed.
#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

//OPT_KEY( Integer, peptidesize)
//OPT_KEY( Real, energycutoff)
//OPT_KEY( Integer, relaxrnds)
//OPT_KEY( IntegerVector, dpositions)

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//utility::vector1 <core::Size> emptyintlist;

	//NEW_OPT( peptidesize, "The number of residues in the peptide, including terminal cysteines.  Default 9.", 9);
	//NEW_OPT( energycutoff, "The energy above which structures are discarded.  Only used if a value is specified by the user.", 1000.0);
	//NEW_OPT( relaxrnds, "The number of rounds of relaxation with the fastrelax mover.  Default 3.", 3);
	//NEW_OPT( dpositions, "The positions at which a D-alanine should be placed.  Default empty list.", emptyintlist);
}

void constrain_loop_ends (
	core::pose::Pose &mypose,
	const core::Size resC, //The last residue of the loop (provides C)
	const core::Size resN, //The anchor residue (provides N)
	core::scoring::ScoreFunctionOP sfxn_constrained
) {
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	using namespace core::id;
	mypose.remove_constraints();
	sfxn_constrained->set_weight(atom_pair_constraint, 1.0);
	sfxn_constrained->set_weight(dihedral_constraint, 1.0);
	sfxn_constrained->set_weight(angle_constraint, 1.0);

	const std::string Hstring = (mypose.residue(resN).has("H")?"H":"CD");

	//Make the linkage between the N- and C-termini:
	mypose.conformation().declare_chemical_bond(resC, "C", resN, "N"); //Declare a chemical bond between the N and C termini.
	{ 
		//Peptide bond length constraint:
		mypose.add_constraint (
			new AtomPairConstraint (
				AtomID( mypose.residue(resC).atom_index("C") , resC ) ,
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				new HarmonicFunc( 1.3288, 0.01)
			)
		);

		//Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		mypose.add_constraint (
			new DihedralConstraint (
				AtomID( mypose.residue(resC).atom_index("O") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index(Hstring) , resN) ,
				new CircularHarmonicFunc( PI, 0.02)
			)
		);
		mypose.add_constraint (
			new DihedralConstraint (
				AtomID( mypose.residue(resC).atom_index("CA") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index("CA") , resN) ,
				new CircularHarmonicFunc( PI, 0.02)
			)
		);

		//Peptide bond angle constraints:
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index("CA") , resN) ,
				new CircularHarmonicFunc( CNCa_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index(Hstring) , resN) ,
				new CircularHarmonicFunc( CNH_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("CA") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				new CircularHarmonicFunc( CaCN_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("O") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				new CircularHarmonicFunc( OCN_ANGLE/180.0*PI, 0.02)
			)
		);

	}
	
}


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


	printf("Starting test_cyclic_rama.cc\n");
	printf("Pilot app created 5 February 2014 by Vikram K. Mulligan, Baker Laboratory.\n");
	register_options(); //Parse the input flags and set global variables appropriately.
	devel::init(argc, argv); //Initialize Rosetta

	if(!option[in::file::s].user()) {
		printf("Error!  The user must specify an input file with the -in:file:s flag.  Crashing.\n"); fflush(stdout);
		exit(1);
	}

	//Default scorefunction:
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.

	//Movers:
	protocols::relax::FastRelax frlx(sfxn, 5);
	//protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn);
	// AMW: VKM prefers dfpmin for short peptides.
	core::optimization::MinimizerOptions minoptions("dfpmin_armijo_nonmonotone", 0.000000001, true, false, false);
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_bb(true); //Allow backbone motion
	mm->set_chi(true); //Allow side-chain motion
	mm->set_jump(false); //Don't bother with "jumps" (relative motion of different molecules relative one another)
	//core::optimization::AtomTreeMinimizer minimizer; //A torsion-space minimizer
	//core::optimization::CartesianMinimizer cminimizer; //A cartesian-space minimizer

	core::pose::Pose mypose;
	core::import_pose::pose_from_file (mypose, option[in::file::s]()[1]);
	const core::Size nres = mypose.size();
	printf("Import complete.\n"); fflush(stdout);

	core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(mypose, nres);
	constrain_loop_ends(mypose, nres, 1, sfxn);
	printf("Setup complete.\n"); fflush(stdout);

	/*char outfile[256];
	core::pose::Pose mypose2=mypose;
	for(core::Size itors=0; itors<36; itors++) {
		mypose2.conformation().set_torsion( core::id::TorsionID(1, core::id::BB, 2), (Real)itors*10.0);
		mypose2.update_residue_neighbors();
		sprintf(outfile, "out%02lu.pdb", itors);
		mypose2.dump_pdb(outfile);
	}*/

	printf ("\nResidue\tTerm\tLower\tUpper\n");
	for(core::Size ir=1; ir<=nres; ir++) {
		printf("%lu\t", ir);
		if(mypose.residue(ir).is_terminus()) printf("true\t"); else printf("false\t");
		printf("%lu\t%lu\n", mypose.residue(ir).residue_connection_partner(1), mypose.residue(ir).residue_connection_partner(2) );
	}
	printf("\n"); fflush(stdout);

	mypose.dump_scored_pdb("initial.pdb", (*sfxn));
	printf("Wrote initial.pdb.  Relaxing...\n"); fflush(stdout);

	frlx.apply(mypose);

	mypose.dump_scored_pdb("final.pdb", (*sfxn));
	printf("Wrote final.pdb.\n"); fflush(stdout);

	printf("***JOB COMPLETED.***\n");
	return 0;
}

