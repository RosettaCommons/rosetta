/*
	test_D_scoring.cc
	A pilot app to test whether two mirror-image peptides are scored identically or not.
	Created by Vikram K. Mulligan, Baker Laboratory.

	File history:
		--Created 23 Jun 2013.
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
//#include <protocols/simple_moves/MinMover.hh>
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

	numeric::random::RandomGenerator RG( 923749 ); //Random generator and seed

	printf("Starting test_D_scoring.cc\n");
	printf("Pilot app created 23 Jun 2013 by Vikram K. Mulligan, Baker Laboratory.\n");
	register_options(); //Parse the input flags and set global variables appropriately.
	devel::init(argc, argv); //Initialize Rosetta

	//Default scorefunction:
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.

	//Movers:
	protocols::relax::FastRelax frlx(sfxn, 5);
	protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn);
	//core::optimization::MinimizerOptions minoptions("dfpmin_armijo_nonmonotone", 0.000000001, true, false, false);
	core::optimization::MinimizerOptions minoptions("dfpmin", 0.000000001, true, false, false);
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_bb(true); //Allow backbone motion
	mm->set_chi(true); //Allow side-chain motion
	mm->set_jump(false); //Don't bother with "jumps" (relative motion of different molecules relative one another)
	core::optimization::AtomTreeMinimizer minimizer; //A torsion-space minimizer
	core::optimization::CartesianMinimizer cminimizer; //A cartesian-space minimizer

	core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	const string sequence = "APACDEFHIKLMNQRSTVWY"; //Use all amino acids at least once.  I put proline first so that I can put all the residues in an alpha-helical conformation.
	const string sequence2 = "A[DALA]P[DPRO]A[DALA]C[DCYS]D[DASP]E[DGLU]F[DPHE]H[DHIS]I[DILE]K[DLYS]L[DLEU]M[DMET]N[DASN]Q[DGLN]R[DARG]S[DSER]T[DTHR]V[DVAL]W[DTRP]Y[DTYR]"; //Use all amino acids at least once.  I put proline first so that I can put all the residues in an alpha-helical conformation.
	//const string sequence = "GAPACDEFHIKLA[DALA]MNQRSTVWYG"; //Use all amino acids at least once.  I put proline first so that I can put all the residues in an alpha-helical conformation.
	//const string sequence2 = "GA[DALA]P[DPRO]A[DALA]C[DCYS]D[DASP]E[DGLU]F[DPHE]H[DHIS]I[DILE]K[DLYS]L[DLEU]AM[DMET]N[DASN]Q[DGLN]R[DARG]S[DSER]T[DTHR]V[DVAL]W[DTRP]Y[DTYR]G"; //Use all amino acids at least once.  I put proline first so that I can put all the residues in an alpha-helical conformation.
	//const string sequence = "GAAAAAEAAARAAG";
	//const string sequence2 = "GA[DALA]A[DALA]A[DALA]A[DALA]A[DALA]E[DGLU]A[DALA]A[DALA]A[DALA]R[DARG]A[DALA]A[DALA]G";
	//const string sequence = "AAAAA";
	//const string sequence2 = "A[DALA]A[DALA]A[DALA]A[DALA]A[DALA]";

	core::pose::Pose Lpose;
	Lpose.clear();
	core::pose::make_pose_from_sequence(Lpose, sequence, *standard_residues, true);

	//Set backbone dihedrals:
	printf("Setting backbone dihedrals for L-peptide.\n");
	for(core::Size ir=1; ir<=Lpose.n_residue(); ir++) {
		if(ir>1) Lpose.set_phi(ir, -60.0);
		if(ir<Lpose.n_residue()) {
			Lpose.set_psi(ir, -45.0);
			Lpose.set_omega(ir, 175.0);
		}
	}
	Lpose.update_residue_neighbors();

	//Repack:
	repack_sc.apply(Lpose);

	//Fast relax:
	frlx.apply(Lpose);
	
	//Replace His8 with HIS_D:
	protocols::simple_moves::MutateResidue mutres(8, "HIS_D");
   mutres.apply(Lpose);

	//Score:
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
	(*sfxn)(Lpose);
	printf("The L-peptide's energy is %.4f\n", Lpose.energies().total_energy());
	sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
	(*sfxn)(Lpose);
	printf("The L-peptide's energy with cart_bonded is %.4f\n", Lpose.energies().total_energy());
	sfxn->set_weight(cart_bonded, 0.0); //Turn off cart_bonded.

	//Output:
	Lpose.dump_scored_pdb("Lpeptide.pdb", *sfxn);

	//Store a copy
	core::pose::Pose Lpose_unrelaxed = Lpose;

	//Minimize:
	printf("Attempting minimization.\n");
	//minimizer.run( Lpose, *mm, *sfxn, minoptions );
	sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
	cminimizer.run( Lpose, *mm, *sfxn, minoptions);
	sfxn->set_weight(cart_bonded, 0.0); //Turn on cart_bonded.
	//frlx.apply(Lpose);
	(*sfxn)(Lpose);
	printf("The L-peptide's energy after minimization is %.4f\n", Lpose.energies().total_energy());

	//Output:
	Lpose.dump_scored_pdb("Lpeptide_min.pdb", *sfxn);
	//sfxn->set_weight(cart_bonded, 0.0); //Turn on cart_bonded.

	//Mirror image pose:
	core::pose::Pose Dpose;
	Dpose.clear();
	core::pose::make_pose_from_sequence(Dpose, sequence2, *standard_residues, false);
	core::pose::add_lower_terminus_type_to_pose_residue(Dpose, 1);
	core::pose::add_upper_terminus_type_to_pose_residue(Dpose, Dpose.n_residue());
	//Replace DHis8 with DHIS_D:
	protocols::simple_moves::MutateResidue mutres2(8, "DHIS_D");
   mutres2.apply(Dpose);

/*	//Set backbone dihedrals:
	printf("Setting backbone dihedrals for D-peptide.\n");
	for(core::Size ir=1; ir<=Dpose.n_residue(); ir++) {
		if(ir>1) Dpose.set_phi(ir, 64.0);
		if(ir<Dpose.n_residue()) {
			Dpose.set_psi(ir, 41.0);
			Dpose.set_omega(ir, 180.0);
		}
	}

	//Set side-chains:
	printf("Settting side-chain dihedrals for D-peptide.\n");
	for(core::Size ir=1; ir<=Lpose.n_residue(); ir++) {
		printf("Res %c%i:", Lpose.residue(ir).name1(), (int)ir);
		fflush(stdout);
		if (Lpose.residue(ir).nchi() < 1) {printf("\n"); continue;}
		for(core::Size ichi=1; ichi<=Dpose.residue(ir).nchi(); ichi++) {
			Dpose.set_chi(ichi, ir, -1.0*Lpose.chi(ichi, ir));
			printf("\tchi%i=%.2f", (int)ichi, Lpose.chi(ichi, ir));
		}
		printf("\n");
	}
*/

	//Set backbone dihedrals:
	/*printf("Setting backbone dihedrals for D-peptide.\n");
	for(core::Size ir=1; ir<=Dpose.n_residue(); ir++) {
		if(ir>1) Dpose.set_phi(ir, 64.0);
		if(ir<Dpose.n_residue()) {
			Dpose.set_psi(ir, 41.0);
			Dpose.set_omega(ir, 180.0);
		}
	}

	//Copy and mirror coordinates of C-terminal glycine:
	for(core::Size ia=1; ia<=Dpose.residue(Dpose.n_residue()).natoms(); ia++) {
		numeric::xyzVector <core::Real> xyzcoord = Lpose_unrelaxed.xyz(id::AtomID(ia, Lpose_unrelaxed.n_residue()));
		xyzcoord[2] *= -1.0;
		Dpose.conformation().set_xyz(id::AtomID(ia, Dpose.n_residue()), xyzcoord);
	}

	Dpose.update_residue_neighbors();

	printf("Attempting repack\n");
	repack_sc.apply(Dpose);*/

	printf("Mirroring pose\n");
	for(core::Size ir=1; ir<=Lpose_unrelaxed.n_residue(); ir++){
		for(core::Size ia=1; ia<=Lpose_unrelaxed.residue(ir).natoms(); ia++) {
			numeric::xyzVector<core::Real> xyzvect = Lpose_unrelaxed.residue(ir).xyz(ia);
			xyzvect[0]*=-1.0;
			Dpose.set_xyz(AtomID(ia,ir), xyzvect);
		}
	}

	Dpose.update_residue_neighbors();

	//Score:
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
	(*sfxn)(Dpose);
	printf("The D-peptide's energy is %.4f\n", Dpose.energies().total_energy());
	sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
	(*sfxn)(Dpose);
	printf("The D-peptide's energy with cart_bonded is %.4f\n", Dpose.energies().total_energy());
	sfxn->set_weight(cart_bonded, 0.0); //Turn off cart_bonded.


	//Output:
	Dpose.dump_scored_pdb("Dpeptide.pdb", *sfxn);

	//Minimize:
	printf("Attempting minimization.\n");
	//minimizer.run( Dpose, *mm, *sfxn, minoptions );
	sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
	cminimizer.run( Dpose, *mm, *sfxn, minoptions);
	sfxn->set_weight(cart_bonded, 0.0); //Turn off cart_bonded.
	//frlx.apply(Dpose);
	(*sfxn)(Dpose);
	printf("The D-peptide's energy after minimization is %.4f\n", Dpose.energies().total_energy());
	
	//Output:
	Dpose.dump_scored_pdb("Dpeptide_min.pdb", *sfxn);
	//sfxn->set_weight(cart_bonded, 0.0); //Turn on cart_bonded.

	printf("\nSequence\nL-peptide\tD-peptide\n");
	for(core::Size ir=1, nres=Dpose.n_residue(); ir<=nres; ir++) {
		printf("%s\t%s\n", Lpose.residue(ir).name().c_str(), Dpose.residue(ir).name().c_str());
	}
	printf("\n");

	printf("***JOB COMPLETED.***\n");
	return 0;
}

