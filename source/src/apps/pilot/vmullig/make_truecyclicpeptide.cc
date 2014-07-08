/*
	make_truecyclicpeptide.cc
	A pilot app to create a true cyclic peptide (with a peptide bond between the N- and C-termini) 
	in a random conformation.  This is used for generating hundreds of thousands of conformations
	in an attempt to enumrate possible conformations exhaustively.
	Created by Vikram K. Mulligan, Baker Laboratory.

	File history:
		--Created 31 May 2013.
		--Added support for D-amino acids, 20 June 2013.
		--Added support for random mixes of ALA and VAL, 21 Aug 2013.
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
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
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
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

OPT_KEY( Integer, peptidesize )
OPT_KEY( Real, energycutoff )
OPT_KEY( Integer, relaxrnds )
OPT_KEY( IntegerVector, dpositions )
OPT_KEY( Real, valfraction )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1 <core::Size> emptyintlist;

	NEW_OPT( peptidesize, "The number of residues in the peptide, including terminal cysteines.  Default 9.", 9);
	NEW_OPT( energycutoff, "The energy above which structures are discarded.  Only used if a value is specified by the user.", 1000.0);
	NEW_OPT( relaxrnds, "The number of rounds of relaxation with the fastrelax mover.  Default 3.", 3);
	NEW_OPT( dpositions, "The positions at which a D-alanine should be placed.  Default empty list.", emptyintlist);
	NEW_OPT( valfraction, "The fraction of positions occupied by a valine.  Default 0.  (Range 0 to 1).", 0.0);
}

/********************************************************************************
	Function to move a residue from one pose to the position of a residue
	from another pose.  This copies the x, y, and z coordinates of each
	atom in the residue.  The reference residue must have an atom with an
	atom name matching each atom in the residue to alter.
********************************************************************************/
void copyresidueposition (
	core::pose::Pose &alteredpose,
	int alteredresidue,
	const core::pose::Pose &referencepose,
	int referenceresidue
) {
	using namespace core;
	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	using namespace std;
	using namespace chemical;
	using namespace conformation;
	using namespace core::pose;
	using namespace scoring::constraints;
	using namespace core::id;

	//Make a copy of the residue to modify:
	core::conformation::ResidueOP new_rsd = new core::conformation::Residue(alteredpose.residue(alteredresidue));

	for(core::Size ia=1; ia<=new_rsd->natoms(); ia++) {
		if (referencepose.residue(referenceresidue).has( new_rsd->atom_name(ia) )) {
			new_rsd->set_xyz(ia, referencepose.residue(referenceresidue).xyz( new_rsd->atom_name(ia) ));
		}
	}
	alteredpose.replace_residue(alteredresidue, *new_rsd, false);
	alteredpose.update_residue_neighbors();

	return;
}

/********************************************************************************
  Function to check whether an integer is in a list:
********************************************************************************/
bool inlist (
	core::Size number,
	utility::vector1 <core::Size> list
){
	if(list.size()==0) return false;
	for(core::Size i=1; i<=list.size(); i++) if(list[i]==number) return true;
	return false;
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

	printf("Starting make_truecyclicpeptide.\n");
	printf("Pilot app created 31 May 2013 by Vikram K. Mulligan, Baker Laboratory.\n");
	register_options(); //Parse the input flags and set global variables appropriately.
	devel::init(argc, argv); //Initialize Rosetta

	//Check that the user hasn't set an impossible value for peptidesize:
	if(option[peptidesize]()<5) {
		printf("Error!  The number of amino acid residues must be greater than or equal to 5.  Crashing gracelessly.\n");
		exit(1);
	}

	//D_positions is a list of all the positions that are d-amino acids.
	utility::vector1 <core::Size> d_positions;
	if(option[dpositions].user()) {
		if(option[dpositions]().size()>0) {
			for(core::Size i=1; i<=option[dpositions]().size(); i++) {
				if(d_positions.size()>0) {
					for(core::Size j=1; j<=d_positions.size(); j++) {
						if(option[dpositions]()[i] == d_positions[j]) {
							printf("Error!  The same residue occurs multiple times in the D-amino acid list.  Crashing gracelessly.\n");
							exit(1);
						}
					}
				}
				d_positions.push_back(option[dpositions]()[i]);
			}
		}
	}

	//Get the output prefix
	std::string outprefix = "out_";
	if(option[out::file::o].user()) outprefix=option[out::file::o]()+"_";

	//Check that the user has specified the number of structures:
	core::Size nstruct=1;
	if(option[out::nstruct].user()) nstruct=option[out::nstruct]();

	//Creating Rosetta's scorefunction.  This can be specified with flags by the user, and will otherwise default to the score12 scorefunction:
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	core::scoring::ScoreFunctionOP sfxn_constrained = core::scoring::get_score_function(); //For final fastrelax; includes constraint linking N- and C-termini.
	sfxn_constrained->set_weight(atom_pair_constraint, 1.0);
	sfxn_constrained->set_weight(dihedral_constraint, 1.0);
	sfxn_constrained->set_weight(angle_constraint, 1.0);

	//Some movers that I'll use.
	protocols::relax::FastRelaxOP frlx = new protocols::relax::FastRelax(sfxn_constrained, option[relaxrnds]()); //A fast relax mover that will do a user-specified number of rounds of fast relaxation.
	core::optimization::MinimizerOptions minoptions("dfpmin_armijo_nonmonotone", 0.00000001, false, false, false); //Options for the Cartesian minimizer
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap; //A move map tells an operator which parts of a structure (or what degrees of freedom) can move and which ones must remain fixed
	mm->set_bb(true); //Allow backbone motion
	mm->set_chi(true); //Allow side-chain motion
	mm->set_jump(false); //Don't bother with "jumps" (relative motion of different molecules relative one another)
	core::optimization::AtomTreeMinimizer minimizer; //A torsion-space minimizer

	//Create a Pose object in which to store our cyclic peptide structure:
	printf("Creating pose.\n");
	core::pose::Pose mypose;

	//Creating the peptide from an input sequence using the standard residues.  Note that this could easily be
	//modified to allow a user-input sequence.
	//TODO -- modify this a bit to allow proline or glycine at user-specified positions.
	core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	std::string myseq;
	if(inlist(1, d_positions)) myseq="A[DALA]";
	else myseq="A";
	for(core::Size ir=2; ir<=option[peptidesize](); ir++) {
		if(inlist(ir, d_positions)) myseq+="A[DALA]";
		else myseq+="A";
	}
	//core::pose::make_pose_from_sequence(mypose, myseq, *standard_residues, false); //The termini are OPEN.

	//Build the pose:
	{	chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( myseq, *standard_residues, false );
		mypose.clear();
		for ( core::Size i = 1; i <= requested_types.size(); ++i )
		{	chemical::ResidueType const & rsd_type = *requested_types[ i ];
			core::conformation::ResidueOP new_rsd( NULL );
			new_rsd = conformation::ResidueFactory::create_residue ( rsd_type );
			if (i>1) mypose.append_residue_by_bond( *new_rsd, true );
			else mypose.append_residue_by_jump(*new_rsd, 1, "", "", true);
		}
	}

	//Make the linkage between the N- and C-termini:
	mypose.conformation().declare_chemical_bond(1, "N", mypose.n_residue(), "C"); //Declare a chemical bond between the N and C termini.
	{ 
		//Peptide bond length constraint:
		mypose.add_constraint (
			new AtomPairConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ) ,
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				new HarmonicFunc( 1.3288, 0.01)
			)
		);

		//Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		mypose.add_constraint (
			new DihedralConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("O") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("H") , 1) ,
				new CircularHarmonicFunc( PI, 0.02)
			)
		);
		mypose.add_constraint (
			new DihedralConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("CA") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("CA") , 1) ,
				new CircularHarmonicFunc( PI, 0.02)
			)
		);

		//Peptide bond angle constraints:
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("CA") , 1) ,
				new CircularHarmonicFunc( CNCa_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("H") , 1) ,
				new CircularHarmonicFunc( CNH_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("CA") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				new CircularHarmonicFunc( CaCN_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("O") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				new CircularHarmonicFunc( OCN_ANGLE/180.0*PI, 0.02)
			)
		);

	}

	//Set all the omega angles:
	for(core::Size ir=1; ir<mypose.n_residue(); ir++) mypose.set_omega(ir, 180.0);
	mypose.update_residue_neighbors();

	//For testing only -- delete me!
	//(*sfxn)(mypose);
	//mypose.dump_scored_pdb("temp.pdb", *sfxn);
	//printf("Res 4 AA = %i\n", (int)mypose.aa(4));
	//printf("aa_dal = %i\n", (int)core::chemical::aa_dal); fflush(stdout);

	//Create a second pose that's just four linked alanines:
	//TODO -- update this if gly or pro are at the 1, 2, n-1, or n positions.
	core::pose::Pose alapose;
	std:string alaseq;
	if(inlist(mypose.n_residue()-1, d_positions)) alaseq+="A[DALA]"; else alaseq+="A";
	if(inlist(mypose.n_residue(), d_positions)) alaseq+="A[DALA]"; else alaseq+="A";
	if(inlist(1, d_positions)) alaseq+="A[DALA]"; else alaseq+="A";
	if(inlist(2, d_positions)) alaseq+="A[DALA]"; else alaseq+="A";
	//REPLACE below
	//core::pose::make_pose_from_sequence(alapose, alaseq, *standard_residues, true); //The termini are CLOSED.
	{	chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( alaseq, *standard_residues, false );
		alapose.clear();
		for ( core::Size i = 1; i <= requested_types.size(); ++i )
		{	chemical::ResidueType const & rsd_type = *requested_types[ i ];
			core::conformation::ResidueOP new_rsd( NULL );
			new_rsd = conformation::ResidueFactory::create_residue ( rsd_type );
			if (i>1) alapose.append_residue_by_bond( *new_rsd, true );
			else {	alapose.append_residue_by_jump(*new_rsd, 1, "", "", true);
				core::pose::add_lower_terminus_type_to_pose_residue(alapose, 1);
			}
		}
	}
	for(core::Size ir=1; ir<alapose.n_residue(); ir++) alapose.set_omega(ir, 180.0); //Set the omega angles.
	
	core::pose::add_upper_terminus_type_to_pose_residue(alapose, alapose.n_residue());

	//Kinematic closure mover:
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kinmover = new protocols::loops::loop_closure::kinematic_closure::KinematicMover;
	kinmover->set_temperature( 1.0 );
	kinmover->set_vary_bondangles( false );
	kinmover->set_sample_nonpivot_torsions( true ); //true
	kinmover->set_rama_check( false );
	kinmover->set_idealize_loop_first(true); //true
	kinmover->set_sfxn(sfxn);
	kinmover->set_pivots(2, (core::Size)(((core::Real)mypose.n_residue()-2.0)/2.0+2.0), mypose.n_residue()-1);

	for(core::Size rep=1; rep<=nstruct; rep++) {
		//make a copy of alapose:
		core::pose::Pose temp_alapose=alapose;

		//If a fraction of residues should be valines, mutate temp_alapose accordingly:
		if(option[valfraction].user()) {
			for(core::Size ir=1; ir<=temp_alapose.n_residue(); ir++) {
				if(RG.uniform() < option[valfraction]()) {
					string aaname = "VAL";
					if( core::chemical::is_D_aa( temp_alapose.residue(ir).aa() ) ) aaname="DVAL";
					protocols::simple_moves::MutateResidue mutres(ir, aaname);
					mutres.apply(temp_alapose);
				}
			}
		}

		//Randomize backbone:
		//TODO -- randomly sample cis-pro if there are prolines.
		for(core::Size ir=1; ir<=4; ir++) {
			if(ir>1) temp_alapose.set_phi(ir, ((core::Real)RG.uniform()-0.5)*360.0);
			if(ir<4) temp_alapose.set_psi(ir, ((core::Real)RG.uniform()-0.5)*360.0);
		}
		temp_alapose.update_residue_neighbors();

		frlx->apply(temp_alapose); //repack and minimize
		//minimizer.run( temp_alapose, *mm, *sfxn, minoptions );

		//Make a copy of mypose
		core::pose::Pose temppose = mypose;

		//If a fraction of residues should be valines, mutate accordingly:
		if(option[valfraction].user()) {
			for(core::Size ir=1; ir<=temppose.n_residue(); ir++) {
				if(ir<3 || ir>temppose.n_residue()-2) { //The first two and last two residues must match temp_alapose
					core::Size ir2=0;
					if(ir<3) ir2=ir+2;
					else ir2=ir-temppose.n_residue()+2;
					if(temp_alapose.residue(ir2).name3()=="DVA") {
						protocols::simple_moves::MutateResidue mutres(ir, "DVAL");
						mutres.apply(temppose);
					} else if (temp_alapose.residue(ir2).name3()=="VAL") {
						protocols::simple_moves::MutateResidue mutres(ir, "VAL");
						mutres.apply(temppose);
					}
				} else { //The rest of the residues must still be mutated
					if(RG.uniform() < option[valfraction]()) {
						string aaname = "VAL";
						if( core::chemical::is_D_aa( temppose.residue(ir).aa() ) ) aaname="DVAL";
						protocols::simple_moves::MutateResidue mutres(ir, aaname);
						mutres.apply(temppose);
					}
				}
			}
		}

		//Put the first two and last two residues of temppose in the same positions as the corresponding residues in temp_alapose.
		copyresidueposition(temppose, 1, temp_alapose, 3); //First residue
		copyresidueposition(temppose, 2, temp_alapose, 4); //Second residue
		copyresidueposition(temppose, temppose.n_residue(), temp_alapose, 2); //Last residue
		copyresidueposition(temppose, temppose.n_residue()-1, temp_alapose, 1); //Second-last residue

		//Update residue neighbours:
		temppose.update_residue_neighbors();

		//Kinematic closure to build the rest of the peptide:
		//printf("Applying kinematic closure.\n"); fflush(stdout); //DELETE ME
		kinmover->apply(temppose);
		temppose.update_residue_neighbors();

		//Check that the closure was successful.  If the gap is too big for the loop, it won't always be a successful closure.
		//printf("Checking successful closure.\n"); fflush(stdout); //DELETE ME
		bool closurefailed=false;
		for (core::Size ir=1; ir<temppose.n_residue(); ir++) {
			numeric::xyzVector<core::Real> bonddeltaC = temppose.residue(ir).xyz("C");
			numeric::xyzVector<core::Real> bonddeltaN = temppose.residue(ir+1).xyz("N");
			core::Real bonddeltalength = sqrt(pow(bonddeltaC[0]-bonddeltaN[0],2.0)+pow(bonddeltaC[1]-bonddeltaN[1],2.0)+pow(bonddeltaC[2]-bonddeltaN[2],2.0));
			if(bonddeltalength>1.4 || bonddeltalength<1.2) {
				closurefailed=true; 
				printf("Peptide bond length between residues %i and %i is %.2f.  Loop closure failed.  Randomizing dihedrals and trying again.\n", ir, ir+1, bonddeltalength);  fflush(stdout);
				rep--;
				break;
			}
		}
		if(closurefailed) continue;

		//Fast relax the structure:
		//printf("Relaxing structure.\n"); fflush(stdout); //DELETE ME
		frlx->apply(temppose);

		//Check that the peptide bond between residues 1 and n has sensible geometry:
		//printf("Checking geometry.\n"); fflush(stdout); //DELETE ME
		{	//Peptide bond length check:
			core::Real bonddeltalength = sqrt(	pow( temppose.residue(1).xyz("N")[0]-temppose.residue(temppose.n_residue()).xyz("C")[0] , 2) +
								pow( temppose.residue(1).xyz("N")[1]-temppose.residue(temppose.n_residue()).xyz("C")[1] , 2) +
								pow( temppose.residue(1).xyz("N")[2]-temppose.residue(temppose.n_residue()).xyz("C")[2] , 2)
							 );
			if(bonddeltalength>1.4 || bonddeltalength<1.2) {
				printf("Peptide bond length between first and last residues is %.2f.  Loop closure failed.  Randomizing dihedrals and trying again.\n", bonddeltalength);
				fflush(stdout);
				rep--;
				continue;
			}
		}
		{	//Peptide dihedral angle check:
			core::Real bonddihedralangle = numeric::dihedral_degrees(
				temppose.residue(temppose.n_residue()).xyz("CA"),
				temppose.residue(temppose.n_residue()).xyz("C"),
				temppose.residue(1).xyz("N"),
				temppose.residue(1).xyz("CA")
			);
			for(core::Size i=1; i<=2; i++) {
				if(bonddihedralangle>190.0 || (bonddihedralangle < 170.0 && bonddihedralangle > -170) || bonddihedralangle < -190) {
					printf("Peptide bond dihedral angle between first and last residues is %.2f degrees.  Loop closure failed.  Randomizing dihedrals and trying again.\n", bonddihedralangle);
					fflush(stdout);
					rep--;
					closurefailed=true;
					break;
				}
				bonddihedralangle = numeric::dihedral_degrees(
					temppose.residue(temppose.n_residue()).xyz("O"),
					temppose.residue(temppose.n_residue()).xyz("C"),
					temppose.residue(1).xyz("N"),
					temppose.residue(1).xyz("H")
				);
			}
		}
		if(closurefailed) continue;
		{	//Peptide bond angle checks:
			core::Real peptidebondangle = numeric::angle_degrees(
				temppose.residue(temppose.n_residue()).xyz("CA"),
				temppose.residue(temppose.n_residue()).xyz("C"),
				temppose.residue(1).xyz("N")
			);
			if(peptidebondangle > CaCN_ANGLE+10.0 || peptidebondangle < CaCN_ANGLE-10.0) {
					printf("Peptide bond between first and last residues has a bad CA-C-N bond angle (%.2f degrees).  Loop closure failed.  Randomizing dihedrals and trying again.\n", peptidebondangle);
					fflush(stdout);
					rep--;
					continue;
			}

			peptidebondangle = numeric::angle_degrees(
				temppose.residue(temppose.n_residue()).xyz("O"),
				temppose.residue(temppose.n_residue()).xyz("C"),
				temppose.residue(1).xyz("N")
			);
			if(peptidebondangle > OCN_ANGLE+10.0 || peptidebondangle < OCN_ANGLE-10.0) {
					printf("Peptide bond between first and last residues has a bad O-C-N bond angle (%.2f degrees).  Loop closure failed.  Randomizing dihedrals and trying again.\n", peptidebondangle);
					fflush(stdout);
					rep--;
					continue;
			}

			peptidebondangle = numeric::angle_degrees(
				temppose.residue(temppose.n_residue()).xyz("C"),
				temppose.residue(1).xyz("N"),
				temppose.residue(1).xyz("CA")
			);
			if(peptidebondangle > CNCa_ANGLE+10.0 || peptidebondangle < CNCa_ANGLE-10.0) {
					printf("Peptide bond between first and last residues has a bad C-N-CA bond angle (%.2f degrees).  Loop closure failed.  Randomizing dihedrals and trying again.\n", peptidebondangle);
					fflush(stdout);
					rep--;
					continue;
			}

			peptidebondangle = numeric::angle_degrees(
				temppose.residue(temppose.n_residue()).xyz("C"),
				temppose.residue(1).xyz("N"),
				temppose.residue(1).xyz("H")
			);
			if(peptidebondangle > CNH_ANGLE+10.0 || peptidebondangle < CNH_ANGLE-10.0) {
					printf("Peptide bond between first and last residues has a bad C-N-H bond angle (%.2f degrees).  Loop closure failed.  Randomizing dihedrals and trying again.\n", peptidebondangle);
					fflush(stdout);
					rep--;
					continue;
			}
			
		}

		//Score the structure and check that the energy isn't ridiculous:
		//printf("Scoring structure.\n"); fflush(stdout); //DELETE ME
		(*sfxn_constrained)(temppose);
		if (option[energycutoff].user() && (temppose.energies().total_energy() > option[energycutoff]())) {
			printf("The final energy is %.1f, which is above the cutoff.  Discarding and trying again.\n", temppose.energies().total_energy()); fflush(stdout);
			rep--;
			continue;
		}

		//printf("Writing output.\n"); fflush(stdout); //DELETE ME
		char outfile[64];
		sprintf(outfile, "%s%04i.pdb", outprefix.c_str(), rep);
		temppose.dump_scored_pdb(outfile, *sfxn_constrained); //Write out the final scored pdb file.
		printf("Wrote %s\n", outfile); fflush(stdout);
	}

	printf("***JOB COMPLETED.***\n");
	return 0;
}

