/****************************************************************************************************
PCA_visualizer.cc

A pilot app to read a PCA file and a PDB, and to generate a series of PDB files animating
the perturbation of the initial PDB file by the PCA vectors.
History:
--File created 22 May 2013 by Vikram K. Mulligan, Baker Laboratory.
--Modified 10 Sept 2013 by VKM to allow beta-amino acids.
--Modified 15 Oct 2013 by VKM to work with jumps.
****************************************************************************************************/

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/cluster/cluster.hh>
#include <protocols/loops/Loops.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <iostream>
#include <string>
#include <deque>
//#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/format.hh>
#include <core/pose/annotated_sequence.hh>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/EulerAngles.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
//#include <protocols/relax/FastRelax.hh>
//#include <numeric/random/random.hh>
//#include <numeric/random/uniform.hh>

#define PI 3.1415926535897932384626433832795028841971693993751058

using ObjexxFCL::format::F;

OPT_KEY (FileVector, PCAfiles) //List of PCA files to open.
OPT_KEY (Integer, outcount) //Number of output files to write per input structure -- default 30.
OPT_KEY (RealVector, amplitudes) //List of oscillation amplitudes.
OPT_KEY (RealVector, frequencies) //List of oscillation frequencies.
OPT_KEY (RealVector, offsets) //List of oscillation offsets.
OPT_KEY (Boolean, use_variances) //Use the variances as multipliers?  Default true.
OPT_KEY (Boolean, do_repack) //Repack side-chains?
OPT_KEY (Boolean, do_minimize) //Minimize side-chains?
OPT_KEY (Boolean, cyclic) //Is this a cyclic peptide?

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1 < std::string > emptystringvector;
	utility::vector1 < core::Real > emptyrealvector;

	NEW_OPT(PCAfiles, "List of PCA files to load.  Must be specified in same order as files listed using the -in:file:s flag.", emptystringvector);
	NEW_OPT(outcount, "Number of output files to write per input structure.  (Or length of animation, in frames.)  Default 30.", 30);
	NEW_OPT(amplitudes, "List of amplitudes of oscillation for the PCA vectors.  If fewer amplitudes are specified than there are PCA vectors, amplitudes of 0 are assumed for all later PCA vectors.", emptyrealvector);
	NEW_OPT(frequencies, "List of frequencies of oscillation for the PCA vectors (in oscillation cycles per animation timecourse).", emptyrealvector);
	NEW_OPT(offsets, "List of offsets of oscillation for the PCA vectors (as a fraction of an oscillation cycle).", emptyrealvector);
	NEW_OPT(use_variances, "If true, amplitudes are premultiplied by variances for each principal component vector.  If false, principal component vectors are used as-is.  Default true.", true);
	NEW_OPT(do_repack, "If true, residues are repacked after each perturbation.  True by default.", true);
	NEW_OPT(do_minimize, "If true, side-chains are minimized after each perturbation.  True by default.", true);
	NEW_OPT(cyclic, "If true, the peptide is treated as cyclic, and extra backbone dihedrals are used.  False by default.", false);
}

//Functions to set the backbone dihedral angles of beta-amino acid peptides:
void betapeptide_setphi (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 1 ), angle );
	return;
}

void betapeptide_settheta (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 2 ), angle );
	return;
}

void betapeptide_setpsi (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 3 ), angle );
	return;
}

void betapeptide_setomega (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 4 ), angle );
	return;
}

/********************************************************************************
Function to read in a PCA file and output a variances vector and a
vector of principal component vectors.
********************************************************************************/
void read_PCAfile (
	const char* filename,
	utility::vector1 <core::Real> &variances,
	utility::vector1 < utility::vector1 < core::Real > > &PCA_vectors,
	bool normalize_variances
) {
	//Clear the storage boxes:
	variances.clear();
	PCA_vectors.clear();

	//Open the PCA file for reading
	FILE* pcafile = fopen(filename, "r");

	//Get the size of each vector, and the number of vectors
	core::Size vectorcount, vectorcomponentcount;
	fscanf(pcafile, "%lu %lu", &vectorcomponentcount, &vectorcount);

	//If there are no principal component vectors for this state, return without doing anything.
	//The sizes of the variancs and PCA_vectors vectors will be 0, in this case.
	if ( vectorcount==0 ) return;

	//Resize the variances and PCA_vectors vectors:
	variances.resize(vectorcount);
	PCA_vectors.resize(vectorcount);

	//Loop through and read the variances vector:
	core::Real realbuffer;
	for ( core::Size i=1; i<=vectorcount; i++ ) {
		fscanf(pcafile, "%lf", &realbuffer);
		variances[i]=realbuffer;
	}

	//Loop through and populate the PCA matrix:
	for ( core::Size i=1; i<=vectorcount; i++ ) {
		PCA_vectors[i].resize(vectorcomponentcount);
		for ( core::Size j=1; j<=vectorcomponentcount; j++ ) {
			fscanf(pcafile, "%lf", &realbuffer);
			PCA_vectors[i][j]=realbuffer;
		}
	}

	//Close the file:
	fclose(pcafile);

	//Normalize the variances vector to the size of the first entry:
	if ( normalize_variances ) { for ( core::Size i=variances.size(); i>=1; i-- ) variances[i]=variances[i]/variances[1];
	}

	return;
}

//Function to fix the termini in the case of cyclic peptides:
void fix_cyclic_termini (core::pose::Pose &mypose) {
	using namespace core;
	using namespace core::id;

	core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(mypose, mypose.size());
	core::Size hindex = mypose.residue(1).atom_index("H");

	//Rebuild the N-terminal proton.  This has to be done in a slightly irritating way because Rosetta doesn't really like the fact
	//that the last residue is connected to the first:
	{
		core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::pose::Pose dialanine;
		std::string diala_seq = "";
		if ( mypose.residue(mypose.size()).has("CM") ) diala_seq+="A[B3A]"; else diala_seq+="A";
		if ( mypose.residue(1).has("CM") ) diala_seq+="A[B3A]"; else diala_seq+="A";

		//core::pose::make_pose_from_sequence(dialanine, "AA", *standard_residues, true); //The termini are OPEN.
		{
			core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
			core::chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( diala_seq, *standard_residues, false );
			for ( core::Size ir=1, numres=requested_types.size(); ir<=numres; ir++ ) {
				core::chemical::ResidueType const & rsd_type = *requested_types[ ir ];
				core::conformation::ResidueOP new_rsd = core::conformation::ResidueFactory::create_residue( rsd_type );
				if ( ir==1 ) dialanine.append_residue_by_jump( *new_rsd, 1, "", "", true );
				else dialanine.append_residue_by_bond(*new_rsd, true);
			}
		}


		core::Real omegaval = numeric::dihedral_degrees(
			mypose.residue(mypose.size()).xyz( (mypose.residue(mypose.size()).has("CM")?"CM":"CA") ),
			mypose.residue(mypose.size()).xyz("C"),
			mypose.residue(1).xyz("N"),
			mypose.residue(1).xyz("CA")
		);

		dialanine.conformation().set_torsion(id::TorsionID(1, id::BB, (dialanine.residue(1).has("CM")?4:3) ), omegaval);
		//dialanine.set_omega(1, omegaval);

		core::id::AtomID_Map< core::id::AtomID > amap;
		core::pose::initialize_atomid_map(amap,dialanine,core::id::AtomID::BOGUS_ATOM_ID());
		amap[AtomID(dialanine.residue(1).atom_index((dialanine.residue(1).has("CM")?"CM":"CA")),1)] = AtomID(mypose.residue(mypose.size()).atom_index((dialanine.residue(1).has("CM")?"CM":"CA")),mypose.size());
		amap[AtomID(dialanine.residue(1).atom_index("C"),1)] = AtomID(mypose.residue(mypose.size()).atom_index("C"),mypose.size());
		amap[AtomID(dialanine.residue(1).atom_index("O"),1)] = AtomID(mypose.residue(mypose.size()).atom_index("O"),mypose.size());
		amap[AtomID(dialanine.residue(2).atom_index("N"),2)] = AtomID(mypose.residue(1).atom_index("N"),1);
		amap[AtomID(dialanine.residue(2).atom_index("CA"),2)] = AtomID(mypose.residue(1).atom_index("CA"),1);
		core::scoring::superimpose_pose( dialanine, mypose, amap );

		mypose.conformation().set_xyz(AtomID(hindex, 1), dialanine.residue(2).xyz("H"));
	}

	mypose.conformation().declare_chemical_bond(1, "N", mypose.size(), "C"); //Declare a chemical bond between the N and C termini.
}

//Function to count the number of backbone dihedral angles in an input structure, and to confirm that these match the PCA file.
bool checkresiduecount (
	const core::pose::Pose &mypose,
	const utility::vector1 < core::Real > &pcavect
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size count=0;

	for ( core::Size ir=1, nres=mypose.size(); ir<=nres; ir++ ) {
		if ( mypose.residue(ir).type().is_beta_aa() ) {
			count+=4; //Beta-amino acid
			if ( mypose.residue(ir).is_lower_terminus() || ir==1 ) count--;
			if ( mypose.residue(ir).is_upper_terminus() || ir==nres ) count-=2;
		} else if ( mypose.residue(ir).type().is_alpha_aa() ) {
			count+=3; //Alpha amino acid
			if ( mypose.residue(ir).is_lower_terminus() || ir==1 ) count--;
			if ( mypose.residue(ir).is_upper_terminus() || ir==nres ) count-=2;
		} else { } //Default case -- do nothing.
	}
	if ( option[cyclic]() ) count+=3; //Add three if the peptide is cyclic.
	else if ( mypose.num_jump() > 0 ) count+=6*mypose.num_jump(); //Add jump data.

	if ( count==pcavect.size() ) return true;

	return false;
}

using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;
using namespace basic::options;

int main( int argc, char * argv [] ) {

	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	using namespace protocols::cluster;
	using namespace std;

	printf("Starting PCA_visualizer.cc\nFile created 22 May 2013 by Vikram K. Mulligan\n\n"); fflush(stdout);

	register_options();
	devel::init(argc, argv);

	//Initial checks:
	if ( option[in::file::s].size()==0 ) {
		printf("Error!  The user must specify input PDB files with the -in::file::s flag.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}
	if ( option[PCAfiles].size()==0 ) {
		printf("Error!  The user must specify input PCA files with the -PCAfiles flag.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}
	if ( option[PCAfiles].size()!=option[in::file::s].size() ) {
		printf("Error!  The number of PDB files must match the number of PCA files.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}
	if ( option[amplitudes].size()==0 ) {
		printf("Error!  At least one amplitude must be specified using the -amplitudes flag.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}
	if ( option[frequencies].size()==0 ) {
		printf("Error!  At least one frequency must be specified using the -frequencies flag.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}
	if ( option[offsets].size()==0 ) {
		printf("Error!  At least one offset must be specified using the -offsets flag.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}
	if ( option[offsets].size()!=option[amplitudes].size() || option[frequencies].size()!=option[amplitudes].size() ) {
		printf("Error!  The lists of amplitudes, frequencies, and offsets must all have the same number of elements.\nCrashing.\n"); fflush(stdout);
		exit(1);
	}

	//Vectors for the variances, the principal component vectors, and the perturbation vector:
	utility::vector1 < core::Real > variances;
	utility::vector1 < utility::vector1 < core::Real > > PCA_vectors;
	utility::vector1 < core::Real > pertvector;

	//Lists of the file names to import:
	utility::vector1 <std::string> pdbfiles = option[in::file::s]();
	utility::vector1 <std::string> pcafiles = option[PCAfiles]();

	//Storage for the imported and perturbed poses:
	core::pose::Pose inputpose;
	core::pose::Pose perturbedpose;

	//List of the frequencies:
	utility::vector1 <core::Real> freqs = option[frequencies]();
	for ( core::Size i=1; i<=freqs.size(); i++ ) freqs[i]=freqs[i]*2.0*PI; //Pre-multiply by 2*PI.

	//List of the offsets:
	utility::vector1 <core::Real> offs = option[offsets]();

	//Side-chains only movemap
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_bb(false);
	mm->set_chi(true);
	mm->set_jump(false);

	//Score function:
	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::get_score_function();

	//Repack and minimize movers:
	protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn); //Create the RepackSidechains mover and set the score function.
	protocols::simple_moves::MinMover minmove(mm, sfxn, "dfpmin_armijo_nonmonotone", 0.0000000001, true);

	for ( core::Size istruct=1; istruct<=pcafiles.size(); istruct++ ) {

		//Import the PCA vectors for this structure:
		read_PCAfile(pcafiles[istruct].c_str(), variances, PCA_vectors, true);

		if ( variances.size()==0 ) {
			printf("No principal component vectors found in file %s.  Moving on to next structure.\n", pcafiles[istruct].c_str() ); fflush(stdout);
			continue;
		} else {
			printf("Read %lu principal component vectors from file %s.\n", variances.size(), pcafiles[istruct].c_str()); fflush(stdout);
		}

		//Import the pose:
		core::import_pose::pose_from_file(inputpose, pdbfiles[istruct], core::import_pose::PDB_file);
		//Check that the number of residues matches the perturbation file:
		if ( !checkresiduecount(inputpose, PCA_vectors[1]) ) {
			printf("Error!  Input pdb %s is a different size than the backbones corresponding to PCA file %s.\nCrashing.\n", pcafiles[istruct].c_str(), pdbfiles[istruct].c_str()); fflush(stdout);
			exit(1);
		}

		if ( option[cyclic]() ) fix_cyclic_termini(inputpose); //Fix the termini if this is a backbone-cyclized peptide.

		//Lists of the amplitudes:
		utility::vector1 <core::Real> amps;
		amps.resize( std::min( option[amplitudes]().size(), variances.size() ) );
		//Premultiply the amplitudes if and only if use_variances is true; otherwise, use them straight.
		for ( core::Size i=1, ampsize=amps.size(); i<=ampsize; i++ ) amps[i] = (option[use_variances] ? (variances[i]*option[amplitudes][i]) : (option[amplitudes][i]) );

		for ( core::Size ioutput=0; ioutput < static_cast<core::Size>(option[outcount]()); ioutput++ ) {
			//Copy the input pose:
			perturbedpose = inputpose;

			//Resize the perturbation vector to match the PCA vectors, and initialize it with zeroes:
			pertvector.clear();
			pertvector.resize(PCA_vectors[1].size(), 0.0);

			//The variable t runs from 0 (inclusive) to 1 (exclusive) over the course of the animation:
			core::Real t = ((core::Real) ioutput)/((core::Real) option[outcount]());

			//Loop through all amplitudes, whichever is smaller.
			for ( core::Size ivect=1; ivect<=amps.size(); ivect++ ) {
				//Loop through all entries in the peturbation vector and generate it from the PCA vectors.
				for ( core::Size j=1; j<=pertvector.size(); j++ ) {
					pertvector[j] += amps[ivect]*sin(freqs[ivect]*(t+offs[ivect]))*PCA_vectors[ivect][j];
				}
			}

			//Once the perturbation vector has been generated, loop through all residues and perturb the backbone:
			core::Size pertindex = 1;
			for ( core::Size ir=1; ir<=perturbedpose.size(); ir++ ) {
				if ( perturbedpose.residue(ir).type().is_beta_aa() ) { //beta-amino acids
					if ( ir>1 && !perturbedpose.residue(ir).is_lower_terminus() ) betapeptide_setphi(perturbedpose, ir, perturbedpose.residue(ir).mainchain_torsion(1)+pertvector[pertindex++]);
					betapeptide_settheta(perturbedpose, ir, perturbedpose.residue(ir).mainchain_torsion(2)+pertvector[pertindex++]);
					if ( ir<perturbedpose.size() && !perturbedpose.residue(ir).is_upper_terminus() ) {
						betapeptide_setpsi(perturbedpose, ir, perturbedpose.residue(ir).mainchain_torsion(3)+pertvector[pertindex++]);
						betapeptide_setomega(perturbedpose, ir, perturbedpose.residue(ir).mainchain_torsion(4)+pertvector[pertindex++]);
					}
				} else if ( perturbedpose.residue(ir).type().is_alpha_aa() ) { //alpha-amino acids:
					if ( ir>1 && !perturbedpose.residue(ir).is_lower_terminus() ) perturbedpose.set_phi(ir, perturbedpose.phi(ir)+pertvector[pertindex++]);
					if ( ir<perturbedpose.size()  && !perturbedpose.residue(ir).is_lower_terminus() ) {
						perturbedpose.set_psi(ir, perturbedpose.psi(ir)+pertvector[pertindex++]);
						perturbedpose.set_omega(ir, perturbedpose.omega(ir)+pertvector[pertindex++]);
					}
				} else { } //Default case -- do nothing.
			}

			//If this is a cyclic peptide, perturb the backbone dihedral angles adjacent to the peptide bond linking the N- and C-termini
			if ( option[cyclic]() ) {
				core::Real ang1 = numeric::dihedral_degrees (
					perturbedpose.residue(perturbedpose.size()).xyz( (perturbedpose.residue(perturbedpose.size()).has("CM")?"CA":"N") ),
					perturbedpose.residue(perturbedpose.size()).xyz( (perturbedpose.residue(perturbedpose.size()).has("CM")?"CM":"CA") ),
					perturbedpose.residue(perturbedpose.size()).xyz("C"),
					perturbedpose.residue(perturbedpose.size()).xyz( "O")
				);
				core::Real ang2 = numeric::dihedral_degrees (
					perturbedpose.residue(1).xyz("H"),
					perturbedpose.residue(1).xyz("N"),
					perturbedpose.residue(1).xyz("CA"),
					perturbedpose.residue(1).xyz( (perturbedpose.residue(1).has("CM")?"CM":"C") )
				);
				perturbedpose.conformation().set_torsion_angle(//Set psi of last residue
					core::id::AtomID(perturbedpose.residue(perturbedpose.size()).atom_index( (perturbedpose.residue(perturbedpose.size()).has("CM")?"CA":"N") ), perturbedpose.size()),
					core::id::AtomID(perturbedpose.residue(perturbedpose.size()).atom_index( (perturbedpose.residue(perturbedpose.size()).has("CM")?"CM":"CA") ), perturbedpose.size()),
					core::id::AtomID(perturbedpose.residue(perturbedpose.size()).atom_index("C"), perturbedpose.size()),
					core::id::AtomID(perturbedpose.residue(perturbedpose.size()).atom_index("O"), perturbedpose.size()),
					(ang1+pertvector[pertindex++])/180.0*PI);
				pertindex++;
				perturbedpose.conformation().set_torsion_angle( //Set phi of first residue
					core::id::AtomID(perturbedpose.residue(1).atom_index("H"), 1),
					core::id::AtomID(perturbedpose.residue(1).atom_index("N"), 1),
					core::id::AtomID(perturbedpose.residue(1).atom_index("CA"), 1),
					core::id::AtomID(perturbedpose.residue(1).atom_index( (perturbedpose.residue(1).has("CM")?"CM":"C") ), 1),
					(ang2+pertvector[pertindex++])/180.0*PI);
			} else if ( perturbedpose.num_jump() > 0 ) { //Otherwise, if this is a multi-chain pose, perturb the jumps:
				for ( core::Size jumpcount=1, numjump=perturbedpose.num_jump(); jumpcount<=numjump; jumpcount++ ) {
					numeric::xyzVector < core::Real > transvect = perturbedpose.jump(jumpcount).get_translation();
					numeric::EulerAngles < core::Real > euler_angles(perturbedpose.jump(jumpcount).get_rotation());
					for ( core::Size j=0; j<3; j++ ) transvect[j]+=pertvector[pertindex++]; //Perturb translation (note zero-based vector!)
					//Perturb rotation:
					euler_angles.phi_degrees( euler_angles.phi_degrees() + pertvector[pertindex++] );
					euler_angles.theta_degrees( euler_angles.theta_degrees() + pertvector[pertindex++] );
					euler_angles.psi_degrees( euler_angles.psi_degrees() + pertvector[pertindex++] );
					//Update the jump:
					core::kinematics::Jump myjump = perturbedpose.jump(jumpcount);
					myjump.set_translation( transvect);
					myjump.set_rotation( euler_angles.to_rotation_matrix() );
					perturbedpose.set_jump(jumpcount, myjump);
				}
			}

			if ( option[do_repack]() ) repack_sc.apply(perturbedpose);
			if ( option[do_minimize]() ) minmove.apply(perturbedpose);

			char outfile [256];
			sprintf(outfile, "PCA_%04lu_%04lu.pdb", istruct, ioutput+1);
			perturbedpose.dump_pdb(outfile);
			printf("Wrote %s\n", outfile); fflush(stdout);
		}
	}

	printf("JOB COMPLETED.\n"); fflush(stdout);
	return 0;
}

