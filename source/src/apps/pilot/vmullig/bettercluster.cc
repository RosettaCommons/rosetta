/****************************************************************************************************
	bettercluster.cc

	A better clustering algorithm for generating conformational libraries for
	explicit multistate design.
	History:
	--File created 6 May 2013 by Vikram K. Mulligan, Baker Laboratory.
	--Modified 3 Jun 2013 to make Cartesian-based clustering much faster.  (No
	more rebuilding poses).
	--Modified 5 Jun 2013 to allow clustering of backbone-cyclized peptides.
	--Modified 28 Aug 2013 to allow N-offset cyclic permutations to be clustered.  (e.g. Offsets
		are incremented by 2 residues, or by 3 residues, or whatever.)
	--Modified 4 Sept 2013 to allow clustering of peptides with beta-amino acid residues:
			--check for beta residues based on first structure loaded (they must be consistent -- DONE).
			--align CM atoms, if present (DONE).
			--update information that is stored for cartesian or dihedral clustering (DONE).
			--update PCA file output --> list of backbone dihedrals must include theta (DONE, though
			other apps must also be updated to READ these properly).
	--Modified 9 Oct 2013 to allow clustering of homooligomers (with swapping of oligomer subunits
		during RMSD calculation).
				TODO:
				--store jumps for PCA analysis
				--calculate RMSD for all possible permutations of homooligomer subunits and keep lowest
	--Modified 21 May 2014 to allow silent file output.
	--Modified 22 May 2014:
		--Added support for constraints files.
		--Added support for a user-defined list of additional atoms to use in RMSD
		calculation.
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
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/relax/FastRelax.hh>
#include <numeric/model_quality/rms.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/annotated_sequence.hh>

#include <numeric/EulerAngles.hh>

//Alglib:
#include "apps/pilot/vmullig/alglib/ap.h"
#include "apps/pilot/vmullig/alglib/alglibmisc.h"
#include "apps/pilot/vmullig/alglib/alglibinternal.h"
#include "apps/pilot/vmullig/alglib/dataanalysis.h"
#include "apps/pilot/vmullig/alglib/linalg.h"
#include "apps/pilot/vmullig/alglib/optimization.h"
#include "apps/pilot/vmullig/alglib/solvers.h"
#include "apps/pilot/vmullig/alglib/statistics.h"
#include "apps/pilot/vmullig/alglib/specialfunctions.h"
#include "apps/pilot/vmullig/alglib/ap.cpp"
#include "apps/pilot/vmullig/alglib/alglibmisc.cpp"
#include "apps/pilot/vmullig/alglib/alglibinternal.cpp"
#include "apps/pilot/vmullig/alglib/dataanalysis.cpp"
#include "apps/pilot/vmullig/alglib/linalg.cpp"
#include "apps/pilot/vmullig/alglib/optimization.cpp"
#include "apps/pilot/vmullig/alglib/solvers.cpp"
#include "apps/pilot/vmullig/alglib/statistics.cpp"
#include "apps/pilot/vmullig/alglib/specialfunctions.cpp"

//For silent file output:
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/simple_moves/ConstraintSetMover.hh>

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

using ObjexxFCL::format::F;
using namespace protocols::cluster;

OPT_KEY( Boolean, v_prerelax ) //Should imported structers be subjected to a round of fast relaxation?
OPT_KEY( Integer, v_relaxrounds ) //The number of fastrelax rounds to apply (default 1)
OPT_KEY( String, v_clusterby ) //The basis for clustering.  bb_cartesian or bb_dihedral
OPT_KEY( Real, v_clusterradius ) //The cluster radius.
OPT_KEY( Boolean, v_weightbyenergy ) //Weight by energy when calculating cluster centers?  True by default.
OPT_KEY( Real, v_kbt ) //The value of k_B*T to use when weighting by energy.  1.0 by default.
OPT_KEY( Boolean, v_CB ) //Should beta carbons be used in Cartesian clustering?
OPT_KEY( IntegerVector, v_ignoreresidue) //Residues to ignore in alignments.
OPT_KEY( Integer, v_limit_structures_per_cluster) //Structures to output per cluster.
OPT_KEY( Boolean, v_cyclic) //Is there a peptide bond between the N- and C-termini?  False by default.
OPT_KEY( Boolean, v_cluster_cyclic_permutations) //Should cyclic permutations be clustered together?  False by default.
OPT_KEY( Integer, v_cyclic_permutation_offset) //1 by default, meaning that every cyclic permutation is clustered if v_cluster_cyclic_permutations is true.  Values X > 1 mean that cyclic permutations shifted by X residues will be clustered.
OPT_KEY( Boolean, v_mutate_to_ala) //Mutate the output structure to a chain of alanines?
OPT_KEY( IntegerVector, v_disulfide_positions) //Positions that are disulfide-bonded.
OPT_KEY( Boolean, v_homooligomer_swap) //Try swapping homoolgigomer chains when calculating RMSDs?
OPT_KEY( Boolean, v_silentoutput) //Use silent files as output?  Default false.
OPT_KEY( File, v_cst_file) //Constraint file.  Default unused.
OPT_KEY( StringVector, v_extra_rms_atoms) //Extra atoms to use in the RMSD calculation.

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1<core::Size> empty_vector;

	NEW_OPT ( v_prerelax, "Should imported structures be subjected to a round of fast relaxation?  Default false.", false);
	NEW_OPT ( v_relaxrounds, "The number of fastrelax rounds to apply.  Default 1.", 1);
	NEW_OPT ( v_clusterby, "What should I use as the basis for clustering?  Options are \"bb_cartesian\" (xyz coordinates of backbone atoms) and \"bb_dihedral\" (phi, psi, omega angles).  Default \"bb_cartesian\".", "bb_cartesian");
	NEW_OPT ( v_clusterradius, "The radius for clustering, in Angstroms for Cartesian clustering and degrees for dihedral clustering.  Default 1.0.", 1.0 );
	NEW_OPT ( v_weightbyenergy, "Should structures be weighted by exp(-E/(k_B*T)) when calculating cluster centers?  True by default.", true);
	NEW_OPT ( v_kbt, "The value of k_b*T to use when weighting by exp(-E/(k_B*T)).  1.0 by default.", 1.0);
	NEW_OPT ( v_CB, "If clustering by backbone Cartesian coordinates, should beta carbons be included?  Default true.", true);
	NEW_OPT ( v_ignoreresidue, "List of residues to ignore in alignments for clustering.  Default empty list.", empty_vector);
	NEW_OPT ( v_limit_structures_per_cluster, "Maximum number of structures to output per cluster.  Default no limit.", 0);
	NEW_OPT ( v_cyclic, "If true, constraints are added to make a peptide bond between the N- and C-termini.  If false (default), the termini are free.  Default false.", false);
	NEW_OPT ( v_cluster_cyclic_permutations, "If true, all cyclic permutations are tried when comparing two structures for clustering.  Requires \"-v_cyclic true\".  Default false.", false);
	NEW_OPT ( v_cyclic_permutation_offset, "1 by default, meaning that every cyclic permutation is clustered if v_cluster_cyclic_permutations is true.  Values X > 1 mean that cyclic permutations shifted by X residues will be clustered.", 1);
	NEW_OPT ( v_mutate_to_ala, "If true, the input structures will be converted to a chain of alanines (L- or D-) before scoring.  Default false.", false);
	NEW_OPT ( v_disulfide_positions, "A space-separated list of positions that are disulfide-bonded.  For example, -v_disulfide_positions 3 8 6 23 would mean that residues 3 and 8 are disulfide-bonded, as are residues 6 and 23.  Default empty list.", empty_vector);
	NEW_OPT ( v_homooligomer_swap, "If the structures contain multiple chains with identical sequence, setting this to true will test all permutations of chains when clustering.  Default false.", false);
	NEW_OPT ( v_silentoutput, "Write output to a silent file instead of to separate PDBs.  This will create two files: one that only contains the first member of each cluster, and one that contains everything.  Default false.", false);
	NEW_OPT ( v_cst_file, "A user-specified constraints file.  Default unused.", "" );
	NEW_OPT ( v_extra_rms_atoms, "A list of additional atoms to use in the RMSD calculation, in the format \"residue:atomname residue:atomname residue:atomname\".  For example, \"-v_extra_rms_atoms 7:SG 12:CG 12:CD 12:CE 12:NZ 14:OG\".  Default empty list.", "");
}

using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;
using namespace basic::options;

void parse_extra_atom_list ( utility::vector1 < core::id::NamedAtomID > &extra_atom_list )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(!option[v_extra_rms_atoms].user()) return; //Do nothing if there's no list of extra atoms.
	extra_atom_list.clear();

	utility::vector1 <std::string> tags = option[v_extra_rms_atoms]();

	//Parse each tag:
	printf("The following additional atoms will be used in the RMSD calculation:\n");
	for(core::Size i=1, imax=tags.size(); i<=imax; ++i) { //Loop through each tag
		core::Size colonposition = tags[i].find( ':' );
		std::string resstring = tags[i].substr( 0, colonposition );
		core::Size res = static_cast<core::Size>( atoi( resstring.c_str() ) ); //The residue number
		std::string atomname = tags[i].substr( colonposition + 1); //The atom name
		printf("\tResidue %lu, Atom %s\n", res, atomname.c_str());
		extra_atom_list.push_back( core::id::NamedAtomID( atomname, res ) );
	}

	printf("\n");

	return;
}

bool
use_in_rmsd_offset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno2,
	core::Size atomno,
	const core::Size pose1_offset,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	signed int resno=resno2-pose1_offset;
	if(resno<=0) resno+=pose1.n_residue();

	if(option[v_ignoreresidue].size()>0) {
		for(core::Size i=1; i<=option[v_ignoreresidue].size(); i++) if (resno2 == static_cast<core::Size>(option[v_ignoreresidue][i])) return false; //If this residue is to be ignored, return false.
	}

	if(pose1.residue(resno).has( "N") && pose2.residue(resno2).has( "N") && pose1.residue(resno).atom_index( "N")==atomno) return true;
	if(pose1.residue(resno).has("CA") && pose2.residue(resno2).has("CA") && pose1.residue(resno).atom_index("CA")==atomno) return true;
	if(pose1.residue(resno).has( "C") && pose2.residue(resno2).has( "C") && pose1.residue(resno).atom_index( "C")==atomno) return true;
	if(pose1.residue(resno).has( "O") && pose2.residue(resno2).has( "O") && pose1.residue(resno).atom_index( "O")==atomno
			&& !pose1.residue(resno).is_upper_terminus() &&!pose2.residue(resno2).is_upper_terminus()) return true;
	if(pose1.residue(resno).has("CM") && pose2.residue(resno2).has("CM") && pose1.residue(resno).atom_index("CM")==atomno) return true;
	if(option[v_CB]()) {
		if(pose1.residue(resno).has( "CB") && pose2.residue(resno2).has( "CB") && pose1.residue(resno).atom_index( "CB")==atomno) return true;
	}

	//Check extra atoms:
	for(core::Size i=1, imax=extra_atom_list.size(); i<=imax; ++i) {
		if(resno!=static_cast<signed int>(extra_atom_list[i].rsd())) continue;
		if(pose1.residue(resno).has( extra_atom_list[i].atom() ) && pose2.residue(resno2).has( extra_atom_list[i].atom() ) && pose1.residue(resno).atom_index( extra_atom_list[i].atom() )==atomno) return true;
	}

	return false;
}

bool
use_in_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno,
	core::Size atomno,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(option[v_ignoreresidue].size()>0) {
		for(core::Size i=1; i<=option[v_ignoreresidue].size(); i++) if (resno == static_cast<core::Size>(option[v_ignoreresidue][i])) return false; //If this residue is to be ignored, return false.
	}

	if(pose1.residue(resno).has( "N") && pose2.residue(resno).has( "N") && pose1.residue(resno).atom_index( "N")==atomno) return true;
	if(pose1.residue(resno).has("CA") && pose2.residue(resno).has("CA") && pose1.residue(resno).atom_index("CA")==atomno) return true;
	if(pose1.residue(resno).has( "C") && pose2.residue(resno).has( "C") && pose1.residue(resno).atom_index( "C")==atomno) return true;
	if(pose1.residue(resno).has( "O") && pose2.residue(resno).has( "O") && pose1.residue(resno).atom_index( "O")==atomno
			&& !pose1.residue(resno).is_upper_terminus() &&!pose2.residue(resno).is_upper_terminus()) return true;
	if(pose1.residue(resno).has("CM") && pose2.residue(resno).has("CM") && pose1.residue(resno).atom_index("CM")==atomno) return true;
	if(option[v_CB]()) {
		if(pose1.residue(resno).has( "CB") && pose2.residue(resno).has( "CB") && pose1.residue(resno).atom_index( "CB")==atomno) return true;
	}

	//Check extra atoms:
	for(core::Size i=1, imax=extra_atom_list.size(); i<=imax; ++i) {
		if(resno!=extra_atom_list[i].rsd()) continue;
		if(pose1.residue(resno).has( extra_atom_list[i].atom() ) && pose2.residue(resno).has( extra_atom_list[i].atom() ) && pose1.residue(resno).atom_index( extra_atom_list[i].atom() )==atomno) return true;
	}

	return false;
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

//Function to check whether this is a beta-amino acid:
bool
is_beta_aminoacid(
	const core::conformation::Residue &rsd
) {
	return rsd.type().is_beta_aa();
	//if (rsd.has("CM")) return true; //For now, this is the test for a beta-residue.
	//return false; //False otherwise.
}

//Variant of the above function that takes a pose and a position.
bool
is_beta_aminoacid(
	const core::pose::Pose &mypose,
	const core::Size position
) {
	return is_beta_aminoacid(mypose.residue(position));
}

//Function to count the number of backbone dihedral angles to be stored:
core::Size count_bb_dihedrals (const core::pose::Pose &mypose)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size ndihedrals = 0;
	for(core::Size ir=1, nres=mypose.n_residue(); ir<=nres; ir++) {
			if(!mypose.residue(ir).type().is_alpha_aa() && !mypose.residue(ir).type().is_beta_aa()) continue;
			ndihedrals+=3;
			if(mypose.residue(ir).is_lower_terminus() || ir==1) ndihedrals--;
			if(mypose.residue(ir).is_upper_terminus() || ir==nres) ndihedrals-=2;
			if(is_beta_aminoacid(mypose, ir)) ndihedrals++;
		}
		if(option[v_cyclic]()) ndihedrals+=3; //There are three more backbone dihedral angles if this is a cyclic peptide
		return(ndihedrals);
}

//Function to determine whether a value is in a list
bool is_in_list (
	const core::Size val,
	const utility::vector1 <core::Size> &vallist
) {
	if(vallist.size()>0) {
		for(core::Size i=1, listlength=vallist.size(); i<=listlength;i++){
			if (vallist[i]==val) return true;
		}
	}
	return false;
}

//Overloaded variant of below function to ensure that every angle lies between -180 and 180, ignoring certain entries
void check_angle_bounds(
	utility::vector1 <core::Real> anglevector,
	utility::vector1 <core::Size> skiplist //list of entries to ignore
) {
	for(core::Size i=1; i<=anglevector.size(); i++) {
		if(!is_in_list(i,skiplist)) {
			if(anglevector[i] > 180.0) anglevector[i] = fmod(anglevector[i]+180.0,360.0)-180.0;
			if(anglevector[i] < -180.0) anglevector[i] = -(fmod(-anglevector[i]+180.0,360.0)-180.0);
		}
		//if(anglevector[i]<180.0 && anglevector[i]>-180.0) continue;
		//anglevector[i]=fmod(anglevector[i]+180.0, 360.0)-180.0; //THIS ISN'T RIGHT.
	}
	
	return;
}

//Function to ensure that every angle lies between -180 and 180:
void check_angle_bounds(
	utility::vector1 <core::Real> anglevector
) {
	utility::vector1 <core::Size> skiplist;
	check_angle_bounds(anglevector, skiplist);
	
	return;
}


//Function to set disulfides:
void make_disulfides (core::pose::Pose outputpose) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(option[v_disulfide_positions].user() && option[v_disulfide_positions]().size()>0) {
		utility::vector1 < std::pair < core::Size, core::Size > > disulfpairs;
		for(core::Size i=1; i<=option[v_disulfide_positions]().size(); i+=2) 
			disulfpairs.push_back(std::pair< core::Size, core::Size > ( option[v_disulfide_positions]()[i], option[v_disulfide_positions]()[i+1] ));
		outputpose.conformation().fix_disulfides(disulfpairs);
	} else {
		outputpose.conformation().detect_disulfides();
	}
	return;
}

void storeposedata(
	const core::pose::Pose &pose,
	utility::vector1 < core::Real > &posedata,
	FArray2D < core::Real > &alignmentdata, //Only for Cartesian clustering: x,y,z coordinates of atoms to be used for alignment.
	const core::Size clustermode,	
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	posedata.clear();
	const core::Size nres = pose.n_residue();

	//Count the number of atoms in the alignment data:
	if(clustermode==1) {
		core::Size alignmentdatasize = 0;
		for(core::Size ir=1; ir<=nres; ir++) {
			if(is_in_list(ir, option[v_ignoreresidue]())) continue; //Skip this residue if it's in the ignore list.
			if(!pose.residue(ir).type().is_alpha_aa() && !pose.residue(ir).type().is_beta_aa()) continue; //Skip this if not an alpha or beta amino acid
			for(core::Size ia=1,iamax=pose.residue(ir).natoms(); ia<=iamax; ia++) {
				if(use_in_rmsd(pose,pose,ir,ia, extra_atom_list)) alignmentdatasize++;
			}
		}
		alignmentdata.redimension(3, alignmentdatasize);
	}
	

/*	//I might want to change the following to something that just counts atoms to be stored...	
	core::Size nglycines=0;
	if(option[v_CB]()) { //Count glycines (and anything else that lacks a CB) if CBs are used:
		for(core::Size ir=1; ir<=pose.n_residue(); ir++)
			if(!pose.residue(ir).has("CB")) nglycines++;
	}

	//Count beta-amino acids:
	core::Size nbeta=0;
	for(core::Size ir=1, nres=pose.n_residue(); ir<=nres; ir++)
		if(is_beta_aminoacid(pose, ir)) nbeta++;

	if(clustermode==1) alignmentdata.redimension(3, (pose.n_residue()-option[v_ignoreresidue].size())*(option[v_CB]()?5:4) - nglycines*(option[v_CB]()?1:0) + nbeta); //Assumes that the user hasn't specified the same residue multiple times!
*/

	core::Size icount=0;

	for(core::Size ir=1, nres=pose.n_residue(); ir<=nres; ir++) { //Loop through all residues
		switch(clustermode) { //Depending on the clustering mode do different things
			case 1: //backbone Cartesian coordinate clustering
				for(core::Size ia=1, natom=pose.residue(ir).natoms(); ia<=natom; ia++) {
					posedata.push_back(pose.residue(ir).atom(ia).xyz()[0]);
					posedata.push_back(pose.residue(ir).atom(ia).xyz()[1]);
					posedata.push_back(pose.residue(ir).atom(ia).xyz()[2]);
					bool useresidue=true;
					if(!pose.residue(ir).type().is_alpha_aa() && !pose.residue(ir).type().is_beta_aa()) useresidue=false; //Skip if this isn't an alpha or beta amino acid
					else {
						if(option[v_ignoreresidue].size() > 0) {
							for(core::Size i=1; i<=option[v_ignoreresidue].size(); i++) {
								if(ir==static_cast<core::Size>(option[v_ignoreresidue][i])) {
									useresidue=false;
									break;
								}
							}
						}
					}
					if(useresidue && use_in_rmsd(pose, pose, ir, ia, extra_atom_list)) { //Store data for Cartesian alignment:
						icount++;
						alignmentdata(1, icount)=pose.residue(ir).atom(ia).xyz()[0];
						alignmentdata(2, icount)=pose.residue(ir).atom(ia).xyz()[1];
						alignmentdata(3, icount)=pose.residue(ir).atom(ia).xyz()[2];
						//printf("%f\t%f\t%f\n", pose.residue(ir).atom(ia).xyz()[0],pose.residue(ir).atom(ia).xyz()[1],pose.residue(ir).atom(ia).xyz()[2]); fflush(stdout); //DELETE ME
					}
				}
				break;
			case 2: //backbone dihedral clustering
				if(is_beta_aminoacid(pose, ir)) { //Beta amino acids:
					if(ir>1 && !pose.residue(ir).is_lower_terminus()) posedata.push_back(pose.residue(ir).mainchain_torsion(1)); //Phi for a beta-amino acid
					posedata.push_back(pose.residue(ir).mainchain_torsion(2)); //Theta for a beta-amino acid
					if(ir<nres && !pose.residue(ir).is_upper_terminus()) {
						posedata.push_back(pose.residue(ir).mainchain_torsion(3)); //Psi for a beta-amino acid
						posedata.push_back(pose.residue(ir).mainchain_torsion(4)); //Omega for a beta-amino acid
					}
				} else if (pose.residue(ir).type().is_alpha_aa()) { //Alpha amino acids:
					if(ir>1 && !pose.residue(ir).is_lower_terminus()) posedata.push_back(pose.phi(ir));
					if(ir<nres && !pose.residue(ir).is_upper_terminus()) {
						posedata.push_back(pose.psi(ir));
						posedata.push_back(pose.omega(ir));
					}
				} else {} //Default case -- neither alpha nor beta amino acid, in which case we do nothing.
				break;
			default:
				break;
		}
	}

	if (clustermode==2) { //Additional data to be stored in the backbone dihedral clustering case:
		if(option[v_cyclic]()) { //Store the additional dihedral angles (psi_n, omega_n->1, phi_1) for backbone clustering
			posedata.push_back (numeric::dihedral_degrees( //Psi of last residue
				pose.residue(pose.n_residue()).xyz( (is_beta_aminoacid(pose,pose.n_residue())?"CA":"N") ),
				pose.residue(pose.n_residue()).xyz( (is_beta_aminoacid(pose,pose.n_residue())?"CM":"CA") ),
				pose.residue(pose.n_residue()).xyz("C"),
				pose.residue(1).xyz("N")
			));
			posedata.push_back (numeric::dihedral_degrees( //Omega of peptide bond between last and first residues
				pose.residue(pose.n_residue()).xyz( (is_beta_aminoacid(pose,pose.n_residue())?"CM":"CA") ),
				pose.residue(pose.n_residue()).xyz("C"),
				pose.residue(1).xyz("N"),
				pose.residue(1).xyz("CA")
			));
			posedata.push_back (numeric::dihedral_degrees( //Phi of first residue
				pose.residue(pose.n_residue()).xyz("C"),
				pose.residue(1).xyz("N"),
				pose.residue(1).xyz("CA"),
				pose.residue(1).xyz( (is_beta_aminoacid(pose,1)?"CM":"C") )
			));
		}
		for (core::Size ir=1; ir<=pose.n_residue(); ir++) { //Store side chain dihedrals, too, for rebuilding structures later
			if(pose.residue(ir).nchi()==0) continue;
			for(core::Size ichi=1; ichi<=pose.residue(ir).nchi(); ichi++) {
				posedata.push_back(pose.residue(ir).chi(ichi));
				//printf("Stored residue %lu chi %lu as posedata index %lu.\n", ir, ichi, posedata.size()); fflush(stdout); //DELETE ME
			}
		}
		check_angle_bounds(posedata);
	}

	return;
}

//This function is overloaded (see above).
void storeposedata(
	const core::pose::Pose &pose,
	utility::vector1 < utility::vector1 <core::Real> > &posedata,
	utility::vector1 < FArray2D < core::Real > > &alignmentdata,
	const core::Size clustermode,	
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	utility::vector1 <core::Real> posedata_tostore;
	FArray2D < core::Real > alignmentdata_tostore;
	storeposedata(pose, posedata_tostore, alignmentdata_tostore, clustermode, extra_atom_list);
	posedata.push_back(posedata_tostore); //Not as efficient as it could be, but that's okay.
	alignmentdata.push_back(alignmentdata_tostore); //Stores an empty FArray2D if not clustermode 1.
	return;
}

//Function to return the difference between two values.  If these are angles, it figures out the angular difference.
//This assumes that all angles are between -180.0 and 180.0.
//UPDATE THIS IF ADDITIONAL CLUSTER MODES ARE ADDED.
core::Real elementdiff (
	const core::Real &n1,
	const core::Real &n2,
	const core::Size &clustmode
) {
	if (clustmode==1) return n1-n2;
	else if (clustmode==2) {
		core::Real diff = n1-n2;
		if (diff<-180.0) diff+=360.0;
		else if (diff>180.0) diff-=360.0;
		return diff;
		//return (core::Real) (fmod((double)n1-(double)n2 + 180.0, 360.0)-180.0); //RETURN TO THIS!  I don't think it's right.
	}
	return 0.0;
}

//Function to create a pose from posedata
void pose_from_posedata (
	const core::pose::Pose &inputpose,
	core::pose::Pose &outputpose,
	const core::Size clustermode,
	const utility::vector1<core::Real> &posedata,
	const bool set_chi=true
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	outputpose=inputpose; //Copy the input pose.
	make_disulfides(outputpose);

	core::Size counter=0;
	numeric::xyzVector <core::Real> atomxyz;
	switch(clustermode) {
		case 1:
			for(core::Size ir=1; ir<=outputpose.n_residue(); ir++) {
				for(core::Size ia=1; ia<=outputpose.residue(ir).natoms(); ia++) {
					counter+=3;
					atomxyz[0]=posedata[counter-2];
					atomxyz[1]=posedata[counter-1];
					atomxyz[2]=posedata[counter];
					outputpose.conformation().set_xyz(id::AtomID(ia, ir), atomxyz);
				}
			}
			break;
		case 2: //Dihedral clustering case.
			counter=1;
			for(core::Size i=1, nres=outputpose.n_residue(); i<=nres; i++) {
				if(is_beta_aminoacid(outputpose, i)) {//beta-amino acids:
					if(i>1 && !inputpose.residue(i).is_lower_terminus()) betapeptide_setphi(outputpose, i, posedata[counter++]);
					betapeptide_settheta(outputpose, i, posedata[counter++]);
					if(i<nres && !inputpose.residue(i).is_upper_terminus()) {
						betapeptide_setpsi(outputpose, i, posedata[counter++]);
						betapeptide_setomega(outputpose, i, posedata[counter++]);
					}
				}else if(outputpose.residue(i).type().is_alpha_aa()){ //alpha-amino acids:
					if(i>1 && !inputpose.residue(i).is_lower_terminus()) outputpose.set_phi(i, posedata[counter++]);
					if(i<nres && !inputpose.residue(i).is_upper_terminus()) {
						outputpose.set_psi(i, posedata[counter++]);
						outputpose.set_omega(i, posedata[counter++]);
					}
				} else {} //Other residues -- do nothing, for now.
			}

			{	//Also need to set side-chain dihedrals:
				if(option[v_cyclic]()) {
					counter+=3;
					outputpose.conformation().set_torsion_angle(//Set psi of last residue
						core::id::AtomID(outputpose.residue(outputpose.n_residue()).atom_index( (is_beta_aminoacid(outputpose,outputpose.n_residue())?"CA":"N") ), outputpose.n_residue()),
						core::id::AtomID(outputpose.residue(outputpose.n_residue()).atom_index( (is_beta_aminoacid(outputpose,outputpose.n_residue())?"CM":"CA") ), outputpose.n_residue()),
						core::id::AtomID(outputpose.residue(outputpose.n_residue()).atom_index("C"), outputpose.n_residue()),
						core::id::AtomID(outputpose.residue(outputpose.n_residue()).atom_index("O"), outputpose.n_residue()),
						(posedata[counter-3]-180.0)/180.0*PI);
					outputpose.conformation().set_torsion_angle( //Set phi of first residue
						core::id::AtomID(outputpose.residue(1).atom_index("H"), 1),
						core::id::AtomID(outputpose.residue(1).atom_index("N"), 1),
						core::id::AtomID(outputpose.residue(1).atom_index("CA"), 1),
						core::id::AtomID(outputpose.residue(1).atom_index( (is_beta_aminoacid(outputpose,1)?"CM":"C") ), 1),
						(posedata[counter-1]-180.0)/180.0*PI);
				}

				core::Size ir=1;
				core::Size ichi=0;
				core::Size ijump=0;
				if(set_chi) {
					//printf("posedata.size()=%lu\n", posedata.size()); fflush(stdout); //DELETE ME
					while(counter<=posedata.size()) {
						ichi++;
						//printf("ir=%lu ichi=%lu counter=%lu\n", ir, ichi, counter); fflush(stdout); //DELETE ME
						if(ichi>outputpose.residue(ir).nchi()) {//If we're done with this residue, go on to the next.
							ichi=0;
							ir++;
							continue;
						}

						outputpose.set_chi(ichi, ir, posedata[counter++]);
						if(counter>posedata.size()) break; //Should be redundant.
					}
				}
				while(counter<=posedata.size()) { //Need to set jumps
					ijump++;
					numeric::xyzVector < core::Real > transvect = outputpose.jump(ijump).get_translation();
					numeric::EulerAngles < core::Real > euler_angles(outputpose.jump(ijump).get_rotation());
					for(core::Size j=0; j<3; j++) transvect[j]=posedata[counter++];
					euler_angles.phi_degrees( posedata[counter++] );
					euler_angles.theta_degrees( posedata[counter++] );
					euler_angles.psi_degrees( posedata[counter++] );
					//Update the jump:
					core::kinematics::Jump myjump = outputpose.jump(ijump);
					myjump.set_translation( transvect);
					myjump.set_rotation( euler_angles.to_rotation_matrix() );
					outputpose.set_jump(ijump, myjump);
				}
			}
			break;
		default:
			break;
	}

	make_disulfides(outputpose);
	outputpose.update_residue_neighbors();

	return;
}

//Function to return the number of alignment atoms in a given residue:
core::Size alignment_atoms_in_res (
	const core::pose::Pose &pose,
	const core::Size position
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(!pose.residue(position).type().is_alpha_aa() && !pose.residue(position).type().is_beta_aa()) return 0;

	core::Size numatoms = 4;
	if(pose.residue(position).is_upper_terminus()) numatoms--; //"O" not included if upper terminus
	if(option[v_CB]() && pose.residue(position).has("CB")) numatoms++;
	if(is_beta_aminoacid(pose, position)) numatoms++;
	return numatoms;
}

//Function to return the number of alignment torsion angles stored for a given residue:
core::Size alignment_torsions_in_res (
	const core::pose::Pose &pose,
	const core::Size position
) {
	if(!pose.residue(position).type().is_alpha_aa() && !pose.residue(position).type().is_beta_aa()) return 0;

	core::Size numtors = 3;
	if(is_beta_aminoacid(pose, position)) numtors++;
	if(position==1) numtors--; //No phi for first position.
	if(position==pose.n_residue()) numtors-=2; //No psi or omega for last position.
	return numtors;
}

//Function to calculate the factorial of a number (recursively calls itself):
core::Size factorial(core::Size val) 
{   if (val == 0) return 1;
    return val * factorial(val - 1);
}

//Function to get the first index in the alignment vector (vector of x,y,z atom coordinates) corresponding to an amino acid:
core::Size get_start_index (const core::Size aa_number, const core::pose::Pose &refpose)
{
	core::Size counter = 1;
	for(core::Size ir=1, nres=refpose.n_residue(); ir<=nres; ir++)	{
		if(ir == aa_number) return counter;
		counter += alignment_atoms_in_res(refpose,ir);
	}
	return 0;
}

//Function to get the first amino acid in a chain:
core::Size get_start_aa (const core::Size chainnumber, const core::pose::Pose &refpose)
{
	for(core::Size ir=1, nres=refpose.n_residue(); ir<=nres; ir++) if(static_cast<core::Size>(refpose.chain(ir)) == chainnumber) return ir;
	return 0;
}

//Function to swap chains in a pose
void swapchains(
	core::pose::Pose &currentpose, //Input and output
	const core::Size permutation_number //Input -- the chain perturbation number
) {
	if(permutation_number == 1) return; //Do nothing if this is the first perturbation (no rearrangement of the list of chains).

	const core::Size chaincount = currentpose.num_jump()+1;
	if(chaincount < 2) return; //Do nothing if there's only one chain.

	utility::vector1 < core::Size > chainlist; //List of chains -- order will be perturbed.
	utility::vector1 < core::Size > startaa_list; //List of starting amino acids of chains
	utility::vector1 < core::Size > startindex_list; //List of the starting index in the parentvect for each chain
	chainlist.resize(chaincount);
	startaa_list.resize(chaincount);
	for(core::Size i=1; i<=chainlist.size(); i++) chainlist[i]=i; //Initialize the list

	//Reorder the list of chains:
	if(permutation_number>1) {
		for(core::Size i=2; i<=permutation_number; i++) {
			std::next_permutation(chainlist.begin(), chainlist.end());
		}
	}

	//DELETE THE FOLLOWING -- only for testing!
	//printf("\nPermutation:\t"); //DELETE ME
	//for(core::Size i=1; i<=chainlist.size(); i++) printf("%lu\t", chainlist[i]); //DELETE ME

	//Get the starting aa for each chain:
	for(core::Size i=1; i<=chainlist.size(); i++) startaa_list[i]=get_start_aa(chainlist[i], currentpose);

	//A pose into which residues will be copied:
	core::pose::Pose newpose;

	//Copy the pose, swapping chains around:
	for(core::Size i=1; i<=chainlist.size(); i++) {
		for (core::Size ir=startaa_list[i], nres=currentpose.n_residue(); ir<=nres && static_cast<core::Size>(currentpose.chain(ir)) == chainlist[i]; ir++)
		{
			if(ir==startaa_list[i]) {
				newpose.append_residue_by_jump(currentpose.residue(ir), 1, "", "", true);
			} else {
				newpose.append_residue_by_bond(currentpose.residue(ir));
			}
		}
	}

	newpose.update_residue_neighbors();

	//newpose.dump_pdb("temp.pdb"); //DELETE ME!

	currentpose=newpose;

	return;
}

//Function to swap chains in an alignment vector (vector of x,y,z atom coordinates):
void swapalignmentvector (
	const FArray2D < core::Real > &parentvect,
	FArray2D < core::Real > &swappedvect,
	const core::pose::Pose &refpose,
	core::Size permutation_number
) {
	swappedvect = parentvect; //Copy the parent vector

	utility::vector1 < core::Size > chainlist; //List of chains, which will be rearranged in each permutation's order
	utility::vector1 < core::Size > startaa_list; //List of starting amino acids of chains
	utility::vector1 < core::Size > startindex_list; //List of the starting index in the parentvect for each chain

	for(core::Size i=1, imax=refpose.num_jump()+1; i<=imax; i++) chainlist.push_back(i);
	if(permutation_number>1) {
		for(core::Size i=2; i<=permutation_number; i++) {
			std::next_permutation(chainlist.begin(), chainlist.end());
		}
	}

	//Store the starting amino acid numbers for each chain, in the perturbed chain order:
	startaa_list.resize(chainlist.size());
	for(core::Size i=1; i<=startaa_list.size(); i++) startaa_list[i]=get_start_aa(chainlist[i], refpose);

	//Store the starting indices for each chain, in the perturbed chain order:
	startindex_list.resize(chainlist.size());
	for(core::Size i=1; i<=startindex_list.size(); i++) startindex_list[i]=get_start_index(startaa_list[i], refpose);

	//DELETE THE FOLLOWING!  ONLY FOR TESTING!
	//printf("\nPermutation:\t"); //DELETE ME
	//for(core::Size i=1; i<=chainlist.size(); i++) printf("%lu\t", chainlist[i]); //DELETE ME
	//printf("\n"); //DELETE ME
	//printf("Starting_aa:\t"); //DELETE ME
	//for(core::Size i=1; i<=startaa_list.size(); i++) printf("%lu\t", startaa_list[i]); //DELETE ME
	//printf("\n"); //DELETE ME
	//printf("Start_index:\t"); //DELETE ME
	//for(core::Size i=1; i<=startindex_list.size(); i++) printf("%lu\t", startindex_list[i]); //DELETE ME
	//printf("\n\n"); //DELETE ME
	//fflush(stdout); //DELETE ME

	core::Size counter = 1;
	for(core::Size ichain=1, nchain=startindex_list.size(); ichain<=nchain; ichain++) {
		core::Size counter2 = startindex_list[ichain]; //Starting index in the parent vector for this chain
		for(core::Size ir=startaa_list[ichain], nres=refpose.n_residue(); ir<=nres && static_cast<core::Size>(refpose.chain(ir))==chainlist[ichain]; ir++) { //Loop through all residues in the current chain
			for(core::Size iatom=1, natom=alignment_atoms_in_res(refpose, ir); iatom<=natom; iatom++) { //Loop through all alignment atoms in the current residue
				swappedvect(1, counter) = parentvect (1, counter2); //X
				swappedvect(2, counter) = parentvect (2, counter2); //Y
				swappedvect(3, counter) = parentvect (3, counter2); //Z
				counter++; counter2++; //Increment both counters to point at the next atom.
			}
		}
	}

	return;
}

//Function to swap chains in a dihedral vector
void swapvector (
	const utility::vector1 < core::Real > &parentvect,
	utility::vector1 < core::Real > &swappedvect,
	const core::pose::Pose &/*refpose*/,
	core::Size /*permutation_number*/
) {

	//TODO!!!
	//Currently there is a check to ensure that this function is never called.

	swappedvect=parentvect;

	return;
}

//Function to calculate the RMSD between two poses, based on whatever is in "posedata" for the second only.
//This assumes that a pose already exists for the first, and that this is provided in refpose.
core::Real calc_dist(
	const utility::vector1 <core::Real> &vect1, //Backbone dihedral vector 1
	const utility::vector1 <core::Real> &vect2, //Backbone dihedral vector 2
	const core::Size clustmode,
	const FArray2D < core::Real > &alignmentvect1, //Alignment vector 1 (for Cartesian clustering)
	const FArray2D < core::Real > &alignmentvect2, //Alignment vector 2 (for Cartesian clustering)
	const core::Size nresidues,
	const core::pose::Pose &firstpose //Used for reference only!
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Real accumulator = 0.0;

	if(alignmentvect1.size2()!=alignmentvect2.size2()) {
		printf("Internal error!  Alignment vector size mismatch!  Crashing!!!\n"); fflush(stdout);
		exit(1);
	}

	if(clustmode==1) accumulator=numeric::model_quality::rms_wrapper(alignmentvect1.size2(), alignmentvect1, alignmentvect2); //Calculate the rms, if we're doing Cartesian clustering
	else if (clustmode==2) { //Backbone dihedral clustering
		bool includeme=true;
		for(core::Size ir=1,index=1; ir<=nresidues; ir++) {
			//The number of torsion angles in the current residue:
			core::Size curres_torsioncount = 3;
			if(firstpose.residue(ir).is_lower_terminus() || ir==1) curres_torsioncount--;
			if(firstpose.residue(ir).is_upper_terminus() || ir==nresidues) curres_torsioncount-=2;
			if(is_beta_aminoacid(firstpose, ir)) curres_torsioncount++;

			if(option[v_ignoreresidue].size()>0) {
				includeme=true;
				for(core::Size j=1; j<=(core::Size)option[v_ignoreresidue].size(); j++) {
					if( ir == static_cast<core::Size>(option[v_ignoreresidue][j])) {
						includeme=false;
						index+=curres_torsioncount; //Increment the dihedral angle index to skip this residue.
						break;
					}
				}
			}

			if(includeme) { //If the current residue, ir, is to be included in the rms calculation

				for(core::Size i=1; i<=curres_torsioncount; i++) {
					accumulator+=pow(elementdiff(vect1[index],vect2[index], clustmode), 2.0);
					index++;				
				}

				if(option[v_cyclic]() && ir==nresidues) { //If this is a backbone-cyclized peptide, then there are additional dihedral angles stored in vect1 and vect2.
					for(core::Size j=1; j<=3; j++) accumulator+=pow(elementdiff(vect1[index],vect2[index], clustmode), 2.0);
				}
			}
		}
		accumulator=sqrt(accumulator);
	}

	return accumulator;
}

//Overloaded form of calc_dist for cyclic permutations and for homooligomer swapping:
core::Real calc_dist(
	const utility::vector1 <core::Real> &vect1,
	const utility::vector1 <core::Real> &vect2,
	const core::Size clustmode,
	const FArray2D < core::Real > &alignmentvect1,
	const FArray2D < core::Real > &alignmentvect2,
	const core::Size nresidues,
	const core::pose::Pose &firstpose, //Used for reference only
	core::Size &offset, //Used only for cyclic permutations (OUTPUT)
	core::Size &permutation //Used only for homooligomer permutations (OUTPUT)
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Real dist=0.0;

	if(option[v_cluster_cyclic_permutations]()) {
		if(clustmode==1) { //Clustering by backbone atom Cartesian coordinates
			//const core::Size atomsperresidue = (option[v_CB]() ? 5 : 4) ;
			FArray2D < core::Real > avect2 (3, alignmentvect2.size2());
			for(core::Size ioffset=0; ioffset<nresidues; ioffset+=option[v_cyclic_permutation_offset]()) { //Loop through all possible offsets
				core::Size offset_index_1 = 1; //The index of the first atom (i.e. alignment atom 1 of the OLD first residue), when offset by ioffset residues.
				//Figure out starting value of offset_index_1.  This means that we're putting the atoms from the END of the sequence ahead of the first atom:
				if(ioffset>0) {
					for(core::Size ir=nresidues; ir>nresidues-ioffset; ir--) {
						offset_index_1 += alignment_atoms_in_res(firstpose, ir);
					}
				}

				for(core::Size ir=1, index=1; ir<=nresidues; ir++) { //For each offset, loop through all residues
					//The number of alignment atoms stored for this residue:
					core::Size atoms_in_this_res = alignment_atoms_in_res(firstpose,ir);

					//Offset the alignment vectors:
					for(core::Size ia=1; ia<=atoms_in_this_res; ia++) {
						avect2(1, offset_index_1) = alignmentvect2(1, index);
						avect2(2, offset_index_1) = alignmentvect2(2, index);
						avect2(3, offset_index_1) = alignmentvect2(3, index);
						index++; //Move to next atom.
						offset_index_1++; //Move this with index.
						//printf("%.2f, %.2f, %.2f\n", avect2(1,ia), avect2(2, ia), avect2(3, ia)); //DELETE ME
					}
					if(ir==(nresidues-ioffset)) offset_index_1=1; //Reset this if we've reached the wrap-around point.
				} //Looping through all residues

				core::Real curdist = calc_dist(vect1, vect2, clustmode, alignmentvect1, avect2, nresidues, firstpose);

				//printf("offset %i\t dist %.2f\n", ioffset, curdist); fflush(stdout); //DELETE ME
				if (ioffset==0 || curdist<dist ) { //If this offset has yielded the lowest RMS encountered so far.
					dist = curdist;
					offset = ioffset;
				}
			}
		} else if (clustmode==2) { //Clustering by backbone dihedrals
			//Make copies of vect1 and vect2
			utility::vector1 <core::Real> v2 = vect2;
			for(core::Size ioffset=0; ioffset<nresidues; ioffset+=option[v_cyclic_permutation_offset]()) {
				core::Size offset_index_1=2;
				if(ioffset>0) {
					for(core::Size ir=nresidues; ir>nresidues-ioffset; ir--) {
						offset_index_1 += alignment_torsions_in_res(firstpose, ir);
						if(ir==nresidues) offset_index_1+=3; //Extra 3 torsions stored at end of list.
					}
				}

				for(core::Size ir=1, index=1; ir<=nresidues; ir++) { //Loop through all residues
					core::Size torsions_in_this_res = alignment_torsions_in_res(firstpose, ir);
					if(ir==nresidues) torsions_in_this_res+=3; //Three additional torsions are stored for a backbone-cyclized peptide.
					for(core::Size itors=1; itors<=torsions_in_this_res; itors++) {
						v2[offset_index_1] = vect2[index];
						offset_index_1++;
						index++;
						if(ir==(nresidues-ioffset+1) && itors==1) offset_index_1=1; //Reset this at the wrap-around point.
					}
				}

				core::Real curdist = calc_dist(vect1, v2, clustmode, alignmentvect1, alignmentvect2, nresidues, firstpose);
				//printf("Offset %i\t Dist %.2f\n", ioffset, curdist); fflush(stdout); //DELETE ME

				if (ioffset==0 || curdist<dist ) { //If this offset has yielded the lowest distance encountered so far.
					dist = curdist;
					offset = ioffset;
				}
			}
		}
	} //if(option[v_cluster_cyclic_permutations]())
	else if (option[v_homooligomer_swap]() && firstpose.num_jump() > 0) { //consider all permutations of chains
		core::Size numchains = firstpose.num_jump()+1;
		for(core::Size i=1; i<=factorial(numchains); i++) {
			if(clustmode==1) {
				FArray2D < core::Real > alignmentvect2_swapped;
				swapalignmentvector(alignmentvect2, alignmentvect2_swapped, firstpose, i);
				core::Real curdist = calc_dist(vect1, vect2, clustmode, alignmentvect1, alignmentvect2_swapped, nresidues, firstpose);
				if(i==1 || curdist < dist) {
					dist=curdist;
					permutation=i;
				}				
			} else if(clustmode==2) {
				utility::vector1 < core::Real > vect2_swapped;
				swapvector(vect2, vect2_swapped, firstpose, i); //TODO -- CURRENTLY DOESN'T WORK!
				core::Real curdist = calc_dist(vect1, vect2_swapped, clustmode, alignmentvect1, alignmentvect2, nresidues, firstpose);
				if(i==1 || curdist < dist) {
					dist=curdist;
					permutation=i;
				}
			}
		}
		
	} //else if (option[v_homooligomer_swap]() && firstpose.num_jump() > 0) {
	else {
		dist = calc_dist(vect1, vect2, clustmode, alignmentvect1, alignmentvect2, nresidues, firstpose);
	}

	return dist;
}

//Function to circularly permutate a vector of backbone dihedral angles:
void slidearound(
	utility::vector1 <core::Real> &dihedrals,
	const core::Size offset,
	const core::Size rescount,
	const core::pose::Pose &pose //for reference only
) {
	utility::vector1 <core::Real> result;

	core::Size bb_dihed_count = 0;
	for(core::Size ir=1; ir<=pose.n_residue(); ir++) bb_dihed_count+=alignment_torsions_in_res(pose, ir);
	bb_dihed_count+=3; //For a true cyclic peptide, there are three additional vectors stored.

	result.resize(bb_dihed_count, 0.0);
	core::Size offset_index = 2;
	if(offset>0) {
		for(core::Size ir=rescount; ir>rescount-offset; ir--) offset_index+=alignment_torsions_in_res(pose, ir);			
		offset_index+=3; //3 additional at end
	}

	for(core::Size ir=1, index=1; ir<=rescount; ir++) {
		for(core::Size itors=1, ntors=alignment_torsions_in_res(pose, ir); itors<ntors; itors++) {
			result[offset_index]=dihedrals[index];
			offset_index++;
			index++;
			if(itors==1 && ir==(rescount-offset+1)) offset_index = 1; //Reset at the wrap-around point.
		}
	}
	dihedrals=result;
	return;
}

//Function to trim a vector of dihedral angles, and to add tranlation and Euler angle vectors for all the jumps in a pose to it.
void trim_and_add_jump_data (
	utility::vector1 < core::Real > &myvect,
	const core::Size numdihedrals,
	const core::pose::Pose &mypose
) {
	myvect.resize(numdihedrals + 6*mypose.num_jump()); //Discard extra information in myvect
	if(mypose.num_jump() > 0) { //If there are jumps, store translation and rotation vectors for each jump in mypose.
		for(core::Size i=numdihedrals+1, counter=1; i<=myvect.size(); counter++) {
			numeric::xyzVector < core::Real > transvect = mypose.jump(counter).get_translation();
			numeric::EulerAngles < core::Real > euler_angles(mypose.jump(counter).get_rotation());
			for(core::Size j=0; j<3; j++) myvect[i++] = transvect[j]; //Note that xyzVectors are zero-based (unlike everything else in Rosetta -- yay for consistency).
			myvect[i++] = euler_angles.phi_degrees();
			myvect[i++] = euler_angles.theta_degrees();
			myvect[i++] = euler_angles.psi_degrees();
			//The following is just for debugging:
			//printf("T=(%.2f,%.2f,%.2f)\nR=(%.2f,%.2f,%.2f)\n", transvect[0], transvect[1], transvect[2], euler_angles.phi_degrees(), euler_angles.theta_degrees(), euler_angles.psi_degrees());
		}
	}

	return;
}

//Fuction to shift the current center of the cluster, weighted by energies.
//This also performs PCA analysis, returning a matrix of PCA vectors.
//The PCA vectors correspond to backbone dihedral angles, with jumps (tx, ty, tz, rx, ry, rz) appended at the end.
//Angles are in degrees, transforms are in angstroms.
void shift_center_and_PCA(
	utility::vector1 <core::Real> &clustcenter,
	utility::vector1 < utility::vector1 < core::Real > > &pca_vector_list,
	utility::vector1 < core::Real > &coeff_list,
	const core::Size clustcenterindex,
	const utility::vector1 < utility::vector1 <core::Real> > &posedata,
	const utility::vector1 <core::Real> &poseenergies,
	const utility::vector1 <core::Size> &statesincluster, //This is the list of states assigned to this cluster, sorted by energies!
	const core::Size currentclusterindex,
	const core::Size clustmode,
	const core::pose::Pose &firstpose,
	const utility::vector1 <core::Size> &cluster_offsets, //Only used if clustering cyclic permutations together
	const utility::vector1 <core::Size> &cluster_oligomer_permutations, //Only used if clustering permutations of homooligomer chains
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	printf("\tShifting center of cluster %lu.\n", currentclusterindex);

	core::Size ndihedrals = count_bb_dihedrals(firstpose); //The number of dihedral angles to use in PCA analysis

	printf("\t\tStarting with structure %lu.\n", clustcenterindex); fflush(stdout);

	//Clear the storage containers:
	pca_vector_list.clear();
	coeff_list.clear();

	core::Real weighting_accumulator = option[v_weightbyenergy]() ? exp(-poseenergies[clustcenterindex]/option[v_kbt]()) : 1.0; //The weighting accumulator (the denominator) starts off with the partition weight of the cluster center.
	utility::vector1 <core::Real> coeff;
	coeff.push_back(weighting_accumulator); //coeff[1] is the coefficient for the first state

	utility::vector1 <core::Real> newclustcenter; //Used for clustmode=1 clustering
	utility::vector1 <core::Real> deltaclustcenter; //The shift in the cluster center

	FArray2D <core::Real> dummyarray; //Need this to satisfy requirements of storeposedata function.

	core::pose::Pose clustcenterpose=firstpose;

	alglib::real_2d_array Dphipsiomega; //A matrix whose columns are the phi, psi, and omega values of each state in the cluster, with the cluster centre subtracted off.  Columns will then be multiplied by the weighting coefficients.
	Dphipsiomega.setlength(statesincluster.size(), ndihedrals+6*clustcenterpose.num_jump()); //This should have a column for each state, and a row for each backbone dihedral angle and jump parameter.

	if(clustmode==1) {
		pose_from_posedata(firstpose, clustcenterpose, clustmode, clustcenter);
		storeposedata(clustcenterpose, newclustcenter, dummyarray, 2, extra_atom_list); //If clustermode is 1, newclustcenter holds the vector of backbone (and side chain) dihedrals at this point.  NOTE THAT THERE'S EXTRA INFORMATION AT THE END OF NEWCLUSTCENTER.
	} else if(clustmode==2) newclustcenter=clustcenter;

	trim_and_add_jump_data (newclustcenter, ndihedrals, clustcenterpose);  //Discard extra data and add jump data.

	utility::vector1 <core::Size> not_angle_list; //A list of entries in the newclustcenter vector that are NOT angles (translations instead of rotations).
	if(clustcenterpose.num_jump()>0) {
		for(core::Size i=ndihedrals+1, counter=1, njump=clustcenterpose.num_jump(); counter<=njump; counter++) {
			not_angle_list.push_back(i);
			not_angle_list.push_back(i+1);
			not_angle_list.push_back(i+2);
			i+=6;
		}
	}

	deltaclustcenter.resize(ndihedrals + 6*clustcenterpose.num_jump());
	for(core::Size i=1; i<=deltaclustcenter.size(); i++) deltaclustcenter[i]=0.0; //Initialize the array
	for(core::Size i=1; i<=newclustcenter.size(); i++) Dphipsiomega[0][i-1] = newclustcenter[i]; //Initialize the first column of Dphipsiomega

	if(statesincluster.size() > 1) { //If there is more than one state in this cluster.
		for (core::Size i=2; i<=statesincluster.size(); i++) { //Loop through all other states in the cluster.

			printf("\t\tCalculating influence of structure %lu.\n", statesincluster[i]); fflush(stdout);

			coeff.push_back(option[v_weightbyenergy]() ? exp(-poseenergies[statesincluster[i]]/option[v_kbt]()) : 1.0);
			weighting_accumulator+=coeff[i]; //The weighting accumulator will ultimately be the demoninator in the partition function.

			if(clustmode==1) {
				utility::vector1 <core::Real> currentdata;
				core::pose::Pose currentpose;
				pose_from_posedata(firstpose, currentpose, clustmode, posedata[statesincluster[i]]);
				if(option[v_homooligomer_swap]()) swapchains( currentpose, cluster_oligomer_permutations[statesincluster[i]] ); //Swap chains around.
				storeposedata(currentpose, currentdata, dummyarray, 2, extra_atom_list); //currentdata is a vector of backbone dihedrals at this point.  NOTE EXTRA DATA AT END -- SIDE-CHAIN DIHEDRALS.

				//If we're clustering cyclic permutations, it's necessary to offset the backbone dihedral vector:
				if(option[v_cluster_cyclic_permutations]()) slidearound(currentdata, cluster_offsets[statesincluster[i]], firstpose.n_residue(), firstpose);
				trim_and_add_jump_data (currentdata, ndihedrals, currentpose); //Discard extra data and add jump data.

				//Add this structure to Dphipsiomega:
				for(core::Size j=1; j<=currentdata.size(); j++) Dphipsiomega[i-1][j-1]=currentdata[j];

				//Add this to deltaclustcenter:
				for(core::Size j=1, diffmode=2; j<=deltaclustcenter.size(); j++) {
					diffmode = 2;
					if(j>ndihedrals && is_in_list(j, not_angle_list)) diffmode=1; //If we're past the dihedral list, check whether this entry is a translation or a rotation.  If it's a translation, diffmode=1.
					deltaclustcenter[j] += (coeff[i]*elementdiff(currentdata[j], newclustcenter[j], diffmode));
				}
			} else if (clustmode==2) {
				utility::vector1 <core::Real> currentdata = posedata[statesincluster[i]];
				//If we're clustering cyclic permutations, it's necessary to offset the backbone dihedral vector:
				if(option[v_cluster_cyclic_permutations]()) slidearound(currentdata, cluster_offsets[statesincluster[i]], firstpose.n_residue(), firstpose);
				for(core::Size j=1; j<=deltaclustcenter.size(); j++) deltaclustcenter[j] += (coeff[i]*elementdiff(currentdata[j], newclustcenter[j], clustmode));
				//Add this structure to Dphipsiomega:
				for(core::Size j=1; j<=ndihedrals; j++) Dphipsiomega[i-1][j-1]=currentdata[j];
			}
		} //Looping through states in the cluster

		for(core::Size i=1; i<=deltaclustcenter.size(); i++) newclustcenter[i]+=(deltaclustcenter[i]/weighting_accumulator); //Shift the cluster center.

		check_angle_bounds(newclustcenter, not_angle_list);

		//Subtract clustcenter from Dphipsiomega, and multiply each resultant column by coeff:
		for(core::Size i=1, imax=statesincluster.size(); i<=imax; i++) { //Loop through all states
			for(core::Size j=1, jmax=newclustcenter.size(); j<=jmax; j++) { //Loop through all backbone dihedrals
				Dphipsiomega[i-1][j-1] = coeff[i]/weighting_accumulator*elementdiff(Dphipsiomega[i-1][j-1],newclustcenter[j], (is_in_list(j, not_angle_list)?1:2));
				//Dphipsiomega[i-1][j-1] = coeff[i]*(Dphipsiomega[i-1][j-1]-newclustcenter[j]);
				//printf("%.4f\t", Dphipsiomega[i-1][j-1]); //DELETE ME
			}
			//printf("\n"); fflush(stdout); //DELETE ME
		}

		//PCA analysis:
		printf("\tPerforming principal component analysis for cluster %lu.\n", currentclusterindex); fflush(stdout);
		alglib::ae_int_t PCAresult = 0;
		alglib::real_1d_array variances;
		variances.setlength(ndihedrals+6*clustcenterpose.num_jump());
		alglib::real_2d_array PCmatrix;
		PCmatrix.setlength(ndihedrals+6*clustcenterpose.num_jump(),ndihedrals+6*clustcenterpose.num_jump());
		alglib::pcabuildbasis(Dphipsiomega, statesincluster.size(), ndihedrals+6*clustcenterpose.num_jump(), PCAresult, variances, PCmatrix);

		//Store the variances and principal component vectors for output from this function.  Only the first N-1 should be nonzero, where N is the number of states in the cluster.
		for(core::Size i=1; i<std::min(statesincluster.size(), ndihedrals+6*clustcenterpose.num_jump()); i++) {
			coeff_list.push_back(variances[i-1]);
			//printf("%.8f\n", variances[i-1]);fflush(stdout); //DELETE ME
			utility::vector1 < core::Real > PCA_vector;
			for(core::Size j=1; j<=ndihedrals+6*clustcenterpose.num_jump(); j++) PCA_vector.push_back(PCmatrix[j-1][i-1]);
			pca_vector_list.push_back(PCA_vector);
		}

		if(clustmode==1) {
			pose_from_posedata(firstpose, clustcenterpose, 2, newclustcenter, false);
			storeposedata(clustcenterpose, newclustcenter, dummyarray, clustmode, extra_atom_list);
		}

		clustcenter=newclustcenter;

	}

	return;
}

//Function to write out a PCA file:
//Note that a PCA file is a list corresponding to backbone dihedrals, followed by entries for each jump.  Each jump has six entries: translation(x,y,z) and rotation(phi, theta, psi).
void output_PCA (
	const utility::vector1 < utility::vector1 < core::Real > > &pca_vector_list,
	const utility::vector1 < core::Real > &coeff_list,
	const core::Size &clustnumber
) {
	FILE* outputfile;
	char filename [128];
	sprintf(filename, "PCA_%lu.txt", clustnumber);
	outputfile = fopen(filename, "w");
	//printf("\tOpened %s for write.\n", filename); fflush(stdout); //DELETE ME

	if(coeff_list.size() > 0) {
		fprintf (outputfile, "%lu\t%lu\n", pca_vector_list[1].size(), coeff_list.size());
		for(core::Size i=1; i<=coeff_list.size(); i++) fprintf(outputfile, "%.14e\t", coeff_list[i]);
		fprintf(outputfile, "\n");
		for(core::Size i=1; i<=pca_vector_list.size(); i++) {
			for(core::Size j=1; j<=pca_vector_list[i].size(); j++) fprintf(outputfile, "%.14e\t", pca_vector_list[i][j]);
			fprintf(outputfile, "\n");
		}
	} else {
		fprintf (outputfile, "0\t0\n");
	}
		
	fclose(outputfile);
	printf("\tWrote %s.\n", filename); fflush(stdout);
	return;
}

//Function to sort the list of states in a cluster by energies
void sortclusterlist(
	utility::vector1 <core::Size> &statelist,
	utility::vector1 <core::Real> &poseenergies
) {
	//This won't be the world's most efficient sort.  I'll use a selection sort algorithm:
	core::Real lowestE = 0.0;
	core::Size lowestEentry = 0;
	core::Size buffer = 0;

	for(core::Size i=1; i<statelist.size(); i++) {
		for(core::Size j=i; j<=statelist.size(); j++) {
			if(j==i || poseenergies[statelist[j]]<lowestE) {
				lowestE = poseenergies[statelist[j]];
				lowestEentry=j;
			}
		}
		if(lowestEentry!=i) {
			buffer=statelist[i];
			statelist[i]=statelist[lowestEentry];
			statelist[lowestEentry]=buffer;
		}
	} 

	return;
}

//Function to add cyclic constraints to a pose:
void addcyclicconstraints (core::pose::Pose &mypose) {
	using namespace scoring::constraints;
	using namespace scoring::func;
	using namespace core::id;

	core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(mypose, mypose.n_residue());

	const std::string hstring = (mypose.residue(1).has("H") ? "H" : "CD");

	if(mypose.residue(1).has("H")) { //If this is not a proline or D-proline
		core::Size hindex = mypose.residue(1).atom_index("H");

		//Rebuild the N-terminal proton.  This has to be done in a slightly irritating way because Rosetta doesn't really like the fact
		//that the last residue is connected to the first:
		{
			core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
			core::pose::Pose dialanine;
			core::pose::make_pose_from_sequence(dialanine, "AA", *standard_residues, true); //The termini are OPEN.
			core::Real omegaval = numeric::dihedral_degrees(
					mypose.residue(mypose.n_residue()).xyz("CA"),
					mypose.residue(mypose.n_residue()).xyz("C"),
					mypose.residue(1).xyz("N"),
					mypose.residue(1).xyz("CA")
				);
			dialanine.set_omega(1, omegaval);
	
			core::id::AtomID_Map< core::id::AtomID > amap;
			core::pose::initialize_atomid_map(amap,dialanine,core::id::BOGUS_ATOM_ID);
			amap[AtomID(dialanine.residue(1).atom_index("CA"),1)] = AtomID(mypose.residue(mypose.n_residue()).atom_index("CA"),mypose.n_residue());
			amap[AtomID(dialanine.residue(1).atom_index("C"),1)] = AtomID(mypose.residue(mypose.n_residue()).atom_index("C"),mypose.n_residue());
			amap[AtomID(dialanine.residue(1).atom_index("O"),1)] = AtomID(mypose.residue(mypose.n_residue()).atom_index("O"),mypose.n_residue());
			amap[AtomID(dialanine.residue(2).atom_index("N"),2)] = AtomID(mypose.residue(1).atom_index("N"),1);
			amap[AtomID(dialanine.residue(2).atom_index("CA"),2)] = AtomID(mypose.residue(1).atom_index("CA"),1);
			core::scoring::superimpose_pose( dialanine, mypose, amap );

			mypose.conformation().set_xyz(AtomID(hindex, 1), dialanine.residue(2).xyz("H"));
		}
	}

	mypose.conformation().declare_chemical_bond(1, "N", mypose.n_residue(), "C"); //Declare a chemical bond between the N and C termini.

	{//Peptide bond length constraint:
		FuncOP harmfunc1 = new HarmonicFunc( 1.3288, 0.02);
		ConstraintCOP distconst1 = new AtomPairConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ) ,
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				harmfunc1
			);
		mypose.add_constraint (distconst1);
	}

	{	//Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1 = new CircularHarmonicFunc( PI, 0.02);
		ConstraintCOP dihedconst1 = new DihedralConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("O") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index(hstring) , 1) ,
				circharmfunc1
			);
		ConstraintCOP dihedconst2 = new DihedralConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ) , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("CA") , 1) ,
				circharmfunc1
			);

		mypose.add_constraint (dihedconst1);
		mypose.add_constraint (dihedconst2);
	}

	{	//Peptide bond angle constraints:
		FuncOP circharmfunc2a = new CircularHarmonicFunc( CNCa_ANGLE/180.0*PI, 0.02);
		FuncOP circharmfunc2b = new CircularHarmonicFunc( CNH_ANGLE/180.0*PI, 0.02);
		FuncOP circharmfunc2c = new CircularHarmonicFunc( CaCN_ANGLE/180.0*PI, 0.02);
		FuncOP circharmfunc2d = new CircularHarmonicFunc( OCN_ANGLE/180.0*PI, 0.02);

		ConstraintCOP angleconst1 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("CA") , 1) ,
				circharmfunc2a
			);
		ConstraintCOP angleconst2 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index(hstring) , 1) ,
				circharmfunc2b
			);
		ConstraintCOP angleconst3 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ) , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				circharmfunc2c
			);
		ConstraintCOP angleconst4 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("O") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				circharmfunc2d
			);

		mypose.add_constraint (angleconst1);
		mypose.add_constraint (angleconst2);
		mypose.add_constraint (angleconst3);
		mypose.add_constraint (angleconst4);
	}

	return;
}

//Function to align one pose to another with an offset in the residue count.
void align_with_offset (
	core::pose::Pose &pose1,
	const core::pose::Pose &pose2, //the target pose -- doesn't change.
	const core::Size offset,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	using namespace core::id;

	AtomID_Map< AtomID > amap;
	signed int ir;
	core::pose::initialize_atomid_map(amap,pose1, BOGUS_ATOM_ID);
	for(core::Size ir2 = 1; ir2 <= (core::Size)pose2.n_residue(); ++ir2) {
		ir= ir2-offset;
		if(ir<=0) ir+=pose1.n_residue(); 
		for(core::Size ia = 1; ia <= (core::Size)pose2.residue(ir2).nheavyatoms(); ++ia) {
			if(use_in_rmsd_offset(pose1,pose2,ir2,ia, offset, extra_atom_list)) {
				amap[AtomID(ia,ir)] = AtomID(ia,ir2);
			}
		}
	}
	core::scoring::superimpose_pose( pose1, pose2, amap );

	return;
}

//Function to mutate a pose to a chain of alanines -- EXCEPT CYSTEINES:
void mutate_to_ala(core::pose::Pose &mypose)
{
	for(core::Size ir=1; ir<=mypose.n_residue(); ir++) { //Loop through all residues
		if(!mypose.residue(ir).type().is_beta_aa() && !mypose.residue(ir).type().is_alpha_aa()) continue; //Skip non-amino acid residues
		if(mypose.residue(ir).name1()=='C') continue; //Skip cysteines
		std::string aaname = "ALA";
		if( core::chemical::is_canonical_D_aa( mypose.residue(ir).aa() ) ) aaname = "DALA"; //If it's a D-amino acid, mutate to D-alanine
		else if ( is_beta_aminoacid(mypose, ir) ) aaname = "B3A"; //If it's a beta-amino acid, mutate to beta-3-alanine.
		protocols::simple_moves::MutateResidue mutres(ir, aaname); //Mutate residue mover
		mutres.apply(mypose); //Apply the mutation
	}

	mypose.update_residue_neighbors();

	return;
}

//Function to check whether a pose contains at least one beta-amino acid residue:
bool contains_beta (const core::pose::Pose &mypose) {
	for(core::Size ir=1, nres=mypose.n_residue(); ir<=nres; ir++)
		if(is_beta_aminoacid(mypose, ir)) return true;
	return false;
}

//Function to check whether two poses have beta amino acid residues at matching positions:
bool betas_match (
	const core::pose::Pose &pose1,
	const core::pose::Pose &pose2
) {
	for(core::Size ir=1,nres=pose1.n_residue(); ir<=nres; ir++) 
		if( is_beta_aminoacid(pose1,ir) != is_beta_aminoacid(pose2,ir) ) return false;
	return true;
}

//Function to add user-specified constraints (specified with a CST file) to a pose:
void add_user_constraints (
	core::pose::Pose &mypose
) {
	using namespace protocols::simple_moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(!option[v_cst_file].user()) return; //Do nothing if no constraint file is specified.

	ConstraintSetMoverOP cst_maker = new ConstraintSetMover();
	std::string cstfile = option[v_cst_file]();

	cst_maker->constraint_file(cstfile);
	cst_maker->add_constraints(true); //Add constraints to anything else already there.
	cst_maker->apply(mypose);	

	return;
}

//MAIN
int main( int argc, char * argv [] ) {

	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	using namespace protocols::cluster;
	using namespace basic::options::OptionKeys::cluster;
	using namespace std;

	printf("Starting bettercluster.cc\nFile created 6 May 2013 by Vikram K. Mulligan\n"); fflush(stdout);

	register_options();
	devel::init(argc, argv);
	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::getScoreFunction();
	if(option[v_cst_file].user()) { //If a constraints file has been specified by the user, turn on the atom_pair, angle, and dihedral constraint weights unless otherwise on.
		if(sfxn->get_weight(atom_pair_constraint) < 1.0e-6) sfxn->set_weight(atom_pair_constraint, 1.0);
		if(sfxn->get_weight(angle_constraint) < 1.0e-6) sfxn->set_weight(angle_constraint, 1.0);
		if(sfxn->get_weight(dihedral_constraint) < 1.0e-6) sfxn->set_weight(dihedral_constraint, 1.0);
	}

	//Parse the clustering mode:
	string clusterby = option[v_clusterby]();
	core::Size clustermode =0;
	if (clusterby=="bb_cartesian") {
		clustermode = 1;
		printf("Clustering by Cartesian coordinates of backbone atoms.\n");
	} else if (clusterby=="bb_dihedral") {
		clustermode = 2;
		printf("Clustering by backbone dihedral angles.\n");
	} else {
		printf("Error!  \"%s\" is an invalid input for the -v_clusterby flag.  Use the --help flag to list accepted inputs.  Crashing.\n", clusterby.c_str());
		fflush(stdout); exit(1);
	}

	//Check whether the v_mutate_to_ala flag has been set:
	if(option[v_mutate_to_ala]())
		printf("The -v_mutate_to_ala flag was set.  Input structures will be mutated to a chain of (alpha-D-, alpha-L-, or beta-3-) alanines (with the exception of cysteine residues).\n");

	//Check that the user hasn't specified the same residue multiple times in -v_ignoreresidue:
	if(option[v_ignoreresidue]().size()>1) {
		for(core::Size i=2; i<=option[v_ignoreresidue]().size(); i++) {
			for(core::Size j=1; j<i; j++) {
				if(option[v_ignoreresidue]()[i] == option[v_ignoreresidue]()[j]) {
					printf("Error!  The same residue must not be specified multiple times with the -v_ignoreresidue flag.  Crashing.\n");
					fflush(stdout); exit(1);
				}
			}
		}
	}

	//Check whether the v_homooligomer_swap flag has been set:
	if(option[v_homooligomer_swap]()) {
		if(clusterby!="bb_cartesian") {
			printf("Error!  Multi-chain structures can currently only be scored with the \"-v_clusterby bb_cartesian\" flag!  Crashing.\n");
			fflush(stdout); exit(1);
		}
		if(option[v_cyclic]()) {
			printf("Error!  The -v_homooligomer_swap flag is not currently compatible with the -v_cyclic flag!  Crashing.\n");
			fflush(stdout); exit(1);
		}
		printf("The -v_homooligomer_swap flag was specified.  When calculating RMSD values, all permutations of chains will be considered in the alignment.\n");
	}

	//Parse the clustering radius:
	const core::Real R_cluster = option[v_clusterradius]();
	if (R_cluster < 1.0e-12) {
		printf("Error!  The clustering radius must be greater than zero.  Crashing.\n");
		fflush(stdout); exit(1);
	} else {
		string unitsstring;
		if(clustermode==1) unitsstring = "Angstroms";
		else if (clustermode==2) unitsstring = "degrees";
		printf("Using a cluster radius of %.2f %s.\n", R_cluster, unitsstring.c_str());
	}

	//Parse v_kbt:
	const core::Real kbt = option[v_kbt]();
	if(kbt<1.0e-12) {
		printf("Error!  The value of k_B*T specified with the -v_kbt flag must be greater than zero.  Crashing.\n");
		fflush(stdout); exit(1);
	} else {
		if(option[v_weightbyenergy]()) {
			printf("Weighting structures by energy when calculating cluster centers.  Setting k_B*T=%.2f\n", kbt);
		}
	}

	//Alter the scoring function for v_cyclic:
	if(option[v_cyclic]()) {
		printf("Setting constraint weights for a peptide bond between the N- and C-termini (-v_cyclic flag).\n");
		sfxn->set_weight(atom_pair_constraint, 1.0);
		sfxn->set_weight(dihedral_constraint, 1.0);
		sfxn->set_weight(angle_constraint, 1.0);
	}

	//Checks releated to v_cluster_cyclic_permutations:
	if(option[v_cluster_cyclic_permutations]()) {
		if(!option[v_cyclic]()) {
			printf("Error!  The -v_cluster_cyclic_permutations flag cannot be used without the -v_cyclic flag.  Crashing gracelessly.\n");
			fflush(stdout); exit(1);
		}
		printf("Cyclic permutations will be tried when aligning structures during clustering (-v_cluster_cyclic_permutations flag).\n");
	}


	//Parse the user-specified list of additional atoms to use in the RMSD calculation:
	utility::vector1<core::id::NamedAtomID> extra_atom_list;
	parse_extra_atom_list(extra_atom_list); //Does nothing if no list provided.

	fflush(stdout);
	
	protocols::relax::FastRelax frlx(sfxn, option[v_relaxrounds]());

	//printf("ping1\n"); fflush(stdout); //DELETE ME
	//core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	//printf("ping2\n"); fflush(stdout); //DELETE ME
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	core::Size count = 0;
	core::Real lowestE = 0; //Lowest energy encountered so far.
	core::Size lowestE_index = 0; //The number of the pose with the lowest energy encountered.
	utility::vector1 <core::Real> poseenergies; //Vector of energies of the poses.
	utility::vector1 < utility::vector1 <core::Real> > posedata; //Vector of vectors to store the data that will be used for clustering.
	utility::vector1 < FArray2D < core::Real > > alignmentdata; //Vector of array of x,y,z coordinates of atoms to be used for alignment in Cartesian clustering.
	utility::vector1 < core::Size > cluster_assignments; //List of which cluster each structure is assigned to.
	utility::vector1 < core::Size > cluster_offsets; //List of offsets for cyclic permutations when clustering (measured in number of amino acid positions we've offset by).
	utility::vector1 < core::Size > cluster_oligomer_permutations; //List of permutations of oligomers if the user has specified the option to do that.

	printf("Scoring energies of all input structures.\n"); fflush(stdout);
	core::pose::Pose firstpose;
	bool firstpose_hasbeta = false; //Are there beta-amino acids in the first input pose?

	while( input.has_another_pose() ) { //For every input structure
		count++;

		if (count % 100 == 0) {printf("."); fflush(stdout); }
		cluster_assignments.push_back(0); //Initially, every structure is assigned to cluster 0 (unassigned).
		if(option[v_cluster_cyclic_permutations]()) cluster_offsets.push_back(0); //Initially, we assume that each structure may be aligned without any cyclic permutation, so these are all zero.
		if(option[v_homooligomer_swap]()) cluster_oligomer_permutations.push_back(1); //Initially, we assume the first permutation of homooligomer subunits for each structure.

		core::pose::Pose pose; //Create the pose
		input.fill_pose( pose ); //Import it

		if(pose.num_jump() > 0 && clusterby == "bb_dihedral") {
			printf("Error!  Backbone dihedral clustering is not currently compatible with multi-chain PDB files.  Crashing.\n"); //TODO -- make this compatible!
			fflush(stdout); exit(1);
		}

		if(pose.num_jump() > 0 && option[v_homooligomer_swap]()) {
			const std::string chain1seq = pose.chain_sequence(1);
			for(core::Size i=2, imax=pose.num_jump()+1; i<=imax; i++) {
				if (pose.chain_sequence(i) != chain1seq) {
					printf("Error!  When using the v_homooligomer_swap option, the lengths and sequences of all chains in the input structures must be identical.  Crashing.\n");
					fflush(stdout); exit(1);
				}
			}
		}

		if(option[v_cyclic]()) {
			addcyclicconstraints(pose); //If this is a cyclic peptide, add constraints for the terminal peptide bond:
			if(option[v_cluster_cyclic_permutations]() && clustermode==1 && option[v_CB]()) {
				for(core::Size ir=1; ir<=pose.n_residue(); ir++) {
					if(!pose.residue(ir).has("CB")) { //Error if we don't have the same set of atoms on which to cluster in each residue.
						printf("Error!  When using Cartesian clustering and clustering cyclic permutations, all residues in the input structures must have the same number of atoms to be used in the alignment.  If beta carbons are included, the sequence cannot contain glycine.\nCrashing.\n");
						fflush(stdout); exit(1);
					}
				}
			}
		}

		add_user_constraints(pose); //This checks for a user-specified CST file, and does nothing if there isn't one.  (Adds constraints if there is one.)

		if(option[v_prerelax]()) {
			make_disulfides(pose); //Add user-specified disulfide bonds.
			frlx.apply(pose); //Relax the pose if the user has specified that it be relaxed.
		}

		if(option[v_mutate_to_ala]()) mutate_to_ala(pose); //Mutate the pose to a chain of alanines if necesssary.

		make_disulfides(pose); //Add user-specified disulfide bonds.

		(*sfxn)(pose); //Score the input pose
		if(count==1) {
			firstpose=pose; //Store the first pose.
			firstpose_hasbeta = contains_beta(firstpose); //Check whether this pose contains beta-amino acids.
		}
		else if (firstpose_hasbeta && !betas_match(firstpose, pose)) { //If this is not the first pose AND the first pose had beta-amino acids, check that this pose has beta-amino acids at the same positions.
			printf("Error!  If input structures contain beta-amino acids, they must have the same number of beta-amino acid residues, and at the same positions.\nCrashing!\n");
			fflush(stdout);
			exit(1);
		}
		poseenergies.push_back(pose.energies().total_energy()); //Store the pose energy

		//Store the pose data that will be used for clustering:
		storeposedata(pose, posedata, alignmentdata, clustermode, extra_atom_list);

		if(count==1 || pose.energies().total_energy() < lowestE) {
			lowestE=pose.energies().total_energy();
			lowestE_index = count;
		}
	}

	printf("\nClustering, starting with lowest-energy structure as the center of the first cluster.\n"); fflush(stdout);
	core::Size unclustered_count = count; //The number of structures that remain to be clustered
	core::Size cluster_count=0; //The number of clusters created
	utility::vector1 < utility::vector1 <core::Real> > clustcenter_list; //Vector of vectors to store the data for the cluster centres
	utility::vector1 < utility::vector1 <core::Size> > clusterlist_sortedbyenergy; //A vector storing lists of the states assigned to each cluster, sorted by energy.

	while(unclustered_count > 0) {
		//Make a new cluster
		cluster_count++;

		//The cumulative weighting of the points considered so far
		//core::Real weighting_accumulator = option[v_weightbyenergy]() ? exp(-lowestE/option[v_kbt]()) : 1.0;
		//core::Real currentweighting = 0.0;

		//Assign the lowest energy structure to the next cluster
		cluster_assignments[lowestE_index] = cluster_count;
		utility::vector1 <core::Size> cluster_sortedbyenergy;
		cluster_sortedbyenergy.push_back(lowestE_index);
		printf ("Started cluster %lu and added structure %lu to it.\n", cluster_count, lowestE_index);
		unclustered_count--;
		utility::vector1 <core::Real> clustcenter = posedata[lowestE_index]; //Set the center of the current cluster

		//Make a list of unassigned candidate structures:
		utility::vector1 <core::Size> candidatelist;
		for(core::Size istruct=1; istruct<=count; istruct++) {
			if(cluster_assignments[istruct]==0) candidatelist.push_back(istruct);
		}
		printf("\tMade list of %lu unassigned structures.\n", candidatelist.size()); fflush(stdout);

		//Assign members of the candidate list to the current cluster if they fall within R_cluster of the cluster center.
		core::Real currentdist=0.0;
		core::Size currentcyclicoffset=0; //Only used for calculating cyclic permutations
		core::Size current_oligomer_permutation=0; //Only used for calculating permutations when swapping around homodimers.
		for(core::Size istruct=1; istruct<=candidatelist.size(); istruct++) {

			currentdist=calc_dist(clustcenter, posedata[candidatelist[istruct]], clustermode, alignmentdata[lowestE_index], alignmentdata[candidatelist[istruct]], firstpose.n_residue(), firstpose, currentcyclicoffset, current_oligomer_permutation);	

			if (currentdist < R_cluster) {
				printf("\tAdding structure %lu (%.6f from the cluster center).\n", candidatelist[istruct], currentdist); fflush(stdout);
				cluster_assignments[candidatelist[istruct]]=cluster_count;
				cluster_sortedbyenergy.push_back(candidatelist[istruct]);
				if(option[v_cluster_cyclic_permutations]()) cluster_offsets[candidatelist[istruct]]=currentcyclicoffset;
				if(option[v_homooligomer_swap]()) cluster_oligomer_permutations[candidatelist[istruct]]=current_oligomer_permutation;
			}
		} //for loop through all candidates

		sortclusterlist(cluster_sortedbyenergy, poseenergies); //Sort the list of states in this cluster by energy
		clusterlist_sortedbyenergy.push_back(cluster_sortedbyenergy); //Add the list of states in this cluster to the list of states in each cluster
		utility::vector1 < utility::vector1 < core::Real > > pca_vector_list; //A place to store the PCA vectors
		utility::vector1 < core::Real > coeff_list; //A place to store coefficients (amplitudes) for each PCA vector
		shift_center_and_PCA(clustcenter, pca_vector_list, coeff_list, lowestE_index, posedata, poseenergies,
			cluster_sortedbyenergy, cluster_count, clustermode, firstpose, cluster_offsets,
			cluster_oligomer_permutations, extra_atom_list);
		output_PCA(pca_vector_list, coeff_list, cluster_count);
		clustcenter_list.push_back(clustcenter); //Store the current cluster center

		bool juststartedsearch=true;
		for(core::Size i=1; i<=poseenergies.size(); i++) {
			if(cluster_assignments[i]!=0) continue; //Continue if this structure has been assigned
			if(juststartedsearch || poseenergies[i]<lowestE) { //If this is the first unassigned encountered OR the lowest energy unassigned encountered so far
				juststartedsearch=false;
				lowestE=poseenergies[i];
				lowestE_index=i;
			}
		}

		if (juststartedsearch) break; //If this is still true, there were no unassigned structures.
	}

	//Outputs:
	printf("Cluster\tStructure\tFile_out\n");

	for (core::Size i=1; i<=clusterlist_sortedbyenergy.size(); i++) {
		core::pose::Pose pose1;
		for(core::Size j=1; j<=clusterlist_sortedbyenergy[i].size(); j++) {
			core::pose::Pose temppose;
			if(option[v_limit_structures_per_cluster]()==0 || j<=static_cast<core::Size>(option[v_limit_structures_per_cluster]())) {
				pose_from_posedata (firstpose, temppose, clustermode, posedata[ clusterlist_sortedbyenergy[i][j] ]);
				if(j>1 && option[v_homooligomer_swap]()) swapchains( temppose, cluster_oligomer_permutations[clusterlist_sortedbyenergy[i][j]] ); //Swap chains around.
				if(j==1) pose1=temppose;
				else align_with_offset(temppose, pose1, (option[v_cluster_cyclic_permutations]() ? cluster_offsets[ clusterlist_sortedbyenergy[i][j] ] : 0), extra_atom_list );
			}
			char outfile[128];
			if(option[v_silentoutput]()) {
				char curstructtag[128];
				sprintf(curstructtag, "c.%lu.%lu", i, j);
				io::silent::SilentFileData outsilentfiledata;
				io::silent::SilentStructOP outsilentstruct = io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
				outsilentstruct->fill_struct(temppose, curstructtag);
				if(j==1) {
					sprintf(outfile, "clusters_firstmember.out");
					outsilentfiledata.write_silent_struct((*outsilentstruct), outfile);
				}
				if(option[v_limit_structures_per_cluster]()==0 || j<=static_cast<core::Size>(option[v_limit_structures_per_cluster]())) {
					if(option[v_limit_structures_per_cluster]()==0) sprintf(outfile, "clusters_allmembers.out");
					else sprintf(outfile, "clusters_first_%lu_members.out", static_cast<core::Size>(option[v_limit_structures_per_cluster]()) );
					outsilentfiledata.write_silent_struct((*outsilentstruct), outfile);
				}
				sprintf(outfile, "c.%lu.%lu", i, j);											
			} else {
				sprintf(outfile, "c.%lu.%lu.pdb", i, j);
				if(option[v_limit_structures_per_cluster]()==0 || j<=static_cast<core::Size>(option[v_limit_structures_per_cluster]())) temppose.dump_pdb(outfile);
			}
			printf("%lu\t%lu\t%s\n", i, clusterlist_sortedbyenergy[i][j], outfile);
		}
	}

	/*printf("\nWriting cluster centers.\n"); fflush(stdout);
	for(core::Size i=1; i<=clustcenter_list.size(); i++) {
		core::pose::Pose cenpose;
		pose_from_posedata(firstpose, cenpose, clustermode, clustcenter_list[i]);
		char outfile[64];
		sprintf(outfile, "center_%lu.pdb", i);
		cenpose.dump_pdb(outfile);
	}*/

	printf("\n*****JOB COMPLETED*****\n"); fflush(stdout);
	return 0;
}

