// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/threefold_linker/TBMB_Helper.cc
/// @brief A derived class of the ThreefoldLinkerMoverHelper base class, used to set up
/// the 1,3,5-tris(bromomethyl)benzene (TBMB) cross-linker.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/threefold_linker/TBMB_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>

// Protocols headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateAngleConstraint.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraint.hh>
#include <protocols/simple_moves/ClearConstraintsMover.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.threefold_linker.TBMB_Helper" );

namespace protocols {
namespace cyclic_peptide {
namespace threefold_linker {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
TBMB_Helper::TBMB_Helper() //:
//TODO initialize data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
TBMB_Helper::TBMB_Helper( TBMB_Helper const &src ) :
	ThreefoldLinkerMoverHelper( src )
	//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
TBMB_Helper::~TBMB_Helper(){}

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Given a pose and a selection of exactly three residues, add the TBMB linker,
/// align it crudely to the selected residues, and set up covalent bonds.
void
TBMB_Helper::add_linker_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	core::Size cys1, cys2, cys3;
	ThreefoldLinkerMoverHelper::get_sidechain_indices( selection, cys1, cys2, cys3 );

	runtime_assert_string_msg( pose.residue_type(cys1).aa() == core::chemical::aa_cys || pose.residue_type(cys1).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::add_linker_asymmetric(): The first residue selected by the ResidueSelector is not an L- or D-cysteine.");
	runtime_assert_string_msg( pose.residue_type(cys2).aa() == core::chemical::aa_cys || pose.residue_type(cys2).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::add_linker_asymmetric(): The second residue selected by the ResidueSelector is not an L- or D-cysteine.");
	runtime_assert_string_msg( pose.residue_type(cys3).aa() == core::chemical::aa_cys || pose.residue_type(cys3).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::add_linker_asymmetric(): The third residue selected by the ResidueSelector is not an L- or D-cysteine.");

	//Mutate the cysteines to CYX:
	protocols::simple_moves::MutateResidue mut1( cys1, pose.residue_type(cys1).is_d_aa() ? "CYX" : "DCYX" );
	protocols::simple_moves::MutateResidue mut2( cys2, pose.residue_type(cys2).is_d_aa() ? "CYX" : "DCYX" );
	protocols::simple_moves::MutateResidue mut3( cys3, pose.residue_type(cys3).is_d_aa() ? "CYX" : "DCYX" );
	mut1.set_preserve_atom_coords(true);
	mut2.set_preserve_atom_coords(true);
	mut3.set_preserve_atom_coords(true);
	mut1.apply(pose);
	mut2.apply(pose);
	mut3.apply(pose);

	//Create the TBMB and put it into a pose of its own:
	core::chemical::ResidueTypeSetCOP standard_residues( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TBMB") ) );
	core::pose::Pose tbmb_pose;
	tbmb_pose.append_residue_by_jump(*new_rsd, 1);

	//Align the TBMB pose to the original pose:
	std::map <core::id::AtomID, core::id::AtomID> alignment_atoms;
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V1"), 1) ] = core::id::AtomID( pose.residue_type(cys1).atom_index("SG"), cys1);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V2"), 1) ] = core::id::AtomID( pose.residue_type(cys2).atom_index("SG"), cys2);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V3"), 1) ] = core::id::AtomID( pose.residue_type(cys3).atom_index("SG"), cys3);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM1"), 1) ] = core::id::AtomID( pose.residue_type(cys1).atom_index("V1"), cys1);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM2"), 1) ] = core::id::AtomID( pose.residue_type(cys2).atom_index("V1"), cys2);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM3"), 1) ] = core::id::AtomID( pose.residue_type(cys3).atom_index("V1"), cys3);
	core::scoring::superimpose_pose( tbmb_pose, pose, alignment_atoms );

	//Merge the poses:
	pose.append_residue_by_jump( tbmb_pose.residue(1), cys1 );
	core::Size const tbmb_res( pose.total_residue() );

	//Declare covalent bonds:
	add_linker_bonds_asymmetric( pose, cys1, cys2, cys3, tbmb_res );
}

/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
/// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
void
TBMB_Helper::add_linker_bonds_asymmetric(
	core::pose::Pose &pose,
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size const linker_index
) const {
	//Declare covalent bonds:
	protocols::cyclic_peptide::DeclareBond bond1, bond2, bond3;
	bond1.set( linker_index, "CM1", res1, "SG", false );
	bond2.set( linker_index, "CM2", res2, "SG", false );
	bond3.set( linker_index, "CM3", res3, "SG", false );
	bond1.apply(pose);
	bond2.apply(pose);
	bond3.apply(pose);
}


/// @brief Given a pose and a selection of exactly three residues, add the TBMB linker,
/// align it crudely to the selected residues, and set up covalent bonds.
/// @details Version for symmetric poses.
void
TBMB_Helper::add_linker_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	core::Size cys1, cys2, cys3;
	ThreefoldLinkerMoverHelper::get_sidechain_indices( selection, cys1, cys2, cys3 );

	runtime_assert_string_msg( pose.residue_type(cys1).aa() == core::chemical::aa_cys || pose.residue_type(cys1).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::add_linker_asymmetric(): The first residue selected by the ResidueSelector is not an L- or D-cysteine.");
	runtime_assert_string_msg( pose.residue_type(cys2).aa() == core::chemical::aa_cys || pose.residue_type(cys2).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::add_linker_asymmetric(): The second residue selected by the ResidueSelector is not an L- or D-cysteine.");
	runtime_assert_string_msg( pose.residue_type(cys3).aa() == core::chemical::aa_cys || pose.residue_type(cys3).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::add_linker_asymmetric(): The third residue selected by the ResidueSelector is not an L- or D-cysteine.");

	//Mutate the cysteines to CYX:
	protocols::simple_moves::MutateResidue mut1( cys1, pose.residue_type(cys1).is_d_aa() ? "CYX" : "DCYX" );
	mut1.set_preserve_atom_coords(true);
	mut1.apply(pose);

	//Create the TBMB and put it into a pose of its own:
	core::chemical::ResidueTypeSetCOP standard_residues( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TBMBSYMM") ) );
	core::pose::Pose tbmb_pose;
	tbmb_pose.append_residue_by_jump(*new_rsd, 1);

	//Align the TBMB pose to the original pose:
	std::map <core::id::AtomID, core::id::AtomID> alignment_atoms;
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V1"), 1) ] = core::id::AtomID( pose.residue_type(cys1).atom_index("SG"), cys1);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM1"), 1) ] = core::id::AtomID( pose.residue_type(cys1).atom_index("V1"), cys1);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V3"), 1) ] = core::id::AtomID( pose.residue_type(cys2).atom_index("V1"), cys2); //Imperfect alignment, but minimization will take care of it.
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V5"), 1) ] = core::id::AtomID( pose.residue_type(cys3).atom_index("V1"), cys3); //Imperfect alignment, but minimization will take care of it.
	core::scoring::superimpose_pose( tbmb_pose, pose, alignment_atoms );

	//Merge the poses and store indices of the new linker residues:
	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const>( pose.conformation_ptr() ) );
	debug_assert( symmconf != nullptr ); //Should be true.
	core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
	debug_assert( symminfo != nullptr ); //Should also be true.
	core::Size const res_per_subunit( symminfo->num_independent_residues() );
	core::Size const tbmb_res1( res_per_subunit + 1 );
	core::Size const tbmb_res2( 2*res_per_subunit + 2 );
	core::Size const tbmb_res3( 3*res_per_subunit + 3 );
	pose.append_residue_by_jump( tbmb_pose.residue(1), cys1 );
	cys2 += 1;
	cys3 += 2;

	//Declare covalent bonds (a bit redundant -- symmetric copy bonds should be handled properly by the symmetry machinery, but we'll declare them just to be sure):
	protocols::cyclic_peptide::DeclareBond bond1, bond2, bond3, bond4, bond5, bond6;
	bond1.set( tbmb_res1, "CM1", cys1, "SG", false );
	bond2.set( tbmb_res2, "CM1", cys2, "SG", false );
	bond3.set( tbmb_res3, "CM1", cys3, "SG", false );
	bond4.set( tbmb_res1, "C2", tbmb_res2, "C1", false );
	bond5.set( tbmb_res2, "C2", tbmb_res3, "C1", false );
	bond6.set( tbmb_res3, "C2", tbmb_res1, "C1", false );
	bond1.apply(pose);
	bond2.apply(pose);
	bond3.apply(pose);
	bond4.apply(pose);
	bond5.apply(pose);
	bond6.apply(pose);
}

/// @brief Given a selection of exactly three residues that have already been connected to a 1,3,5-tris(bromomethyl)benzene crosslinker,
/// add constraints for the crosslinker.
void
TBMB_Helper::add_linker_constraints_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	//Get indices of residues
	core::Size cys1, cys2, cys3;
	ThreefoldLinkerMoverHelper::get_sidechain_indices( selection, cys1, cys2, cys3 );
	core::Size const tbmb_index( get_linker_index_asymmetric( pose, cys1, cys2, cys3) );

	//Set up distance constraints:
	{ //Begin scope
		std::string const dist_cst_string("HARMONIC 0.0 0.01");
		protocols::cyclic_peptide::CreateDistanceConstraint dist_csts;
		utility::vector1< core::Size > res1(6), res2(6);
		utility::vector1< std::string > atom1(6), atom2(6);
		utility::vector1< std::string > cst_fxn(6);
		res1[1] = cys1; res1[2] = cys2; res1[3] = cys3;
		atom2[1] = "V1"; atom2[2] = "V2"; atom2[3] = "V3";
		res1[4] = cys1; res1[5] = cys2; res1[6] = cys3;
		atom2[4] = "CM1"; atom2[5] = "CM2"; atom2[6] = "CM3";
		for ( core::Size i=1; i<=3; ++i ) {
			atom1[i]="SG"; res2[i]=tbmb_index; cst_fxn[i]=dist_cst_string;
		}
		for ( core::Size i=4; i<=6; ++i ) {
			atom1[i]="V1"; res2[i]=tbmb_index; cst_fxn[i]=dist_cst_string;
		}
		if ( TR.Debug.visible() ) {
			TR.Debug << "R1\tA1\tR2\tA2\tFUNC\n";
			for ( core::Size i=1; i<=6; ++i ) {
				TR.Debug << res1[i] << "\t";
				TR.Debug << atom1[i] << "\t";
				TR.Debug << res2[i] << "\t";
				TR.Debug << atom2[i] << "\t";
				TR.Debug << cst_fxn[i] << "\n";
			}
			TR.Debug << std::endl;
		}
		dist_csts.set( res1, atom1, res2, atom2, cst_fxn );
		dist_csts.apply(pose);
	} //End of scope

	//Set up torsion constraints:
	{ //Begin scope
		protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
		utility::vector1 < core::Size > res1(9), res2(9), res3(9), res4(9);
		utility::vector1 < std::string > atom1(9), atom2(9), atom3(9), atom4(9);
		utility::vector1 < std::string > cst_fxns(9);

		//CM#-SG-CB-CA should be three-well potential:
		res1[1] = tbmb_index; res2[1] = cys1; res3[1] = cys1; res4[1] = cys1;
		atom1[1] = "CM1"; atom2[1] = "SG"; atom3[1] = "CB"; atom4[1] = "CA";
		cst_fxns[1] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
		res1[2] = tbmb_index; res2[2] = cys2; res3[2] = cys2; res4[2] = cys2;
		atom1[2] = "CM2"; atom2[2] = "SG"; atom3[2] = "CB"; atom4[2] = "CA";
		cst_fxns[2] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
		res1[3] = tbmb_index; res2[3] = cys3; res3[3] = cys3; res4[3] = cys3;
		atom1[3] = "CM3"; atom2[3] = "SG"; atom3[3] = "CB"; atom4[3] = "CA";
		cst_fxns[3] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

		//C#-C#-CM#-SG should be two-well potential (above/below plane).
		res1[4] = tbmb_index; res2[4] = tbmb_index; res3[4] = tbmb_index; res4[4] = cys1;
		atom1[4] = "C2"; atom2[4] = "C1"; atom3[4] = "CM1"; atom4[4] = "SG";
		cst_fxns[4] = "AMBERPERIODIC 0 2 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
		res1[5] = tbmb_index; res2[5] = tbmb_index; res3[5] = tbmb_index; res4[5] = cys2;
		atom1[5] = "C4"; atom2[5] = "C3"; atom3[5] = "CM2"; atom4[5] = "SG";
		cst_fxns[5] = "AMBERPERIODIC 0 2 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
		res1[6] = tbmb_index; res2[6] = tbmb_index; res3[6] = tbmb_index; res4[6] = cys3;
		atom1[6] = "C6"; atom2[6] = "C5"; atom3[6] = "CM3"; atom4[6] = "SG";
		cst_fxns[6] = "AMBERPERIODIC 0 2 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

		//C#-CM#-SG-CB should be three-well potential:
		res1[7] = tbmb_index; res2[7] = tbmb_index; res3[7] = cys1; res4[7] = cys1;
		atom1[7] = "C1"; atom2[7] = "CM1"; atom3[7] = "SG"; atom4[7] = "CB";
		cst_fxns[7] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
		res1[8] = tbmb_index; res2[8] = tbmb_index; res3[8] = cys2; res4[8] = cys2;
		atom1[8] = "C3"; atom2[8] = "CM2"; atom3[8] = "SG"; atom4[8] = "CB";
		cst_fxns[8] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
		res1[9] = tbmb_index; res2[9] = tbmb_index; res3[9] = cys3; res4[9] = cys3;
		atom1[9] = "C5"; atom2[9] = "CM3"; atom3[9] = "SG"; atom4[9] = "CB";
		cst_fxns[9] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

		if ( TR.Debug.visible() ) {
			TR.Debug << "R1\tA1\tR2\tA2\tR3\tA3\tR4\tA4\tFUNC\n";
			for ( core::Size i=1; i<=9; ++i ) {
				TR.Debug << res1[i] << "\t";
				TR.Debug << atom1[i] << "\t";
				TR.Debug << res2[i] << "\t";
				TR.Debug << atom2[i] << "\t";
				TR.Debug << res3[i] << "\t";
				TR.Debug << atom3[i] << "\t";
				TR.Debug << res4[i] << "\t";
				TR.Debug << atom4[i] << "\t";
				TR.Debug << cst_fxns[i] << "\n";
			}
			TR.Debug << std::endl;
		}

		tors_csts.set( res1, atom1, res2, atom2, res3, atom3, res4, atom4, cst_fxns );
		tors_csts.apply( pose );
	} //End of scope

}

/// @brief Given a selection of exactly three residues that have already been connected to a 1,3,5-tris(bromomethyl)benzene crosslinker,
/// add constraints for the crosslinker.  This version is for symmetric poses.
void
TBMB_Helper::add_linker_constraints_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const linker_was_added
) const {
	//Get indices of residues
	core::Size cys1, cys2, cys3;
	ThreefoldLinkerMoverHelper::get_sidechain_indices( selection, cys1, cys2, cys3 );
	if ( linker_was_added ) {
		cys2 += 1;
		cys3 += 2;
	}
	core::Size tbmb_index1, tbmb_index2, tbmb_index3;
	get_linker_indices_symmetric( pose, cys1, cys2, cys3, tbmb_index1, tbmb_index2, tbmb_index3 );

	//Set up distance constraints:
	for ( core::Size i=1; i<=3; ++i ) {
		core::Size r1, t1, t2;
		if ( i==1 ) { r1=cys1; t1=tbmb_index1; t2=tbmb_index2; }
		else if ( i==2 ) { r1=cys2; t1=tbmb_index2; t2=tbmb_index3; }
		else { r1=cys3; t1=tbmb_index3; t2=tbmb_index1; }
		std::string const dist_cst_string("HARMONIC 0.0 0.01");
		protocols::cyclic_peptide::CreateDistanceConstraint dist_csts;
		utility::vector1< core::Size > res1(6), res2(6);
		utility::vector1< std::string > atom1(6), atom2(6);
		utility::vector1< std::string > cst_fxn(6);

		res1[1] = r1; res2[1] = t1; atom1[1]="SG"; atom2[1] = "V1";
		res1[2] = r1; res2[2] = t1; atom1[2]="V1"; atom2[2] = "CM1";
		res1[3] = t1; res2[3] = t2; atom1[3]="C1"; atom2[3] = "V5";
		res1[4] = t1; res2[4] = t2; atom1[4]="C2"; atom2[4] = "V6";
		res1[5] = t1; res2[5] = t2; atom1[5]="V3"; atom2[5] = "C1";
		res1[6] = t1; res2[6] = t2; atom1[6]="V4"; atom2[6] = "C2";

		for ( core::Size ii=1; ii<=6; ++ii ) { cst_fxn[ii] = dist_cst_string; }

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nSYMM REPEAT " << i << ":\nR1\tA1\tR2\tA2\tFUNC\n";
			for ( core::Size ii=1; ii<=3; ++ii ) {
				TR.Debug << res1[ii] << "\t";
				TR.Debug << atom1[ii] << "\t";
				TR.Debug << res2[ii] << "\t";
				TR.Debug << atom2[ii] << "\t";
				TR.Debug << cst_fxn[ii] << "\n";
			}
			TR.Debug << std::endl;
		}
		dist_csts.set( res1, atom1, res2, atom2, cst_fxn );
		dist_csts.apply(pose);
	} //End of for loop from 1 to 3 for distance constraints.

	//Set up torsion constraints:
	for ( core::Size i=1; i<=3; ++i ) {
		core::Size r1, t1;
		if ( i==1 ) { r1=cys1; t1=tbmb_index1; }
		else if ( i==2 ) { r1=cys2; t1=tbmb_index2; }
		else { r1=cys3; t1=tbmb_index3; }

		protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
		utility::vector1 < core::Size > res1(3), res2(3), res3(3), res4(3);
		utility::vector1 < std::string > atom1(3), atom2(3), atom3(3), atom4(3);
		utility::vector1 < std::string > cst_fxns(3);

		//CM#-SG-CB-CA should be three-well potential:
		res1[1] = t1; res2[1] = r1; res3[1] = r1; res4[1] = r1;
		atom1[1] = "CM1"; atom2[1] = "SG"; atom3[1] = "CB"; atom4[1] = "CA";
		cst_fxns[1] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

		//C#-C#-CM#-SG should be two-well potential (above/below plane).
		res1[2] = t1; res2[2] = t1; res3[2] = t1; res4[2] = r1;
		atom1[2] = "C2"; atom2[2] = "C1"; atom3[2] = "CM1"; atom4[2] = "SG";
		cst_fxns[2] = "AMBERPERIODIC 0 2 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

		//C#-CM#-SG-CB should be three-well potential:
		res1[3] = t1; res2[3] = t1; res3[3] = r1; res4[3] = r1;
		atom1[3] = "C1"; atom2[3] = "CM1"; atom3[3] = "SG"; atom4[3] = "CB";
		cst_fxns[3] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

		if ( TR.Debug.visible() ) {
			TR.Debug << "\nSYMM SET " << i << ":\nR1\tA1\tR2\tA2\tR3\tA3\tR4\tA4\tFUNC\n";
			for ( core::Size ii=1; ii<=3; ++ii ) {
				TR.Debug << res1[ii] << "\t";
				TR.Debug << atom1[ii] << "\t";
				TR.Debug << res2[ii] << "\t";
				TR.Debug << atom2[ii] << "\t";
				TR.Debug << res3[ii] << "\t";
				TR.Debug << atom3[ii] << "\t";
				TR.Debug << res4[ii] << "\t";
				TR.Debug << atom4[ii] << "\t";
				TR.Debug << cst_fxns[ii] << "\n";
			}
			TR.Debug << std::endl;
		}

		tors_csts.set( res1, atom1, res2, atom2, res3, atom3, res4, atom4, cst_fxns );
		tors_csts.apply( pose );
	} //End of scope for torsion constraints.

}

/// @brief Given indices of three cysteine residues that are already linked to a TBMB, get the index
/// of the TBMB residue.
/// @details Throws an error if the three cysteines are not all linked to the same TBMB residue.
core::Size
TBMB_Helper::get_linker_index_asymmetric(
	core::pose::Pose const &pose,
	core::Size const res1,
	core::Size const res2,
	core::Size const res3
) const {
	std::string const errmsg( "Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::get_linker_index_asymmetric(): " );

	core::Size const nconn1( pose.residue(res1).n_possible_residue_connections() );
	core::Size const nconn2( pose.residue(res2).n_possible_residue_connections() );
	core::Size const nconn3( pose.residue(res3).n_possible_residue_connections() );

	core::Size tbmb_index( pose.residue(res1).residue_connection_partner(nconn1) );

	runtime_assert_string_msg( !pose.residue(tbmb_index).name3().compare("TBM"),
		errmsg + "The residue connected to the side-chain of the first cysteine is not TBMB." );
	runtime_assert_string_msg( pose.residue(res2).residue_connection_partner(nconn2) == tbmb_index,
		errmsg + "The residue connected to the side-chain of the first cysteine is not the same as the residue connected to the side-chain of the second cysteine." );
	runtime_assert_string_msg( pose.residue(res3).residue_connection_partner(nconn3) == tbmb_index,
		errmsg + "The residue connected to the side-chain of the first cysteine is not the same as the residue connected to the side-chain of the second cysteine." );

	return tbmb_index;
}

/// @brief Given indices of three cysteine residues that are already linked to pieces of a linker, get
/// of the indices of the symmetric pieces of the linker.
/// @details Throws an error if a residue is not linked to something.  Must be defined by derived classes.
void
TBMB_Helper::get_linker_indices_symmetric(
	core::pose::Pose const &pose,
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size &linker_index1,
	core::Size &linker_index2,
	core::Size &linker_index3
) const {
	std::string const errmsg( "Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::get_linker_indices_symmetric(): " );

	core::Size const nconn1( pose.residue(res1).n_possible_residue_connections() );
	core::Size const nconn2( pose.residue(res2).n_possible_residue_connections() );
	core::Size const nconn3( pose.residue(res3).n_possible_residue_connections() );

	linker_index1 = pose.residue(res1).residue_connection_partner(nconn1) ;
	linker_index2 = pose.residue(res2).residue_connection_partner(nconn2) ;
	linker_index3 = pose.residue(res3).residue_connection_partner(nconn3) ;

	runtime_assert_string_msg( linker_index1 > 0, errmsg + "The first cysteine is not connected to anything!" );
	runtime_assert_string_msg( linker_index2 > 0, errmsg + "The second cysteine is not connected to anything!" );
	runtime_assert_string_msg( linker_index3 > 0, errmsg + "The third cysteine is not connected to anything!" );

	runtime_assert_string_msg( !pose.residue(linker_index1).name3().compare("TBS"),
		errmsg + "The residue connected to the side-chain of the first cysteine is not a symmetric fragment of TBMB." );
	runtime_assert_string_msg( !pose.residue(linker_index2).name3().compare("TBS"),
		errmsg + "The residue connected to the side-chain of the second cysteine is not a symmetric fragment of TBMB." );
	runtime_assert_string_msg( !pose.residue(linker_index3).name3().compare("TBS"),
		errmsg + "The residue connected to the side-chain of the third cysteine is not a symmetric fragment of TBMB." );
}


/// @brief Given a pose with residues selected to be linked by a 1,3,5-tris(bromomethyl)benzene crosslinker,
/// determine whether the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
bool
TBMB_Helper::filter_by_sidechain_distance_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const &filter_multiplier
) const {

	core::Real const hardcoded_cutoff( 12.0 * filter_multiplier );
	core::Real const hardcoded_cutoff_sq( hardcoded_cutoff*hardcoded_cutoff );

	//Get indices of residues
	core::Size cys1, cys2, cys3;
	ThreefoldLinkerMoverHelper::get_sidechain_indices( selection, cys1, cys2, cys3 );

	runtime_assert_string_msg( pose.residue_type(cys1).aa() == core::chemical::aa_cys || pose.residue_type(cys1).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::filter_by_sidechain_distance_asymmetric(): The first residue selected is not D- or L-cysteine." );
	runtime_assert_string_msg( pose.residue_type(cys2).aa() == core::chemical::aa_cys || pose.residue_type(cys2).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::filter_by_sidechain_distance_asymmetric(): The second residue selected is not D- or L-cysteine." );
	runtime_assert_string_msg( pose.residue_type(cys3).aa() == core::chemical::aa_cys || pose.residue_type(cys3).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::filter_by_sidechain_distance_asymmetric(): The third residue selected is not D- or L-cysteine." );

	core::Real const d1( pose.residue(cys1).xyz("CB").distance_squared( pose.residue(cys2).xyz("CB") ) );
	if ( TR.Debug.visible() ) {
		TR.Debug << "D1: " << sqrt(d1) << " MAX: " << hardcoded_cutoff << std::endl;
	}
	if ( d1 > hardcoded_cutoff_sq ) return true;

	core::Real const d2( pose.residue(cys2).xyz("CB").distance_squared( pose.residue(cys3).xyz("CB") ) );
	if ( TR.Debug.visible() ) {
		TR.Debug << "D2: " << sqrt(d2) << " MAX: " << hardcoded_cutoff << std::endl;
	}
	if ( d2 > hardcoded_cutoff_sq ) return true;
	core::Real const d3( pose.residue(cys2).xyz("CB").distance_squared( pose.residue(cys3).xyz("CB") ) );
	if ( TR.Debug.visible() ) {
		TR.Debug << "D3: " << sqrt(d3) << " MAX: " << hardcoded_cutoff << std::endl;
	}
	if ( d3 > hardcoded_cutoff_sq ) return true;

	return false;
}

/// @brief Given a pose with residues selected to be linked by a 1,3,5-tris(bromomethyl)benzene crosslinker,
/// determine whether the residues are too far apart.  This version is for symmetric poses.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
/// @note This version is faster since it only considers the cys1->cys2 distance, and assumes that, due to symmetry,
/// the cys2->cys3 and cys3->cys1 distances must be identical.
bool
TBMB_Helper::filter_by_sidechain_distance_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const &filter_multiplier
) const {

	core::Real const hardcoded_cutoff( 12.0 * filter_multiplier );
	core::Real const hardcoded_cutoff_sq( hardcoded_cutoff*hardcoded_cutoff );

	//Get indices of residues
	core::Size cys1, cys2, cys3;
	ThreefoldLinkerMoverHelper::get_sidechain_indices( selection, cys1, cys2, cys3 );

	runtime_assert_string_msg( pose.residue_type(cys1).aa() == core::chemical::aa_cys || pose.residue_type(cys1).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::filter_by_sidechain_distance_asymmetric(): The first residue selected is not D- or L-cysteine." );
	runtime_assert_string_msg( pose.residue_type(cys2).aa() == core::chemical::aa_cys || pose.residue_type(cys2).aa() == core::chemical::aa_dcs,
		"Error in protocols::cyclic_peptide::threefold_linker::TBMB_Helper::filter_by_sidechain_distance_asymmetric(): The second residue selected is not D- or L-cysteine." );

	core::Real const d1( pose.residue(cys1).xyz("CB").distance_squared( pose.residue(cys2).xyz("CB") ) );
	if ( TR.Debug.visible() ) {
		TR.Debug << "D1: " << sqrt(d1) << " MAX: " << hardcoded_cutoff << std::endl;
	}
	if ( d1 > hardcoded_cutoff_sq ) return true;

	return false;
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
TBMB_Helper::filter_by_constraints_energy_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const &filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, false, false, filter_multiplier );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
TBMB_Helper::filter_by_constraints_energy_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const linker_was_added,
	core::Real const &filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, true, linker_was_added, filter_multiplier );
}

/********************************************
PRIVATE FUNCTIONS
*********************************************/

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
TBMB_Helper::filter_by_constraints_energy(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const symmetric,
	bool const linker_was_added,
	core::Real const &filter_multiplier
) const {
	core::pose::Pose pose_copy( pose ); //Scoring alters the energies in the pose, so we need a copy here.  We're also going to manipulate it a bit.

	//Remove constraints from pose, and add back only those for the linker:
	protocols::simple_moves::ClearConstraintsMover clear_csts;
	clear_csts.apply(pose_copy);
	if ( symmetric ) {
		add_linker_constraints_symmetric( pose_copy, selection, linker_was_added );
	} else {
		add_linker_constraints_asymmetric( pose_copy, selection );
	}

	core::scoring::ScoreFunctionOP sfxn;
	if ( symmetric ) {
		sfxn = core::scoring::ScoreFunctionOP( core::scoring::symmetry::SymmetricScoreFunctionOP ( new core::scoring::symmetry::SymmetricScoreFunction ) );
	} else {
		sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	}

	sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfxn->set_weight( core::scoring::angle_constraint, 1.0 );
	sfxn->set_weight( core::scoring::dihedral_constraint, 1.0 );

	(*sfxn)( pose_copy );

	core::Real const cst_energy( pose_copy.energies().total_energy() );
	bool const failure( cst_energy > (50.0 * filter_multiplier) );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Constraint energy for crosslinked pose is " << cst_energy << ".  " << (failure ? "Filter failed." : "Filter passes.") << std::endl;
	}

	return failure;
}

} //threefold_linker
} //protocols
} //cyclic_peptide
