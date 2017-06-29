// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/TMA_Helper.cc
/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the trimesic acid (TMA) cross-linker.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/TMA_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/MoveMap.hh>

// Protocols headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateAngleConstraint.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraint.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.TMA_Helper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
TMA_Helper::TMA_Helper() //:
//TODO initialize data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
TMA_Helper::TMA_Helper( TMA_Helper const &src ) :
	CrosslinkerMoverHelper( src )
	//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
TMA_Helper::~TMA_Helper(){}

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Given a pose and a selection of exactly three residues, add the TMA linker,
/// align it crudely to the selected residues, and set up covalent bonds.
void
TMA_Helper::add_linker_asymmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	core::Size res1, res2, res3;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res1, res2, res3 );

	runtime_assert_string_msg( is_allowed_residue_type( pose.residue_type(res1) ),
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::add_linker_asymmetric(): The first residue selected by the ResidueSelector is not of the following types: " + list_allowed_residue_types() + ".");
	runtime_assert_string_msg( is_allowed_residue_type( pose.residue_type(res2) ),
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::add_linker_asymmetric(): The second residue selected by the ResidueSelector is not of the following types: " + list_allowed_residue_types() + ".");
	runtime_assert_string_msg( is_allowed_residue_type( pose.residue_type(res3) ),
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::add_linker_asymmetric(): The third residue selected by the ResidueSelector is not of the following types: " + list_allowed_residue_types() + ".");

	//Add the sidechain-conjugation type to the position:
	add_sidechain_conjugation( pose, res1 );
	add_sidechain_conjugation( pose, res2 );
	add_sidechain_conjugation( pose, res3 );

	// Make a temporary TMA pose, and place a TMA in it, aligning it crudely
	// to the connector residues.
	core::pose::Pose tma_pose;
	place_tma_asymmetric(tma_pose, pose, res1, res2, res3);

	//Merge the poses:
	pose.append_residue_by_jump( tma_pose.residue(1), res1 );
	core::Size const tma_res( pose.total_residue() );

	//Declare covalent bonds:
	add_linker_bonds_asymmetric( pose, res1, res2, res3, tma_res );
}

/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
/// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
void
TMA_Helper::add_linker_bonds_asymmetric(
	core::pose::Pose &pose,
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size const linker_index
) const {
	//Declare covalent bonds:
	protocols::cyclic_peptide::DeclareBond bond1, bond2, bond3;
	bond1.set( linker_index, "CM1", res1, get_sidechain_amide_name( pose.residue_type(res1) ), false );
	bond2.set( linker_index, "CM2", res2, get_sidechain_amide_name( pose.residue_type(res2) ), false );
	bond3.set( linker_index, "CM3", res3, get_sidechain_amide_name( pose.residue_type(res3) ), false );
	bond1.apply(pose);
	bond2.apply(pose);
	bond3.apply(pose);
}


/// @brief Given a pose and a selection of exactly three residues, add the TMA linker,
/// align it crudely to the selected residues, and set up covalent bonds.
/// @details Version for symmetric poses.
void
TMA_Helper::add_linker_symmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {

	core::Size res1, res2, res3;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res1, res2, res3 );

	runtime_assert_string_msg( is_allowed_residue_type( pose.residue_type(res1) ),
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::add_linker_symmetric(): The first residue selected by the ResidueSelector is not of the following types: " + list_allowed_residue_types() + ".");
	runtime_assert_string_msg( is_allowed_residue_type( pose.residue_type(res2) ),
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::add_linker_symmetric(): The second residue selected by the ResidueSelector is not of the following types: " + list_allowed_residue_types() + ".");
	runtime_assert_string_msg( is_allowed_residue_type( pose.residue_type(res3) ),
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::add_linker_symmetric(): The third residue selected by the ResidueSelector is not of the following types: " + list_allowed_residue_types() + ".");

	//Add the sidechain-conjugation type to the position:
	add_sidechain_conjugation( pose, res1 );

	//Place ASYMMETRIC trimesic acid into a temporary pose (will be used for placing the symmetric).
	core::pose::Pose asymm_tma_pose;
	place_tma_asymmetric( asymm_tma_pose, pose, res1, res2, res3 );

	//Place SYMMETRIC trimesic acid into a temporary pose, and align to the ASYMMETRIC.
	core::pose::Pose symm_tma_pose;
	place_tma_symmetric( symm_tma_pose, asymm_tma_pose );

	//Merge the poses and store indices of the new linker residues:
	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const>( pose.conformation_ptr() ) );
	debug_assert( symmconf != nullptr ); //Should be true.
	core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
	debug_assert( symminfo != nullptr ); //Should also be true.
	core::Size const res_per_subunit( symminfo->num_independent_residues() );
	core::Size const tmares1( res_per_subunit + 1 );
	core::Size const tmares2( 2*res_per_subunit + 2 );
	pose.append_residue_by_jump( symm_tma_pose.residue(1), res1 );
	res2 += 1;
	res3 += 2;

	//Declare covalent bonds:
	add_linker_bonds_symmetric(pose, res1, tmares1, tmares2);

}

/// @brief Given a pose and a TMA linker, add bonds between the TMA and the residues that coordinate the linker.
/// @details Called by add_linker_symmetric().  Version for symmetric poses.
void
TMA_Helper::add_linker_bonds_symmetric(
	core::pose::Pose &pose,
	core::Size const res1,
	core::Size const linker_index1,
	core::Size const linker_index2
) const {
	protocols::cyclic_peptide::DeclareBond bond1, bond2;
	bond1.set( res1, get_sidechain_amide_name(pose.residue_type(res1)), linker_index1, "CM1", false );
	bond2.set( linker_index2, "C1", linker_index1, "C2", false );
	bond1.apply(pose);
	bond2.apply(pose);
}


/// @brief Given a selection of exactly three residues that have already been connected to a trimesic acid crosslinker,
/// add constraints for the crosslinker.
void
TMA_Helper::add_linker_constraints_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	//Get indices of residues
	core::Size res1, res2, res3;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res1, res2, res3 );
	core::Size const tma_index( get_linker_index_asymmetric( pose, res1, res2, res3) );

	//Set up distance constraints:
	{ //Begin scope
		std::string const dist_cst_string("HARMONIC 0.0 0.02");
		protocols::cyclic_peptide::CreateDistanceConstraint dist_csts;
		utility::vector1< core::Size > resindex1(6), resindex2(6);
		utility::vector1< std::string > atom1(6), atom2(6);
		utility::vector1< std::string > cst_fxn(6);

		// 1, 2, and 3: amide nitrogen of sidechain overlaps with TMA virtual atoms.
		resindex1[1] = res1; resindex1[2] = res2; resindex1[3] = res3;
		atom1[1] = get_sidechain_amide_name( pose.residue_type(res1) ); atom1[2] = get_sidechain_amide_name( pose.residue_type(res2) ); atom1[3] = get_sidechain_amide_name( pose.residue_type(res3) );
		atom2[1] = "V1"; atom2[2] = "V2"; atom2[3] = "V3";
		for ( core::Size i=1; i<=3; ++i ) {
			resindex2[i]=tma_index; cst_fxn[i]=dist_cst_string;
		}

		// 4, 5, and 6: sidechain virtual overlaps with TMA CM1,2,3.
		resindex1[4] = res1; resindex1[5] = res2; resindex1[6] = res3;
		atom2[4] = "CM1"; atom2[5] = "CM2"; atom2[6] = "CM3";
		for ( core::Size i=4; i<=6; ++i ) {
			atom1[i]="V1"; resindex2[i]=tma_index; cst_fxn[i]=dist_cst_string;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "R1\tA1\tR2\tA2\tFUNC\n";
			for ( core::Size i=1; i<=6; ++i ) {
				TR.Debug << resindex1[i] << "\t";
				TR.Debug << atom1[i] << "\t";
				TR.Debug << resindex2[i] << "\t";
				TR.Debug << atom2[i] << "\t";
				TR.Debug << cst_fxn[i] << "\n";
			}
			TR.Debug << std::endl;
		}
		dist_csts.set( resindex1, atom1, resindex2, atom2, cst_fxn );
		dist_csts.apply(pose);
	} //End of scope

	{ //Begin scope -- dihedral constraints for TMA chi values.
		std::string const gaussianwidth("0.2"); //Calculate this.
		std::string const gaussianweight( "100.0" ); //Calculate this.
		std::string const cosineweight("50.0"); //Calculate this.
		utility::vector1< std::string > chi_cst_string(5);
		chi_cst_string[1] = "GAUSSIANFUNC 0 " + gaussianwidth + " NOLOG WEIGHT " + gaussianweight;
		chi_cst_string[2] = "GAUSSIANFUNC 3.141592654 " + gaussianwidth + " NOLOG WEIGHT " + gaussianweight;
		chi_cst_string[3] = "GAUSSIANFUNC -3.141592654 " + gaussianwidth + " NOLOG WEIGHT " + gaussianweight;
		chi_cst_string[4] = "AMBERPERIODIC 3.141592654 2.0 " + cosineweight;
		chi_cst_string[5] = "CONSTANTFUNC -31.47852"; //To ensure that the minimum is zero.

		utility::vector1< std::string > atom1(3), atom2(3), atom3(3), atom4(3);
		atom1[1] = "C2"; atom2[1] = "C1"; atom3[1] = "CM1"; atom4[1] = get_sidechain_amide_name( pose.residue_type(res1) );
		atom1[2] = "C2"; atom2[2] = "C3"; atom3[2] = "CM2"; atom4[2] = get_sidechain_amide_name( pose.residue_type(res2) );
		atom1[3] = "C4"; atom2[3] = "C5"; atom3[3] = "CM3"; atom4[3] = get_sidechain_amide_name( pose.residue_type(res3) );

		utility::vector1 < core::Size > r1(3), r2(3), r3(3), r4(3);
		for ( core::Size i=1; i<=3; ++i ) {
			r1[i] = tma_index; r2[i] = tma_index; r3[i] = tma_index;
		}
		r4[1] = res1; r4[2] = res2; r4[3] = res3;

		for ( core::Size i=1; i<=5; ++i ) {
			protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
			utility::vector1< std::string > cst_strings(3);
			cst_strings[1] = cst_strings[2] = cst_strings[3] = chi_cst_string[i];
			tors_csts.set( r1, atom1, r2, atom2, r3, atom3, r4, atom4, cst_strings );
			tors_csts.apply(pose);
		}

	}

	{ //Begin scope -- dihedral constraints for new amide bonds (planar)
		std::string const cst_string( "AMBERPERIODIC 3.141592654 2.0 400.0" );
		utility::vector1< std::string > cst_strings( 3, cst_string );
		utility::vector1< core::Size > r1(3), r2(3), r3(3), r4(3);
		for ( core::Size i=1; i<=3; ++i ) {
			r1[i] = tma_index; r2[i] = tma_index;
		}
		r3[1] = r4[1] = res1;
		r3[2] = r4[2] = res2;
		r3[3] = r4[3] = res3;

		utility::vector1< std::string > atom1(3), atom2(3), atom3(3), atom4(3);
		atom1[1] = "C1"; atom2[1] = "CM1"; atom3[1] = get_sidechain_amide_name( pose.residue_type(res1) ); atom4[1] = get_last_sidechain_carbon_name( pose.residue_type(res1) );
		atom1[2] = "C3"; atom2[2] = "CM2"; atom3[2] = get_sidechain_amide_name( pose.residue_type(res2) ); atom4[2] = get_last_sidechain_carbon_name( pose.residue_type(res2) );
		atom1[3] = "C5"; atom2[3] = "CM3"; atom3[3] = get_sidechain_amide_name( pose.residue_type(res3) ); atom4[3] = get_last_sidechain_carbon_name( pose.residue_type(res3) );

		protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
		tors_csts.set( r1, atom1, r2, atom2, r3, atom3, r4, atom4, cst_strings );
		tors_csts.apply(pose);
	}
}

/// @brief Given a selection of exactly three residues that have already been connected to a trimesic acid crosslinker,
/// add constraints for the crosslinker.  This version is for symmetric poses.
void
TMA_Helper::add_linker_constraints_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const linker_was_added
) const {

	//Get indices of residues
	core::Size res1, res2, res3;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res1, res2, res3 );
	if ( linker_was_added ) {
		res2 += 1;
		res3 += 2;
	}
	core::Size tma_index1, tma_index2, tma_index3;
	get_linker_indices_symmetric( pose, res1, res2, res3, tma_index1, tma_index2, tma_index3 );

	std::string const dist_cst_string("HARMONIC 0.0 0.02");

	//Set up distance constraints (redundantly, repeating for all symmetry repeats)
	for ( core::Size isymm(1); isymm<=3; ++isymm ) { //Begin scope
		core::Size r1, t1, t2;
		if ( isymm==1 ) { r1=res1; t1=tma_index1; t2=tma_index2; }
		else if ( isymm==2 ) { r1=res2; t1=tma_index2; t2=tma_index3; }
		else { r1=res3; t1=tma_index3; t2=tma_index1; }
		protocols::cyclic_peptide::CreateDistanceConstraint dist_csts;
		utility::vector1< core::Size > resindex1(4), resindex2(4);
		utility::vector1< std::string > atom1(4), atom2(4);
		utility::vector1< std::string > cst_fxn(4);

		// 1: amide nitrogen of sidechain overlaps with TMA virtual atoms.
		resindex1[1] = r1;
		resindex2[1] = t1;
		atom1[1] = get_sidechain_amide_name( pose.residue_type(r1) );
		atom2[1] = "V1";
		cst_fxn[1] = dist_cst_string;

		// 2: Virtual on amide nitrogen overlaps with TMA CM1.
		resindex1[2] = r1;
		resindex2[2] = t1;
		atom1[2] = "V1";
		atom2[2] = "CM1";
		cst_fxn[2] = dist_cst_string;

		// 3: Virtual on C2 overlaps with next C1.
		resindex1[3] = t1;
		resindex2[3] = t2;
		atom1[3] = "V3";
		atom2[3] = "C1";
		cst_fxn[3] = dist_cst_string;

		// 4: Virtual on next C1 overlaps with this C2.
		resindex1[4] = t1;
		resindex2[4] = t2;
		atom1[4] = "C2";
		atom2[4] = "V2";
		cst_fxn[4] = dist_cst_string;

		dist_csts.set( resindex1, atom1, resindex2, atom2, cst_fxn );
		dist_csts.apply(pose);
	} //End of scope

	{ //Begin scope -- dihedral constraints for TMA chi values.
		std::string const gaussianwidth("0.2"); //Calculate this.
		std::string const gaussianweight( "100.0" ); //Calculate this.
		std::string const cosineweight("50.0"); //Calculate this.
		utility::vector1< std::string > chi_cst_string(5);
		chi_cst_string[1] = "GAUSSIANFUNC 0 " + gaussianwidth + " NOLOG WEIGHT " + gaussianweight;
		chi_cst_string[2] = "GAUSSIANFUNC 3.141592654 " + gaussianwidth + " NOLOG WEIGHT " + gaussianweight;
		chi_cst_string[3] = "GAUSSIANFUNC -3.141592654 " + gaussianwidth + " NOLOG WEIGHT " + gaussianweight;
		chi_cst_string[4] = "AMBERPERIODIC 3.141592654 2.0 " + cosineweight;
		chi_cst_string[5] = "CONSTANTFUNC -31.47852"; //To ensure that the minimum is zero.
		for ( core::Size isymm(1); isymm<=3; ++isymm ) {

			utility::vector1< std::string > atom1(1), atom2(1), atom3(1), atom4(1);
			atom1[1] = "C2"; atom2[1] = "C1"; atom3[1] = "CM1"; atom4[1] = get_sidechain_amide_name( pose.residue_type(res1) );

			utility::vector1 < core::Size > r1(1), r2(1), r3(1), r4(1);
			if ( isymm == 1 ) { r1[1] = tma_index1; r2[1] = tma_index1; r3[1] = tma_index1; r4[1] = res1; }
			else if ( isymm == 2 ) { r1[1] = tma_index2; r2[1] = tma_index2; r3[1] = tma_index2; r4[1] = res2; }
			else { r1[1] = tma_index3; r2[1] = tma_index3; r3[1] = tma_index3; r4[1] = res3; }

			for ( core::Size i=1; i<=5; ++i ) {
				protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
				utility::vector1< std::string > cst_strings(1);
				cst_strings[1] = chi_cst_string[i];
				tors_csts.set( r1, atom1, r2, atom2, r3, atom3, r4, atom4, cst_strings );
				tors_csts.apply(pose);
			}
		}
	}

	{ //Begin scope -- dihedral constraints for new amide bonds (planar)
		std::string const cst_string( "AMBERPERIODIC 3.141592654 2.0 400.0" );
		utility::vector1< std::string > cst_strings( 3, cst_string );
		utility::vector1< core::Size > r1(3), r2(3), r3(3), r4(3);
		r1[1] = r2[1] = tma_index1;
		r1[2] = r2[2] = tma_index2;
		r1[3] = r2[3] = tma_index3;
		r3[1] = r4[1] = res1;
		r3[2] = r4[2] = res2;
		r3[3] = r4[3] = res3;

		utility::vector1< std::string > atom1(3), atom2(3), atom3(3), atom4(3);
		atom1[1] = "C1"; atom2[1] = "CM1"; atom3[1] = get_sidechain_amide_name( pose.residue_type(res1) ); atom4[1] = get_last_sidechain_carbon_name( pose.residue_type(res1) );
		atom1[2] = "C1"; atom2[2] = "CM1"; atom3[2] = get_sidechain_amide_name( pose.residue_type(res2) ); atom4[2] = get_last_sidechain_carbon_name( pose.residue_type(res2) );
		atom1[3] = "C1"; atom2[3] = "CM1"; atom3[3] = get_sidechain_amide_name( pose.residue_type(res3) ); atom4[3] = get_last_sidechain_carbon_name( pose.residue_type(res3) );

		protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
		tors_csts.set( r1, atom1, r2, atom2, r3, atom3, r4, atom4, cst_strings );
		tors_csts.apply(pose);
	} //End of scope for dihedral constraints for new amide bonds

}

/// @brief Given indices of three residues that are already linked to a TMA, get the index
/// of the TMA residue.
/// @details Throws an error if the three residues are not all linked to the same TMA residue.
core::Size
TMA_Helper::get_linker_index_asymmetric(
	core::pose::Pose const &pose,
	core::Size const res1,
	core::Size const res2,
	core::Size const res3
) const {
	core::Size const sidechain_connect_index_1( get_sidechain_connect_index( pose.residue(res1) ) );
	core::Size const sidechain_connect_index_2( get_sidechain_connect_index( pose.residue(res2) ) );
	core::Size const sidechain_connect_index_3( get_sidechain_connect_index( pose.residue(res3) ) );

	runtime_assert_string_msg( sidechain_connect_index_1 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric(): The first residue does not have a sidechain connection." );
	runtime_assert_string_msg( sidechain_connect_index_2 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric(): The second residue does not have a sidechain connection." );
	runtime_assert_string_msg( sidechain_connect_index_3 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric(): The third residue does not have a sidechain connection." );

	core::Size const tma_index( pose.residue(res1).connected_residue_at_resconn( sidechain_connect_index_1 ) );
	runtime_assert_string_msg( tma_index != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric():  The first sidechain is connected to nothing!" );
	runtime_assert_string_msg(
		pose.residue(res2).connected_residue_at_resconn( sidechain_connect_index_2 ) == tma_index &&
		pose.residue(res2).connected_residue_at_resconn( sidechain_connect_index_2 ) == tma_index,
		"Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric():  The same reside is not connected to all three side-chains."
	);

	return tma_index;
}

/// @brief Given indices of three residues that are already linked to pieces of a linker, get
/// of the indices of the symmetric pieces of the linker.
/// @details Throws an error if a residue is not linked to something.  Must be defined by derived classes.
void
TMA_Helper::get_linker_indices_symmetric(
	core::pose::Pose const &pose,
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size &linker_index1,
	core::Size &linker_index2,
	core::Size &linker_index3
) const {
	core::Size const sidechain_connect_index_1( get_sidechain_connect_index( pose.residue(res1) ) );
	core::Size const sidechain_connect_index_2( get_sidechain_connect_index( pose.residue(res2) ) );
	core::Size const sidechain_connect_index_3( get_sidechain_connect_index( pose.residue(res3) ) );

	runtime_assert_string_msg( sidechain_connect_index_1 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_symmetric(): The first residue does not have a sidechain connection." );
	runtime_assert_string_msg( sidechain_connect_index_2 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_symmetric(): The second residue does not have a sidechain connection." );
	runtime_assert_string_msg( sidechain_connect_index_3 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_symmetric(): The third residue does not have a sidechain connection." );

	linker_index1 = pose.residue(res1).connected_residue_at_resconn( sidechain_connect_index_1 );
	runtime_assert_string_msg( linker_index1 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric():  The first sidechain is connected to nothing!" );
	linker_index2 = pose.residue(res2).connected_residue_at_resconn( sidechain_connect_index_2 );
	runtime_assert_string_msg( linker_index2 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric():  The second sidechain is connected to nothing!" );
	linker_index3 = pose.residue(res3).connected_residue_at_resconn( sidechain_connect_index_3 );
	runtime_assert_string_msg( linker_index3 != 0, "Error in protocols::cyclic_peptide::crosslinker::TMA_Helper::get_linker_index_asymmetric():  The third sidechain is connected to nothing!" );
}


/// @brief Given a pose with residues selected to be linked by a trimesic acid crosslinker,
/// determine whether the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
bool
TMA_Helper::filter_by_sidechain_distance_asymmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const & filter_multiplier
) const {
	//Get indices of residues
	core::Size res1, res2, res3;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res1, res2, res3 );

	core::Real const tma_dist( 8.0 ); //The max distance between atoms bonded to different CM atoms in TMA (slightly padded -- I measure 7.2 A in Avogadro).
	core::Real const maxlen1( get_max_sidechain_length( pose.residue_type(res1) ) );
	core::Real const maxlen2( get_max_sidechain_length( pose.residue_type(res2) ) );
	core::Real const maxlen3( get_max_sidechain_length( pose.residue_type(res3) ) );

	core::Real const one_two_dist( (pose.residue(res1).xyz("CA") - pose.residue(res2).xyz("CA")).length() );
	core::Real const one_two_dist_max( ( maxlen1 + maxlen2 + tma_dist ) * filter_multiplier );
	if ( one_two_dist > one_two_dist_max ) return true;

	core::Real const two_three_dist( (pose.residue(res2).xyz("CA") - pose.residue(res3).xyz("CA")).length() );
	core::Real const two_three_dist_max( ( maxlen2 + maxlen3 + tma_dist ) * filter_multiplier );
	if ( two_three_dist > two_three_dist_max ) return true;

	core::Real const one_three_dist( (pose.residue(res1).xyz("CA") - pose.residue(res3).xyz("CA")).length() );
	core::Real const one_three_dist_max( ( maxlen1 + maxlen3 + tma_dist ) * filter_multiplier );
	if ( one_three_dist > one_three_dist_max ) return true;

	return false; //All checks passed; residures aren't too far apart.
}

/// @brief Given a pose with residues selected to be linked by a trimesic acid crosslinker,
/// determine whether the residues are too far apart.  This version is for symmetric poses.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
/// @note This version is faster since it only considers the cys1->cys2 distance, and assumes that, due to symmetry,
/// the cys2->cys3 and cys3->cys1 distances must be identical.
bool
TMA_Helper::filter_by_sidechain_distance_symmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const & filter_multiplier
) const {
	//Get indices of residues
	core::Size res1, res2, res3;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res1, res2, res3 );

	debug_assert( pose.residue_type(res1).name3() == pose.residue_type(res2).name3() );
	debug_assert( pose.residue_type(res1).name3() == pose.residue_type(res3).name3() );

	core::Real const tma_dist( 8.0 ); //The max distance between atoms bonded to different CM atoms in TMA (slightly padded -- I measure 7.2 A in Avogadro).
	core::Real const maxlen( get_max_sidechain_length( pose.residue_type(res1) ) );

	core::Real const one_two_dist( (pose.residue(res1).xyz("CA") - pose.residue(res2).xyz("CA")).length() );
	core::Real const one_two_dist_max( ( 2*maxlen + tma_dist ) * filter_multiplier );
	if ( one_two_dist > one_two_dist_max ) return true;

	return false; //All checks passed; residures aren't too far apart.
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
TMA_Helper::filter_by_constraints_energy_asymmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const & filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, false, false, 5*filter_multiplier );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
TMA_Helper::filter_by_constraints_energy_symmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const linker_was_added,
	core::Real const & filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, true, linker_was_added, 5*filter_multiplier );
}

/// @brief Optional steps that the helper can apply before every relaxation round.
/// @details Overrides default (doing nothing) to update positions of amide bond dependent atoms.
void
TMA_Helper::pre_relax_round_update_steps(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	bool const /*whole_structure*/,
	bool const symmetric,
	bool const linker_was_added
) const {
	update_tma_amide_bond_dependent_atoms(pose, selection, symmetric, linker_was_added);
}

/// @brief Optional steps that the helper can apply after every relaxation round.
/// @details Overrides default (doing nothing) to update positions of amide bond dependent atoms.
void
TMA_Helper::post_relax_round_update_steps(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	bool const /*whole_structure*/,
	bool const symmetric,
	bool const linker_was_added
) const {
	update_tma_amide_bond_dependent_atoms(pose, selection, symmetric, linker_was_added);
}


/********************************************
PRIVATE FUNCTIONS
*********************************************/

/// @brief Given a residue type, return true for D or L Lys, Orn, DAB, or DAP/DPP, and false otherwise.
/// @details Returns true only for the residue types that I can connect to TMA.
bool
TMA_Helper::is_allowed_residue_type(
	core::chemical::ResidueType const &type
) const {
	//Check for lysine by aa enum:
	if ( type.aa() == core::chemical::aa_lys || type.aa() == core::chemical::aa_dly ) return true;

	//If not lysine, need to do string comparisons:
	std::string const &basetype_name( type.base_name() );

	if ( basetype_name == "ORN" || //L-ornithine
			basetype_name == "DORN" || //D-ornithine
			basetype_name == "DAB" || //L-2,4-diaminobutyric acid
			basetype_name == "DDAB" || //D-2,4-diaminobutyric acid
			basetype_name == "DPP" || //L-2,3-diaminopropionic acid
			basetype_name == "DDPP" //D-2,3-diaminopropionic acid
			) {
		return true;
	}

	//At this point, comparison has failed.
	return false;
}

/// @brief Returns a string listing the allowed residue types that can connect to TMA.
/// @details Intended to be human-readable.  Used for error messages.
std::string
TMA_Helper::list_allowed_residue_types() const {
	return "L- or D-lysine, L- or D-ornithine, L- or D-2,4-diaminobutyric acid, or L- or D-2,3-diaminopropionic acid";
}

/// @brief Given a pose and a position that is one of the allowed types, add the sidechain-conjugation
/// type or patch to the position.
void
TMA_Helper::add_sidechain_conjugation(
	core::pose::Pose &pose,
	core::Size const position
) const {
	debug_assert( is_allowed_residue_type( pose.residue_type( position ) ) ); //Should be guaranteed true by this point.
	protocols::simple_moves::ModifyVariantTypeMover add_sidechain_conjugated_type;
	add_sidechain_conjugated_type.set_additional_type_to_add( "SIDECHAIN_CONJUGATION" );
	core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
	selector->append_index( position );
	add_sidechain_conjugated_type.set_residue_selector(selector);
	add_sidechain_conjugated_type.apply(pose);
}

/// @brief Given a residue type, get the name of the sidechain amide atom.
/// @details Only works for the types allowed to connect to TMA (LYS, ORN, DAP/DPP, DAB and their D-versions).
std::string
TMA_Helper::get_sidechain_amide_name(
	core::chemical::ResidueType const &type
) const {
	debug_assert( is_allowed_residue_type( type ) ); //Should be guaranteed true by this point.

	if ( type.aa() == core::chemical::aa_lys || type.aa() == core::chemical::aa_dly ) {
		return "NZ";
	}

	std::string const &basename( type.base_name() );
	if ( basename == "ORN" || basename == "DORN" ) {
		return "NE";
	}
	if ( basename == "DAB" || basename == "DDAB" ) {
		return "ND";
	}
	if ( basename == "DPP" || basename == "DDPP" ) {
		return "NG";
	}

	return "";
}

/// @brief Given a residue type, get the name of the last carbon in the side chain.
/// @details Only works for the types allowed to connect to TMA (LYS, ORN, DAP/DPP, DAB and their D-versions).
std::string
TMA_Helper::get_last_sidechain_carbon_name(
	core::chemical::ResidueType const &type
) const {
	debug_assert( is_allowed_residue_type( type ) ); //Should be guaranteed true by this point.

	if ( type.aa() == core::chemical::aa_lys || type.aa() == core::chemical::aa_dly ) {
		return "CE";
	}

	std::string const &basename( type.base_name() );
	if ( basename == "ORN" || basename == "DORN" ) {
		return "CD";
	}
	if ( basename == "DAB" || basename == "DDAB" ) {
		return "CG";
	}
	if ( basename == "DPP" || basename == "DDPP" ) {
		return "CB";
	}

	return "";
}

/// @brief Given a residue in a pose, get the index of its first sidechain connection.
///
core::Size
TMA_Helper::get_sidechain_connect_index(
	core::conformation::Residue const &res
) const {
	core::Size returnval(1);
	if ( res.has_lower_connect() ) ++returnval;
	if ( res.has_upper_connect() ) ++returnval;
	if ( returnval > res.n_possible_residue_connections() ) returnval = 0;
	return returnval;
}

/// @brief Given one of the residue types that TMA can connect to, return the max distance between the CA
/// and the NZ.  (Slightly padded).
core::Real
TMA_Helper::get_max_sidechain_length(
	core::chemical::ResidueType const &restype
) const {
	if ( restype.aa() == core::chemical::aa_lys || restype.aa() == core::chemical::aa_dly ) { return 7.0; }
	std::string const &basename( restype.base_name() );
	if ( basename == "ORN" || basename == "DORN" ) { return 6.0; }
	if ( basename == "DAB" || basename == "DDAB" ) { return 5.0; }
	if ( basename == "DPP" || basename == "DDPP" ) { return 4.0; }
	utility_exit_with_message( "Error in protocols::cyclic_peptide::TMA_Helper::get_max_sidechain_length(): The residue type passed to this funcion was not an allowed type that could connect to TMA." );
	return 0.0;
}

/// @brief Given a pose and three connector residue indices, create a new pose with TMA placed somewhere in space that
/// crudely aligns it to the connector residues.
void
TMA_Helper::place_tma_asymmetric(
	core::pose::Pose &tma_pose,
	core::pose::Pose const &pose,
	core::Size const res1,
	core::Size const res2,\
	core::Size const res3
) const {
	tma_pose.clear();

	//Create the TMA and put it into a pose of its own:
	core::chemical::ResidueTypeSetCOP standard_residues( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TMA") ) );
	tma_pose.append_residue_by_jump(*new_rsd, 1);

	//Align the TMA pose to the original pose:
	std::map <core::id::AtomID, core::id::AtomID> alignment_atoms;

	numeric::xyzVector<core::Real> v1( (pose.residue(res1).xyz("CB") - pose.residue(res1).xyz("CA")).normalize() );
	numeric::xyzVector<core::Real> v2( (pose.residue(res2).xyz("CB") - pose.residue(res2).xyz("CA")).normalize() );
	numeric::xyzVector<core::Real> v3( (pose.residue(res3).xyz("CB") - pose.residue(res3).xyz("CA")).normalize () );

	//Make copies of the alignment residues and put them in a new pose for alignment
	core::conformation::ResidueOP rescopy1( pose.residue(res1).clone() );
	core::conformation::ResidueOP rescopy2( pose.residue(res2).clone() );
	core::conformation::ResidueOP rescopy3( pose.residue(res3).clone() );
	rescopy1->set_xyz( get_sidechain_amide_name( pose.residue_type(res1) ), pose.residue(res1).xyz("CA") + v1 * 0.75 * get_max_sidechain_length( pose.residue_type( res1 ) ) );
	rescopy2->set_xyz( get_sidechain_amide_name( pose.residue_type(res2) ), pose.residue(res2).xyz("CA") + v2 * 0.75 * get_max_sidechain_length( pose.residue_type( res2 ) ) );
	rescopy3->set_xyz( get_sidechain_amide_name( pose.residue_type(res3) ), pose.residue(res3).xyz("CA") + v3 * 0.75 * get_max_sidechain_length( pose.residue_type( res3 ) ) );
	rescopy1->set_xyz( "V1", pose.residue(res1).xyz("CA") + v1 * 0.75 * (get_max_sidechain_length( pose.residue_type( res1 ) ) + 1.0) );
	rescopy2->set_xyz( "V1", pose.residue(res2).xyz("CA") + v2 * 0.75 * (get_max_sidechain_length( pose.residue_type( res2 ) ) + 1.0) );
	rescopy3->set_xyz( "V1", pose.residue(res3).xyz("CA") + v3 * 0.75 * (get_max_sidechain_length( pose.residue_type( res3 ) ) + 1.0) );
	core::pose::Pose temp_pose;
	temp_pose.append_residue_by_jump( *rescopy1, 1 );
	temp_pose.append_residue_by_jump( *rescopy2, 1 );
	temp_pose.append_residue_by_jump( *rescopy3, 1 );
	//temp_pose.dump_pdb( "temp_pose_test.pdb" ); //DELETE ME!!!

	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V1"), 1) ] = core::id::AtomID( temp_pose.residue_type(1).atom_index( get_sidechain_amide_name( temp_pose.residue_type(1) ) ), 1);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V2"), 1) ] = core::id::AtomID( temp_pose.residue_type(2).atom_index( get_sidechain_amide_name( temp_pose.residue_type(2) ) ), 2);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V3"), 1) ] = core::id::AtomID( temp_pose.residue_type(3).atom_index( get_sidechain_amide_name( temp_pose.residue_type(3) ) ), 3);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM1"), 1) ] = core::id::AtomID( temp_pose.residue_type(1).atom_index("V1"), 1);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM2"), 1) ] = core::id::AtomID( temp_pose.residue_type(2).atom_index("V1"), 2);
	alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("CM3"), 1) ] = core::id::AtomID( temp_pose.residue_type(3).atom_index("V1"), 3);
	core::scoring::superimpose_pose( tma_pose, temp_pose, alignment_atoms );
}

/// @brief Given a pose containing only an asymmetric TMA molecule, align a symmetric TMA to
/// the first third of the asymmetric TMA, place that in symm_pose.
void
TMA_Helper::place_tma_symmetric(
	core::pose::Pose &symm_pose,
	core::pose::Pose const &asymm_pose
) const {
	using namespace core::id;

	symm_pose.clear();

	//Create the symmetric TMA and put it into a pose of its own:
	core::chemical::ResidueTypeSetCOP standard_residues( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TMASYMM") ) );
	symm_pose.append_residue_by_jump(*new_rsd, 1);

	//Align the symmetric TMA pose to the original pose:
	std::map <AtomID, AtomID> alignment_atoms;
	alignment_atoms[ AtomID(symm_pose.residue_type(1).atom_index("CM1"), 1) ] = AtomID(asymm_pose.residue_type(1).atom_index("CM1"), 1);
	alignment_atoms[ AtomID(symm_pose.residue_type(1).atom_index("C1"), 1) ] = AtomID(asymm_pose.residue_type(1).atom_index("C1"), 1);
	alignment_atoms[ AtomID(symm_pose.residue_type(1).atom_index("C2"), 1) ] = AtomID(asymm_pose.residue_type(1).atom_index("C2"), 1);

	core::scoring::superimpose_pose( symm_pose, asymm_pose, alignment_atoms );
}

/// @brief Given a pose and a selection of three residues that connect to a TMA residue,
/// update the positions of carbonyl oxygen and amide hydrogen atoms in the amide bonds
/// connecting the TMA to the sidechains.
void
TMA_Helper::update_tma_amide_bond_dependent_atoms(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	bool const symmetric,
	bool const tma_was_added
) const {
	if ( symmetric ) {
		update_tma_amide_bond_dependent_atoms_symmetric(pose, selection, tma_was_added);
	} else {
		update_tma_amide_bond_dependent_atoms_asymmetric(pose, selection);
	}
}

/// @brief Given a pose and a selection of three residues that connect to a TMA residue,
/// update the positions of carbonyl oxygen and amide hydrogen atoms in the amide bonds
/// connecting the TMA to the sidechains.
/// @details Version for asymmetric poses.
void
TMA_Helper::update_tma_amide_bond_dependent_atoms_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const &selection
) const {
	runtime_assert_string_msg( !core::pose::symmetry::is_symmetric(pose), "Error in protocols::cyclic_peptide::TMA_Helper::update_tma_amide_bond_dependent_atoms_asymmetric(): The asymmetric version of this function was called on a symmetric pose." );
	utility::vector1< core::Size > res_indices(3);
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices[1], res_indices[2], res_indices[3] );
	core::Size const tma_index( get_linker_index_asymmetric( pose, res_indices[1], res_indices[2], res_indices[3] ) );

	for ( core::Size i=1; i<=3; ++i ) {
		pose.conformation().rebuild_residue_connection_dependent_atoms( tma_index, i ); //Rebuild the carbonyl oxygen on the TMA.
		core::Size const res_connid( pose.residue( tma_index ).residue_connection_conn_id(i) ); //The connection ID on the residue that connects to the TMA
		pose.conformation().rebuild_residue_connection_dependent_atoms( res_indices[i], res_connid ); //Rebuild the amide proton on the side-chain of the residue connecting to the TMA.
	}
}

/// @brief Given a pose and a selection of three residues that connect to a TMA residue,
/// update the positions of carbonyl oxygen and amide hydrogen atoms in the amide bonds
/// connecting the TMA to the sidechains.
/// @details Version for symmetric poses.
void
TMA_Helper::update_tma_amide_bond_dependent_atoms_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	bool const tma_was_added
) const {
	runtime_assert_string_msg( core::pose::symmetry::is_symmetric(pose), "Error in protocols::cyclic_peptide::TMA_Helper::update_tma_amide_bond_dependent_atoms_symmetric(): The symmetric version of this function was called on an asymmetric pose." );
	utility::vector1< core::Size > res_indices(3), tma_index(3);
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices[1], res_indices[2], res_indices[3] );
	if ( tma_was_added ) {
		res_indices[2] += 1;
		res_indices[3] += 2;
	}
	get_linker_indices_symmetric( pose, res_indices[1], res_indices[2], res_indices[3], tma_index[1], tma_index[2], tma_index[3] );

	pose.conformation().rebuild_residue_connection_dependent_atoms( tma_index[1], 1 ); //Rebuild the carbonyl oxygen on the TMA.
	core::Size const res_connid( pose.residue( tma_index[1] ).residue_connection_conn_id(1) ); //The connection ID on the residue that connects to the TMA
	pose.conformation().rebuild_residue_connection_dependent_atoms( res_indices[1], res_connid ); //Rebuild the amide proton on the side-chain of the residue connecting to the TMA.
}

} //crosslinker
} //protocols
} //cyclic_peptide
