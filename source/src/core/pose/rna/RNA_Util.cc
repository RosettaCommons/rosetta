// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/rna/RNA_Util.cc
/// @author Rhiju Das

// Unit headers
#include <core/chemical/rna/RNA_Util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/types.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <core/chemical/AA.hh>

// Project headers
#include <numeric/constants.hh>

#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <ObjexxFCL/string.functions.hh>

// Utility headers

// C++

using namespace core::chemical::rna;
static const RNA_FittedTorsionInfo torsion_info;

namespace core {
namespace pose {
namespace rna {

bool
is_cutpoint_open( Pose const & pose, Size const i ) {

	if ( i < 1 ) return true; // user may pass zero -- sometimes checking if we are at a chain terminus, and this would corresponde to n-1 with n=1.

	if ( i >= pose.total_residue() ) return true;

	if ( ! pose.fold_tree().is_cutpoint(i) ) return false;

	if ( pose.residue_type( i   ).has_variant_type( chemical::CUTPOINT_LOWER ) ||
			 pose.residue_type( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) return false;

	return true;
}

//////////////////////////////////////////////////////
bool
is_rna_chainbreak( Pose const & pose, Size const i ) {

	static Real const CHAINBREAK_CUTOFF2 ( 2.5 * 2.5 );

	if ( i >= pose.total_residue() ) return true;
	if ( i < 1 ) return true;

	conformation::Residue const & current_rsd( pose.residue( i   ) ) ;
	conformation::Residue const &    next_rsd( pose.residue( i+1 ) ) ;

	//A little inefficient, since atom indices for these backbone
	// atoms should be the same for all RNA residue types. I think.
	Size atom_O3prime = current_rsd.atom_index( " O3'" );
	Size atom_P      =    next_rsd.atom_index( " P  " );
	Real const dist2 =
		( current_rsd.atom( atom_O3prime ).xyz() - next_rsd.atom( atom_P ).xyz() ).length_squared();

	if ( dist2 > CHAINBREAK_CUTOFF2 ) {
		//std::cout << "Found chainbreak at residue "<< i << " .  O3'-P distance: " << sqrt( dist2 ) << std::endl;
		return true;
	}

	if ( pose.pdb_info() ){
		if ( pose.pdb_info()->number( i ) + 1 != pose.pdb_info()->number( i+1 ) ) return true;
		if ( pose.pdb_info()->chain( i ) != pose.pdb_info()->chain( i+1 ) ) return true;
	}

	return false;

}

////////////////////////////////////////////////////////////////////////////
// Following is quite slow, because it works on a residue in the context
//  of a full pose -- a lot of time wasted on refolding everything.
void
fix_sugar_coords_WORKS_BUT_SLOW(
								 utility::vector1< std::string> atoms_for_which_we_need_new_dofs,
								 utility::vector1< utility::vector1< id::DOF_Type > > which_dofs,
								 utility::vector1< Vector > const & non_main_chain_sugar_coords,
								 Pose & pose,
								 Size const i )
{

	using namespace core::id;

	conformation::Residue const & rsd( pose.residue( i ) );

	//Yup, hard-wired...
	kinematics::Stub const input_stub( rsd.xyz( " C3'" ), rsd.xyz( " C3'" ), rsd.xyz( " C4'" ), rsd.xyz( " C5'" ) );

	utility::vector1< Vector > start_vectors;
	utility::vector1< utility::vector1< Real > > new_dof_sets;

	for (Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
		Size const j = rsd.atom_index( non_main_chain_sugar_atoms[ n ] );
		//Save a copy of starting location
		Vector v = rsd.xyz( non_main_chain_sugar_atoms[ n ]  );
		start_vectors.push_back( v );

		//Desired location
		Vector v2 = input_stub.local2global( non_main_chain_sugar_coords[ n ] );
		pose.set_xyz( id::AtomID( j,i ), v2 );
	}

	for (Size n = 1; n <= atoms_for_which_we_need_new_dofs.size(); n++  ) {
		utility::vector1< Real > dof_set;
		Size const j = rsd.atom_index( atoms_for_which_we_need_new_dofs[ n ] );
		for (Size m = 1; m <= which_dofs[n].size(); m++ ) {
			dof_set.push_back( pose.atom_tree().dof( DOF_ID( AtomID( j,i) , which_dofs[ n ][ m ] ) ) );
		}
		new_dof_sets.push_back( dof_set );
	}


	// Now put the ring atoms in the desired spots, but by changing internal DOFS --
	// rest of the atoms (e.g., in base and 2'-OH) will scoot around as well, preserving
	// ideal bond lengths and angles.
	for (Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
		Size const j = rsd.atom_index( non_main_chain_sugar_atoms[ n ] );
		pose.set_xyz( id::AtomID( j,i ), start_vectors[n] );
	}

	for (Size n = 1; n <= atoms_for_which_we_need_new_dofs.size(); n++  ) {
		Size const j = rsd.atom_index( atoms_for_which_we_need_new_dofs[ n ] );

		for (Size m = 1; m <= which_dofs[n].size(); m++ ) {
			pose.set_dof(  DOF_ID( AtomID( j,i) , which_dofs[ n ][ m ] ), new_dof_sets[ n ][ m ] );
		}

	}
}

/////////////////////////////////////////////////////////
void
prepare_scratch_residue(
		 core::conformation::ResidueOP & scratch_rsd,
		 core::conformation::Residue const & start_rsd,
		 utility::vector1< Vector > const & non_main_chain_sugar_coords,
		 Pose const & pose)
{

	for (Size j = 1; j < scratch_rsd->first_sidechain_atom(); j++ ){
		scratch_rsd->set_xyz( j , start_rsd.xyz( j ) );
	}

	//Yup, hard-wired...
	kinematics::Stub const input_stub( scratch_rsd->xyz( " C3'" ), scratch_rsd->xyz( " C3'" ), scratch_rsd->xyz( " C4'" ), scratch_rsd->xyz( " C5'" ) );

	for (Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
		//Desired location
		Size const j = scratch_rsd->atom_index( non_main_chain_sugar_atoms[ n ] );
		Vector v2 = input_stub.local2global( non_main_chain_sugar_coords[ n ] );
		scratch_rsd->set_xyz( j, v2 );
	}

	Size const o2prime_index( scratch_rsd->atom_index( " O2'" ) );
	scratch_rsd->set_xyz( o2prime_index, scratch_rsd->build_atom_ideal( o2prime_index, pose.conformation() ) );

}


/////////////////////////////////////////////////////////
void
fix_sugar_coords(
								 utility::vector1< std::string> atoms_for_which_we_need_new_dofs,
								 utility::vector1< Vector > const & non_main_chain_sugar_coords,
								 Pose & pose,
								 Pose const & reference_pose,
								 Size const i )
{

	using namespace core::id;
	using namespace core::conformation;

	conformation::Residue const & start_rsd( reference_pose.residue( i ) );

	static ResidueOP scratch_rsd( new Residue( start_rsd.type(), false /*dummy arg*/ ) );

	prepare_scratch_residue( scratch_rsd, start_rsd, non_main_chain_sugar_coords, pose );

	for (Size n = 1; n <= atoms_for_which_we_need_new_dofs.size(); n++  ) {
		Size const j = scratch_rsd->atom_index( atoms_for_which_we_need_new_dofs[ n ] );

		// "Don't do update" --> my hack to prevent lots of refolds. I just want information about whether the
		// atom is a jump_atom, what its stub atoms are, etc... in principle could try to use input_stub_atom1_id(), etc.
		kinematics::tree::AtomCOP current_atom ( & reference_pose.atom_tree().atom_dont_do_update( AtomID(j,i) ) );
		kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );

		if ( !input_stub_atom1) continue;
		if ( (input_stub_atom1->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom1->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

 		Real const d = ( scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ) -
										 scratch_rsd->xyz( (current_atom->id()).atomno() ) ).length();
		pose.set_dof( DOF_ID( AtomID( j, i), D), d );

		if ( input_stub_atom1->is_jump() ) continue;

		kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
		if ( !input_stub_atom2) continue;
		if ( (input_stub_atom2->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom2->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

		Real const theta = numeric::angle_radians(
			scratch_rsd->xyz( (current_atom->id()).atomno() ) ,
			scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom2->id()).atomno() ) );

		pose.set_dof( DOF_ID( AtomID( j, i), THETA), numeric::constants::d::pi - theta );

		// I commented out the following because otherwise, O4' at the 5' end of a pose did not get set properly. (RD, Nov. 2010)
		//  but there may be fallout.
		// if ( input_stub_atom2->is_jump() ) continue; //HEY NEED TO BE CAREFUL HERE.

		kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );

		if ( !input_stub_atom3) continue;
		if ( (input_stub_atom3->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom3->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

		Real const phi = numeric::dihedral_radians(
			scratch_rsd->xyz( (current_atom->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom2->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom3->id()).atomno() ) );

		pose.set_dof( DOF_ID( AtomID( j, i), PHI), phi );

	}
}

//////////////////////////////////////////////////////////////////////////////
void
initialize_atoms_for_which_we_need_new_dofs(
					utility::vector1< std::string > & atoms_for_which_we_need_new_dofs,
					Pose const & pose,  Size const i)
{

	using namespace id;
	using namespace conformation;

	//
	// Which way does atom_tree connectivity flow, i.e. is sugar drawn after base,
	// or after backbone?
	// This is admittedly very ugly, and very RNA specific.
	// ... perhaps this will be figured out in the Cartesian Fragment class?
	//
	conformation::Residue const & rsd( pose.residue( i ) );

 	kinematics::tree::AtomCOP c1prime_atom ( & pose.atom_tree().atom( AtomID( rsd.atom_index( " C1'" ), i ) ) );
 	kinematics::tree::AtomCOP o2prime_atom ( & pose.atom_tree().atom( AtomID( rsd.atom_index( " O2'" ), i ) ) );
 	kinematics::tree::AtomCOP c2prime_atom ( & pose.atom_tree().atom( AtomID( rsd.atom_index( " C2'" ), i ) ) );

	if ( (c1prime_atom->parent()->id()).atomno() == first_base_atom_index( rsd ) ) {
		// There's a jump to this residue.
		//std::cout << "RESIDUE WITH JUMP CONNECTIVITY : " <<  i << std::endl;
		atoms_for_which_we_need_new_dofs.push_back( " C2'" );
		atoms_for_which_we_need_new_dofs.push_back( " C3'" );
		atoms_for_which_we_need_new_dofs.push_back( " O4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C5'" );
		atoms_for_which_we_need_new_dofs.push_back( " O3'" );

	} else if ( (c2prime_atom->parent()->id()).atomno() ==  (o2prime_atom->id()).atomno() ) {

		atoms_for_which_we_need_new_dofs.push_back( " C1'" );
		atoms_for_which_we_need_new_dofs.push_back( " C3'" );
		atoms_for_which_we_need_new_dofs.push_back( " O4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C5'" );
		atoms_for_which_we_need_new_dofs.push_back( " O3'" );

	} else {

		atoms_for_which_we_need_new_dofs.push_back( " C1'" );
		atoms_for_which_we_need_new_dofs.push_back( " C2'" );
		atoms_for_which_we_need_new_dofs.push_back( " O4'" );

	}

}


/////////////////////////////////////////////////////////////////////
//* Passing in a "reference" pose is a dirty trick to prevent
// computationally expensive refolds whenever we change a DOF in the pose.
void
apply_non_main_chain_sugar_coords(
    utility::vector1< Vector > const & non_main_chain_sugar_coords,
		Pose & pose,
		Pose const & reference_pose,
		Size const i
 )
{

	using namespace id;

	/////////////////////////////////////////////
	// Save desired torsion values.
	utility::vector1< Real > start_torsions;
	for (Size j = 1; j <= NUM_RNA_TORSIONS; j++) {
		id::TorsionID rna_torsion_id( i, id::BB, j );
		if ( j > NUM_RNA_MAINCHAIN_TORSIONS) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
		start_torsions.push_back( reference_pose.torsion( rna_torsion_id ) );
	}

	/////////////////////////////////////////////
	//What DOFS do I need to get the ring atoms where I want them?
	utility::vector1< std::string > atoms_for_which_we_need_new_dofs;

	initialize_atoms_for_which_we_need_new_dofs( atoms_for_which_we_need_new_dofs,  pose, i );

	fix_sugar_coords( atoms_for_which_we_need_new_dofs, non_main_chain_sugar_coords, pose, reference_pose, i );

	/////////////////////////////////////////////
	// Reapply desired torsion values.
	for (Size j = 1; j <= NUM_RNA_TORSIONS; j++) {
		id::TorsionID rna_torsion_id( i, id::BB, j );
		if ( j > NUM_RNA_MAINCHAIN_TORSIONS) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
		pose.set_torsion( rna_torsion_id, start_torsions[ j ] );
	}

}

////////////////////////////////////////////////////////////////////
//FANG: All these sugar coord stuffs should be deprecated in favor of
//RNA_IdealCoord class? Are there any performance concern for using
//copy_dof_match_atom_name there?
void
apply_ideal_c2endo_sugar_coords(
		Pose & pose,
		Size const i
 )
{

	//Torsion angles associated with a 2'-endo sugar in 1jj2 (large ribosomal subunit xtal structure ).
	pose.set_torsion( id::TorsionID( i, id::BB, 1), 69.404192 );
	pose.set_torsion( id::TorsionID( i, id::BB, 2), -173.031790 );
	pose.set_torsion( id::TorsionID( i, id::BB, 3), 58.877828 );
	pose.set_torsion( id::TorsionID( i, id::BB, 4), 147.202313 );
	pose.set_torsion( id::TorsionID( i, id::BB, 5), -85.360367 );
	pose.set_torsion( id::TorsionID( i, id::BB, 6), -38.381256 );
	pose.set_torsion( id::TorsionID( i, id::CHI, 1), 111.708846 );
	pose.set_torsion( id::TorsionID( i, id::CHI, 2), -36.423711 );
	pose.set_torsion( id::TorsionID( i, id::CHI, 3), 156.438552 );
	pose.set_torsion( id::TorsionID( i, id::CHI, 4), 179.890442 );

	static utility::vector1< Vector > _non_main_chain_sugar_coords;
	_non_main_chain_sugar_coords.push_back( Vector(  0.329122,   -0.190929,   -1.476983  ) );
	_non_main_chain_sugar_coords.push_back( Vector( -0.783512,   -1.142556,   -1.905737  ) );
	_non_main_chain_sugar_coords.push_back( Vector( -1.928054,   -0.731911,   -1.195034  ) );

	apply_non_main_chain_sugar_coords( _non_main_chain_sugar_coords, pose, pose, i );
}

////////////////////////////////////////////////////////////////////
Size
assign_pucker(
	Pose const & pose,
	Size const rsd_id
) {
	Real const delta = pose.torsion( id::TorsionID( rsd_id, id::BB,  4 ) );
	Size const pucker_state = ( delta < torsion_info.delta_cutoff() ) ? NORTH : SOUTH;
	return pucker_state;
}
////////////////////////////////////////////////////////////////////
void
apply_pucker(
	Pose & pose,
	Size const i,
	Size pucker_state, //0 for using the current pucker
	bool const skip_same_state,
	bool const idealize_coord
) {
	assert( pucker_state <= 2 );

	static const RNA_IdealCoord ideal_coord;
	Real delta, nu1, nu2;

	Size const curr_pucker = assign_pucker( pose, i );
	if ( skip_same_state && pucker_state == curr_pucker ) return;

	if ( pucker_state == WHATEVER ) pucker_state = curr_pucker;

	if ( idealize_coord ) {
		ideal_coord.apply_pucker(pose, i, pucker_state);
	} else {
		if (pucker_state == NORTH) {
			delta = torsion_info.delta_north();
			nu2 = torsion_info.nu2_north();
			nu1 = torsion_info.nu1_north();
		} else {
			delta = torsion_info.delta_south();
			nu2 = torsion_info.nu2_south();
			nu1 = torsion_info.nu1_south();
		}
		pose.set_torsion( id::TorsionID( i, id::BB,  4 ), delta );
		pose.set_torsion( id::TorsionID( i, id::CHI, 2 ), nu2 );
		pose.set_torsion( id::TorsionID( i, id::CHI, 3 ), nu1 );
	}
}
////////////////////////////////////////////////////////////////////

//

} //ns rna
} //ns pose
} //ns core
