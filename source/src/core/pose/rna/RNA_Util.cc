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
#include <core/pose/rna/RNA_Util.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Package headers
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <core/chemical/AA.hh>

// Project headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "core.pose.rna.RNA_Util" );

// Utility headers

// C++

using namespace core::chemical::rna;
static const RNA_FittedTorsionInfo torsion_info;

namespace core {
namespace pose {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////
bool
mutate_position( pose::Pose & pose, Size const i, char const & new_seq ){

	using namespace core::conformation;
	using namespace core::chemical;

	if ( new_seq == pose.sequence()[i-1] ) return false;

	ResidueTypeSet const & rsd_set = pose.residue( i ).residue_type_set();

	ResidueSelector residue_selector;
	residue_selector.set_name1( new_seq );
	residue_selector.match_variants( pose.residue(i).type() );
	if ( pose.residue( i ).is_RNA() ) residue_selector.set_property( "RNA" );
	ResidueTypeCOP new_rsd_type( residue_selector.select( rsd_set )[1] );
	ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );

	Real const save_chi = pose.chi(i);
	pose.replace_residue( i, *new_rsd, false );
	pose.set_chi( i, save_chi );

	return true;
}

//////////////////////////////////////////////////////
void
figure_out_reasonable_rna_fold_tree( pose::Pose & pose )
{
	using namespace core::conformation;

	//Look for chainbreaks in PDB.
	Size const nres = pose.total_residue();
	kinematics::FoldTree f( nres );

	Size m( 0 );

	for (Size i=1; i < nres; ++i) {

		// if ( pose.residue( i ).is_protein() && pose.residue( i+1 ).is_protein() ){
		// 	if ( pose.residue_type( i ).has_variant_type( "C_METHYLAMIDATION" ) ){
		// 		f.new_jump( i, i+1, i );
		// 		m++;
		// 		continue;
		// 	}
		// }

		if ( !pose.residue(i).is_RNA() && !pose.residue(i+1).is_RNA() )  continue;

		if ( (  pose.residue(i).is_RNA() && !pose.residue(i+1).is_RNA() ) ||
				 ( !pose.residue(i).is_RNA() &&  pose.residue(i+1).is_RNA() ) ) {
			f.new_jump( i, i+1, i );
			m++;
			continue;
		}

		if ( pose::rna::is_rna_chainbreak( pose, i ) ){

			//std::cout << "CHAINBREAK between " << i << " and " << i+1 << std::endl;

			f.new_jump( i, i+1, i );
			m++;

			Residue const & current_rsd( pose.residue( i   ) ) ;
			Residue const &    next_rsd( pose.residue( i+1 ) ) ;
			//			Size dummy( 0 ), jump_atom1( 0 ), jump_atom2( 0 );
			//rna_basepair_jump_atoms( current_rsd.aa(), jump_atom1, dummy, dummy );
			//rna_basepair_jump_atoms( next_rsd.aa(), jump_atom2, dummy, dummy );
			//f.set_jump_atoms( m, current_rsd.atom_name( jump_atom1 ), next_rsd.atom_name( jump_atom2 ) );

			f.set_jump_atoms( m,
												chemical::rna::chi1_torsion_atom( current_rsd ),
												chemical::rna::chi1_torsion_atom( next_rsd )   );

		}

	}

	pose.fold_tree( f );
}

////////////////////////////////////////////////////////
void
virtualize_5prime_phosphates( pose::Pose & pose ){

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		if ( i==1 || ( pose.fold_tree().is_cutpoint( i-1 ) &&
									 !pose.residue( i-1 ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
									 !pose.residue( i   ).has_variant_type( chemical::CUTPOINT_UPPER ) ) ){
			if ( pose.residue(i).is_RNA() ) {
				pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", i );
			}
		}

	}

}


//////////////////////////////////////////////////////
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

	static bool const use_phenix_geo = basic::options::option[  basic::options::OptionKeys::rna::corrected_geo ]();
	if ( use_phenix_geo ) {
		apply_pucker( pose, i, SOUTH, false /*skip_same_state*/, true /*idealize_coord*/ );
		return;
	}

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
PuckerState
assign_pucker(
	Pose const & pose,
	Size const rsd_id
) {
	Real const delta = pose.torsion( id::TorsionID( rsd_id, id::BB,  4 ) );
	PuckerState const pucker_state = ( delta < torsion_info.delta_cutoff() ) ? NORTH : SOUTH;
	return pucker_state;
}
////////////////////////////////////////////////////////////////////
void
apply_pucker(
	Pose & pose,
	Size const i,
	PuckerState pucker_state, //0 for using the current pucker
	bool const skip_same_state,
	bool const idealize_coord
) {
	assert( pucker_state <= 2 );

	static const RNA_IdealCoord ideal_coord;
	Real delta, nu1, nu2;

	PuckerState const curr_pucker = assign_pucker( pose, i );
	if ( skip_same_state && pucker_state == curr_pucker ) return;

	if ( pucker_state == core::chemical::rna::ANY_PUCKER ) pucker_state = curr_pucker;

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
	//When a CUTPOINT_UPPER is added to 3' chain_break residue, the EXISTENCE of the CUTPOINT_UPPER atoms means that the alpha torsion which previously DOES NOT exist due to the chain_break now exist. The alpha value is automatically defined to the A-form value by Rosetta. However Rosetta does not automatically adjust the OP2 and OP1 atom position to account for this fact. So it is important that the OP2 and OP1 atoms position are correctly set to be consistent with A-form alpha torsion before the CUTPOINT_UPPER IS ADDED Parin Jan 2, 2009
	void
	correctly_position_cutpoint_phosphate_torsions( pose::Pose & current_pose, Size const five_prime_chainbreak ){

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::io::pdb;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(	core::chemical::RNA );

		chemical::AA res_aa = aa_from_name( "RAD" );
		ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set->aa_map( res_aa )[1] ) ) ;

		Size three_prime_chainbreak = five_prime_chainbreak + 1;
		current_pose.prepend_polymer_residue_before_seqpos( *new_rsd, three_prime_chainbreak, true );
		chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

		//Actually just by prepending the residue causes the alpha torsion to automatically be set to -64.0274,
		//so the manual setting below is actually not needed, May 24, 2010.. Parin S.
		//These are the initial value of virtual upper and lower cutpoint atom.
		//Actually only the alpha (id::BB, 1) is important here since it set the position of O3' (LOWER) atom which in turn determines  OP2 and OP1 atom
		current_pose.set_torsion( TorsionID( three_prime_chainbreak + 1, id::BB, 1 ), -64.027359 );

		/* BEFORE AUG 24, 2011
		//Where the hell did I get these numbers from value...by appending with ideal geometry and look at the initalized value? Oct 13, 2009
		current_pose.set_torsion( TorsionID( five_prime_chainbreak + 1, id::BB, 5 ), -151.943 ); //Not Important?
		current_pose.set_torsion( TorsionID( five_prime_chainbreak + 1, id::BB, 6 ), -76.4185 ); //Not Important?
		current_pose.set_torsion( TorsionID( three_prime_chainbreak + 1, id::BB, 1 ), -64.0274 );
		*/

		//RAD.params
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		//RCY.params
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		//RGU.params
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		//URA.parms
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		current_pose.delete_polymer_residue( five_prime_chainbreak + 1 );
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// try to unify all cutpoint addition into this function.
	void
	correctly_add_cutpoint_variants( core::pose::Pose & pose,
																	 Size const res_to_add,
																	 bool const check_fold_tree /* = true*/){

		using namespace core::chemical;

		runtime_assert( res_to_add < pose.total_residue() );
		if ( check_fold_tree ) runtime_assert( pose.fold_tree().is_cutpoint( res_to_add ) );

		remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS, res_to_add );
		remove_variant_type_from_pose_residue( pose, "THREE_PRIME_PHOSPHATE", res_to_add );
		remove_variant_type_from_pose_residue( pose, "C_METHYLAMIDATION", res_to_add );

		remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS, res_to_add + 1 );
		remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, res_to_add + 1 );
		remove_variant_type_from_pose_residue( pose, "FIVE_PRIME_PHOSPHATE", res_to_add + 1 );
		remove_variant_type_from_pose_residue( pose, "N_ACETYLATION", res_to_add + 1);

		if ( pose.residue_type( res_to_add ).is_RNA() ){
			// could also keep track of alpha, beta, etc.
			runtime_assert( pose.residue_type( res_to_add + 1 ).is_RNA() );
			correctly_position_cutpoint_phosphate_torsions( pose, res_to_add );
		}
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, res_to_add   );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, res_to_add + 1 );

	}


} //ns rna
} //ns pose
} //ns core
