// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_Util.cc
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_Util.hh>
#include <core/types.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/kinematics/util.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
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


namespace core {
namespace scoring {
namespace rna {

///////////////////////////////////////////////////////////////////////////////
Size
convert_acgu_to_1234( char const c )
{
	if ( c == 'a') return 1;
	if ( c == 'c') return 2;
	if ( c == 'g') return 3;
	if ( c == 'u') return 4;
	return 0;
}

/////////////////////////////////////////////////////////////////////////////
//This may be used elsewhere -- set up a util.hh?
char get_edge_from_num( Size const num ) {
  if (num == WATSON_CRICK) return 'W';
  if (num == HOOGSTEEN)    return 'H';
  if (num == SUGAR)        return 'S';
  if (num == O2STAR)       return '2';
  if (num == PHOSPHATE)    return 'P';
  return 'X';
}

/////////////////////////////////////////////////////////////////////////////
std::string  //Parin March 7, 2011
get_full_edge_from_num( Size const num ) {
  if (num == WATSON_CRICK) return "WC";
  if (num == HOOGSTEEN)    return "HOOG";
  if (num == SUGAR)        return "SUGAR";
  if (num == O2STAR)       return "O2STAR";
  if (num == PHOSPHATE)    return "PHOS";

	std::cout << "Invalid edge num= " << num << std::endl;
	utility_exit_with_message("Invalid edge num!");
  return "ERROR";

}

/////////////////////////////////////////////////////////////////////////////
//This may be used elsewhere -- set up a util.hh?
char get_orientation_from_num( Size const num ) {
	if (num == 1) return 'A';
	if (num == 2) return 'P';
	return 'X';
}

/////////////////////////////////////////////////////////////////////////////

std::string //Parin March 7, 2011
get_full_orientation_from_num( Size const num ) {
  if (num == 0) return "BLAH";
  if (num == 1) return "ANTI";
  if (num == 2) return "PARA";

	std::cout << "Invalid orientation num= " << num << std::endl;
	utility_exit_with_message("Invalid orientation num!");
  return "ERROR";
}

std::string //Parin April 19, 2011
get_full_LW_orientation_from_num( Size const num ){
  if (num == 0) return "BLAH ";
  if (num == 1) return "CIS  ";
  if (num == 2) return "TRANS";

	std::cout << "Invalid orientation num= " << num << std::endl;
	utility_exit_with_message("Invalid orientation num!");
  return "ERROR";
}

///////////////////////////////////////////////////////////////////////////////
std::string const	first_base_atom( conformation::Residue const & rsd ) {
	//	if (rsd.name1() == 'a' || rsd.name1() == 'g' ) 	return " N9 ";
	//	return " N1 ";
	return rsd.atom_name( first_base_atom_index( rsd ) );
}

///////////////////////////////////////////////////////////////////////////////
bool	is_purine( conformation::Residue const & rsd ) {
	if (rsd.name1() == 'a' || rsd.name1() == 'g' ) 	return true;
	return false;
}


///////////////////////////////////////////////////////////////////////////////
Size first_base_atom_index( conformation::Residue const & rsd ) {
	// HEY MAKE THIS MORE GENERAL? Maybe look at chi1?
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 1 /*chi # 1 must be nucleic acid "chi"*/ );
	return atom_indices[ 3 ]; /* C2' ... C1' ... first base atom ...  chi1 torsion atom*/
}

///////////////////////////////////////////////////////////////////////////////
std::string const	chi1_torsion_atom( conformation::Residue const & rsd ) {
	//	if (rsd.name1() == 'a' || rsd.name1() == 'g' ) 	return " N9 ";
	//	return " N1 ";
	return rsd.atom_name( chi1_torsion_atom_index( rsd ) );
}


///////////////////////////////////////////////////////////////////////////////
Size chi1_torsion_atom_index( conformation::Residue const & rsd ) {
	// HEY MAKE THIS MORE GENERAL? Maybe look at chi1?
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 1 /*chi # 1 must be nucleic acid "chi"*/ );
	return atom_indices[ 4 ]; /* C2' ... C1' ... first base atom ... chi1 torsion atom*/
}


///////////////////////////////////////////////////////////////////////////////
std::string const	default_jump_atom( conformation::Residue const & rsd ) {
	if ( rsd.is_RNA() ){
		if ( !rsd.is_coarse() ){
			return chi1_torsion_atom( rsd );
		} else {
			return " Y  ";
		}
	}
	if ( rsd.name3() == " MG" )		return "MG  ";
	if ( rsd.name3() == "XXX" )	return " Y  ";

	std::cerr << "Residue ??? " << rsd.name3() << std::endl;
	utility_exit_with_message( "Do not know jump atom for this residue" );

	return "????";
}

///////////////////////////////////////////////////////////////////////////////
bool
possibly_canonical( chemical::AA const & aa1,  chemical::AA const & aa2 ) {
	using namespace core::chemical;
	return ( (aa1 == na_rgu && aa2 == na_rcy ) ||
					 (aa1 == na_rcy && aa2 == na_rgu ) ||
					 (aa1 == na_rgu && aa2 == na_ura ) ||
					 (aa1 == na_ura && aa2 == na_rgu ) ||
					 (aa1 == na_rad && aa2 == na_ura ) ||
					 (aa1 == na_ura && aa2 == na_rad )  );
}

///////////////////////////////////////////////////////////////////////////////
//no G-U
bool
possibly_canonical_strict( chemical::AA const & aa1,  chemical::AA const & aa2 ) {
	using namespace core::chemical;
	return ( (aa1 == na_rgu && aa2 == na_rcy ) ||
					 (aa1 == na_rcy && aa2 == na_rgu ) ||
					 (aa1 == na_rad && aa2 == na_ura ) ||
					 (aa1 == na_ura && aa2 == na_rad )  );
}

///////////////////////////////////////////////////////////////////////////////
void
get_watson_crick_base_pair_atoms(
	 chemical::AA const & aa1,
	 chemical::AA const & aa2,
	 std::string & atom1,
	 std::string & atom2 ) {

	using namespace core::chemical;

	if ( aa1==na_rad && aa2==na_ura) {
		atom1 = " N1 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1==na_ura && aa2==na_rad) {
		atom1 = " N3 ";
		atom2 = " N1 ";
		return;
	}
	if ( aa1==na_rgu && aa2==na_rcy) {
		atom1 = " N1 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1==na_rcy && aa2==na_rgu) {
		atom1 = " N3 ";
		atom2 = " N1 ";
		return;
	}
	if ( aa1==na_rgu && aa2==na_ura) {
		atom1 = " O6 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1==na_ura && aa2==na_rgu) {
		atom1 = " N3 ";
		atom2 = " O6 ";
		return;
	}

	atom1 = "XXXX";
	atom2 = "XXXX";
	return;
}

/////////////////////////////////////////////////////////////////////
void
get_watson_crick_base_pair_atoms(
	 chemical::AA const & aa1,
	 chemical::AA const & aa2,
	 utility::vector1< std::string > & atom_ids1,
	 utility::vector1< std::string > & atom_ids2	 )
{

	using namespace chemical;

	atom_ids1.clear();
	atom_ids2.clear();

	if ( aa1==na_rad && aa2==na_ura) {
		atom_ids1.push_back( " N1 ");		atom_ids2.push_back( " H3 ");
		atom_ids1.push_back( " H61");		atom_ids2.push_back( " O4 ");
		return;
	} else if ( aa1==na_rgu && aa2==na_rcy) {
		atom_ids1.push_back( " H1 ");		atom_ids2.push_back( " N3 ");
		atom_ids1.push_back( " H21");		atom_ids2.push_back( " O2 ");
		atom_ids1.push_back( " O6 ");		atom_ids2.push_back( " H41");
		return;
	} else if ( aa1==na_rgu && aa2==na_ura) {
		atom_ids1.push_back( " O6 ");		atom_ids2.push_back( " H3 ");
		atom_ids1.push_back( " H1 ");		atom_ids2.push_back( " O2 ");
		return;
	} else	if ( aa2==na_rad && aa1==na_ura) {
		atom_ids2.push_back( " N1 ");		atom_ids1.push_back( " H3 ");
		atom_ids2.push_back( " H61");		atom_ids1.push_back( " O4 ");
		return;
	} else if ( aa2==na_rgu && aa1==na_rcy) {
		atom_ids2.push_back( " H1 ");		atom_ids1.push_back( " N3 ");
		atom_ids2.push_back( " H21");		atom_ids1.push_back( " O2 ");
		atom_ids2.push_back( " O6 ");		atom_ids1.push_back( " H41");
		return;
	} else if ( aa2==na_rgu && aa1==na_ura) {
		atom_ids2.push_back( " O6 ");		atom_ids1.push_back( " H3 ");
		atom_ids2.push_back( " H1 ");		atom_ids1.push_back( " O2 ");
		return;
	}


}

//////////////////////////////////////////////////////
bool
is_cutpoint_open( core::pose::Pose const & pose, Size const i ) {

	if ( i < 1 ) return true; // user may pass zero -- sometimes checking if we are at a chain terminus, and this would corresponde to n-1 with n=1.

	if ( i >= pose.total_residue() ) return true;

	if ( ! pose.fold_tree().is_cutpoint(i) ) return false;

	if ( pose.residue_type( i   ).has_variant_type( chemical::CUTPOINT_LOWER ) ||
			 pose.residue_type( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) return false;

	return true;
}

//////////////////////////////////////////////////////
bool
is_rna_chainbreak( core::pose::Pose const & pose, Size const i ) {

	static Real const CHAINBREAK_CUTOFF2 ( 2.5 * 2.5 );

	if ( i >= pose.total_residue() ) return true;
	if ( i < 1 ) return true;

	conformation::Residue const & current_rsd( pose.residue( i   ) ) ;
	conformation::Residue const &    next_rsd( pose.residue( i+1 ) ) ;

	//A little inefficient, since atom indices for these backbone
	// atoms should be the same for all RNA residue types. I think.
	Size atom_O3star = current_rsd.atom_index( " O3'" );
	Size atom_P      =    next_rsd.atom_index( " P  " );
	Real const dist2 =
		( current_rsd.atom( atom_O3star ).xyz() - next_rsd.atom( atom_P ).xyz() ).length_squared();

	if ( dist2 > CHAINBREAK_CUTOFF2 ) {
		//std::cout << "Found chainbreak at residue "<< i << " .  O3'-P distance: " << sqrt( dist2 ) << std::endl;
		return true;
	}

	return false;

}

//////////////////////////////////////////////////
utility::vector1< std::string > non_main_chain_sugar_atoms;

void
initialize_non_main_chain_sugar_atoms() {
	static bool init( false );
	if (init) return;

	non_main_chain_sugar_atoms.clear();
	non_main_chain_sugar_atoms.push_back( " C2'" );
	non_main_chain_sugar_atoms.push_back( " C1'" );
	non_main_chain_sugar_atoms.push_back( " O4'" );

	init = true;
}

////////////////////////////////////////////////////////////////////////////
// Following is quite slow, because it works on a residue in the context
//  of a full pose -- a lot of time wasted on refolding everything.
void
fix_sugar_coords_WORKS_BUT_SLOW(
								 utility::vector1< std::string> atoms_for_which_we_need_new_dofs,
								 utility::vector1< utility::vector1< id::DOF_Type > > which_dofs,
								 utility::vector1< Vector > const & non_main_chain_sugar_coords,
								 core::pose::Pose & pose,
								 core::Size const & i )
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
		 core::pose::Pose const & pose)
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

	Size const o2star_index( scratch_rsd->atom_index( " O2'" ) );
	scratch_rsd->set_xyz( o2star_index, scratch_rsd->build_atom_ideal( o2star_index, pose.conformation() ) );

}


/////////////////////////////////////////////////////////
void
fix_sugar_coords(
								 utility::vector1< std::string> atoms_for_which_we_need_new_dofs,
								 utility::vector1< Vector > const & non_main_chain_sugar_coords,
								 core::pose::Pose & pose,
								 core::pose::Pose const & reference_pose,
								 core::Size const & i
								 )
{

	using namespace core::id;
	using namespace core::conformation;
	using namespace core::kinematics;

	conformation::Residue const & start_rsd( reference_pose.residue( i ) );

	static ResidueOP scratch_rsd( new Residue( start_rsd.type(), false /*dummy arg*/ ) );

	prepare_scratch_residue( scratch_rsd, start_rsd, non_main_chain_sugar_coords, pose );

	for (Size n = 1; n <= atoms_for_which_we_need_new_dofs.size(); n++  ) {
		Size const j = scratch_rsd->atom_index( atoms_for_which_we_need_new_dofs[ n ] );

		// "Don't do update" --> my hack to prevent lots of refolds. I just want information about whether the
		// atom is a jump_atom, what its stub atoms are, etc... in principle could try to use input_stub_atom1_id(), etc.
		core::kinematics::tree::AtomCOP current_atom ( & reference_pose.atom_tree().atom_dont_do_update( AtomID(j,i) ) );
		core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );

		if ( !input_stub_atom1) continue;
		if ( (input_stub_atom1->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom1->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

 		Real const d = ( scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ) -
										 scratch_rsd->xyz( (current_atom->id()).atomno() ) ).length();
		pose.set_dof( DOF_ID( AtomID( j, i), D), d );

		if ( input_stub_atom1->is_jump() ) continue;

		core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
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

		core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );

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
					pose::Pose const & pose,  Size const & i
)
{

	using namespace core::id;
	using namespace core::scoring::rna;
	using namespace core::kinematics;
	using namespace conformation;

	//
	// Which way does atom_tree connectivity flow, i.e. is sugar drawn after base,
	// or after backbone?
	// This is admittedly very ugly, and very RNA specific.
	// ... perhaps this will be figured out in the Cartesian Fragment class?
	//
	conformation::Residue const & rsd( pose.residue( i ) );

 	core::kinematics::tree::AtomCOP c1star_atom ( & pose.atom_tree().atom( AtomID( rsd.atom_index( " C1'" ), i ) ) );
 	core::kinematics::tree::AtomCOP o2star_atom ( & pose.atom_tree().atom( AtomID( rsd.atom_index( " O2'" ), i ) ) );
 	core::kinematics::tree::AtomCOP c2star_atom ( & pose.atom_tree().atom( AtomID( rsd.atom_index( " C2'" ), i ) ) );

	if ( (c1star_atom->parent()->id()).atomno() == first_base_atom_index( rsd ) ) {
		// There's a jump to this residue.
		//std::cout << "RESIDUE WITH JUMP CONNECTIVITY : " <<  i << std::endl;
		atoms_for_which_we_need_new_dofs.push_back( " C2'" );
		atoms_for_which_we_need_new_dofs.push_back( " C3'" );
		atoms_for_which_we_need_new_dofs.push_back( " O4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C5'" );
		atoms_for_which_we_need_new_dofs.push_back( " O3'" );

	} else if ( (c2star_atom->parent()->id()).atomno() ==  (o2star_atom->id()).atomno() ) {

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
		core::pose::Pose & pose,
		core::pose::Pose const & reference_pose,
		core::Size const & i
 )
{

	using namespace core::id;
	using namespace core::scoring::rna;

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
	initialize_non_main_chain_sugar_atoms();
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
void
apply_ideal_c2endo_sugar_coords(
		core::pose::Pose & pose,
		core::Size const & i
 )
{

	//Torsion angles associated with a 2'-endo sugar in 1jj2 (large ribosomal subunit xtal structure ).
	pose.set_torsion( core::id::TorsionID( i, id::BB, 1), 69.404192 );
	pose.set_torsion( core::id::TorsionID( i, id::BB, 2), -173.031790 );
	pose.set_torsion( core::id::TorsionID( i, id::BB, 3), 58.877828 );
	pose.set_torsion( core::id::TorsionID( i, id::BB, 4), 147.202313 );
	pose.set_torsion( core::id::TorsionID( i, id::BB, 5), -85.360367 );
	pose.set_torsion( core::id::TorsionID( i, id::BB, 6), -38.381256 );
	pose.set_torsion( core::id::TorsionID( i, id::CHI, 1), 111.708846 );
	pose.set_torsion( core::id::TorsionID( i, id::CHI, 2), -36.423711 );
	pose.set_torsion( core::id::TorsionID( i, id::CHI, 3), 156.438552 );
	pose.set_torsion( core::id::TorsionID( i, id::CHI, 4), 179.890442 );

	utility::vector1< Vector > non_main_chain_sugar_coords;
	non_main_chain_sugar_coords.push_back( Vector(  0.329122,   -0.190929,   -1.476983  ) );
	non_main_chain_sugar_coords.push_back( Vector( -0.783512,   -1.142556,   -1.905737  ) );
	non_main_chain_sugar_coords.push_back( Vector( -1.928054,   -0.731911,   -1.195034  ) );

	apply_non_main_chain_sugar_coords( non_main_chain_sugar_coords, pose, pose, i );


}

///////////////////////////////////////////////////////////////////////////////
// Simple cubic spline.
void
get_fade_correction(
   Real const z,
	 Real const cutoff_lower,
	 Real const cutoff_upper,
	 Real const fade_zone,
	 Real & fade_value,
	 Real & fade_deriv )
{
	assert( fade_zone > 0 );

	fade_value = 1.0;
	fade_deriv = 0.0;

	if (z < cutoff_lower || z > cutoff_upper ){
		fade_value = 0.0;
	} else if ( z < cutoff_lower + fade_zone ) {
		//Check little strip near lower cutoff.
		Real const b = -1.0 * ( z - (cutoff_lower + fade_zone) )/ fade_zone;
		Real const b2 = b*b;
		Real const b3 = b2*b;
		fade_value = ( 2 * b3 - 3 * b2 + 1 );
		fade_deriv = -1.0 * (6 * b2 - 6 * b ) / fade_zone;
	} else if ( z > cutoff_upper - fade_zone ) {
		//Check little strip near upper cutoff.
		Real const b =  ( z - (cutoff_upper - fade_zone) )/ fade_zone;
		Real const b2 = b*b;
		Real const b3 = b2*b;
		fade_value = ( 2 * b3 - 3 * b2 + 1 );
		fade_deriv = (6 * b2 - 6 * b ) / fade_zone;
	}

	return;

}

///////////////////////////////////////////////////////////////////////////////////////////////
//Unify the version in StepWiseRNA_Utill.cc and RNA_CentroidInfo.cc on June 25, 2011
// Also, this copies some code from Phil's dna/base_geometry.cc
//Comments (Parin Sep 23 ,2009)...possible problem if every atoms in the nucleotide is virtual...in that case numatoms=0....will this crash the code??

Vector
get_rna_base_centroid( conformation::Residue const & rsd , bool verbose){

  //SML PHENIX conference
	if ( !rsd.is_RNA() ) {
		if (basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value()){
			return Vector( 0.0, 0.0, 0.0 );
		} else { //if not option
			utility_exit_with_message("non-RNA residue inside get_rna_base_centroid");
		}
	}//if not RNA

  Vector centroid( 0.0 );
  Size numatoms = 0;

	//Consistency check:
	//if(rsd.type().atom_name(rsd.first_sidechain_atom()) !=" O2'") utility_exit_with_message( "rsd.type().atom_name(rsd.first_sidechain_atom()) !=\" O2'\" " );
	//if(rsd.atom_name( rsd.first_sidechain_atom() )!=" O2'") utility_exit_with_message("rsd.atom_name( rsd.first_sidechain_atom() )!=\" O2'\"");

	if( rsd.RNA_type().o2star_index()!=rsd.first_sidechain_atom() ){
		utility_exit_with_message( "rsd.RNA_info().o2star_index()!=rsd.first_sidechain_atom()");
	}

	if(verbose)  std::cout << "Base atoms" << std::endl;

	for ( Size i=rsd.first_sidechain_atom()+1; i<= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2star oxygen.

		if(verbose) std::cout << "atom " << i  << " " << 	"name= " << rsd.type().atom_name(i) << " type= " << rsd.atom_type(i).name()  << " " << rsd.atom_type_index(i) << " " << rsd.atomic_charge(i);

		if(rsd.RNA_type().atom_is_virtual(i)){
			if(verbose) std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}

		if(verbose) std::cout << std::endl;

    centroid += rsd.xyz(i);
    numatoms++;
  }

	if(numatoms==0){//Centroid not well defined in this case...probably because rsd is a virtual residue...just return 0
		Vector dummy_centroid( 0.0 );
		return dummy_centroid;
	}

  centroid /= static_cast< Real >( numatoms );

  return centroid;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Unify the version in StepWiseRNA_Utill.cc and RNA_CentroidInfo.cc on June 25, 2011
numeric::xyzMatrix< core::Real >
get_rna_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ){

	using namespace chemical;


  //SML PHENIX conference
  if ( !rsd.is_RNA() ) {
    if (basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value()){
      return numeric::xyzMatrix< core::Real >::identity();
    } else { //if not option
      utility_exit_with_message("non-RNA residue inside get_rna_base_coordinate_system, abort");
    }
  }//if not RNA

 	Size res_type = rsd.aa();

	Vector x,y,z;

	// Make an axis pointing from base centroid to Watson-Crick edge.
	std::string WC_atom;
	if ( res_type == na_rad ) WC_atom = " N1 ";
	if ( res_type == na_rcy ) WC_atom = " N3 ";
	if ( res_type == na_rgu ) WC_atom = " N1 ";
	if ( res_type == na_ura ) WC_atom = " N3 ";

	Vector const WC_coord (rsd.xyz( WC_atom ) );
	x = WC_coord - centroid;
	x.normalize();

	// Make a perpendicular axis pointing from centroid towards
	// Hoogstein edge (e.g., major groove in a double helix).
	std::string H_atom;
	if ( res_type == na_rad ) H_atom = "N7";
	if ( res_type == na_rcy ) H_atom = "C5";
	if ( res_type == na_rgu ) H_atom = "N7";
	if ( res_type == na_ura ) H_atom = "C5";

  Vector const H_coord (rsd.xyz( H_atom ) );
	y = H_coord - centroid; //not orthonormal yet...
	z = cross(x, y);
	z.normalize(); // Should poSize roughly 5' to 3' if in a double helix.

	y = cross(z, x);
	y.normalize(); //not necessary but doesn't hurt.


 	numeric::xyzMatrix< core::Real > M=numeric::xyzMatrix< core::Real >::cols( x, y, z );
	return M;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Check whether one of the atom beyond to a base and another to a phosphate. Doesn't have to be on the same nucleotide.
bool
Is_base_phosphate_atom_pair( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const atomno_1, Size const atomno_2){

	bool Is_base_phosphate_atom_pair=false;

	if( ( rsd_1.RNA_type().atom_is_phosphate( atomno_1 ) && (rsd_2.RNA_type().is_RNA_base_atom( atomno_2 ) ) ) ) Is_base_phosphate_atom_pair=true;
	if( ( rsd_2.RNA_type().atom_is_phosphate( atomno_2 ) && (rsd_1.RNA_type().is_RNA_base_atom( atomno_1 ) ) ) ) Is_base_phosphate_atom_pair=true;

	if(Is_base_phosphate_atom_pair){ //This Assume that rsd_1 and rsd_2 are the same!!!
		if( rsd_1.seqpos()==rsd_2.seqpos() && (rsd_1.path_distance( atomno_1, atomno_2 ) < 4) ){ //consistency check!
			utility_exit_with_message("Is_base_phosphate_atom_pair but rsd.path_distance( " + ObjexxFCL::string_of(atomno_1) + " , " + ObjexxFCL::string_of(atomno_2) + " ) < 4");
		}
	}

	return Is_base_phosphate_atom_pair;

}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< std::string > const &
get_atoms_involved_in_phosphate_torsion()
{
	static utility::vector1< std::string > atoms_involved_in_phosphate_torsion;
	static bool init( false );

	if ( !init ){

		atoms_involved_in_phosphate_torsion.clear();

		atoms_involved_in_phosphate_torsion.push_back( " P  " );
		atoms_involved_in_phosphate_torsion.push_back( " OP2" );
		atoms_involved_in_phosphate_torsion.push_back( " OP1" );
		atoms_involved_in_phosphate_torsion.push_back( " O5'" );
		atoms_involved_in_phosphate_torsion.push_back( " H5'" );
		atoms_involved_in_phosphate_torsion.push_back( "H5''" );

		init = true;

	}

	return atoms_involved_in_phosphate_torsion;

}


} //ns rna
} //ns scoring
} //ns core
