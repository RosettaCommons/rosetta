// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/RNA_Util.cc
/// @author Rhiju Das

// Unit headers
#include <core/pose/rna/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/types.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularSplineFunc.hh>
#include <core/scoring/methods/chainbreak_util.hh>
#include <core/scoring/rna/RNA_Motif.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>

// Project headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <utility/string_util.hh>

static basic::Tracer TR( "core.pose.rna.RNA_Util" );

// Utility headers

// C++

using namespace core::chemical::rna;
using namespace core::scoring::rna;
using utility::strip;
static const RNA_FittedTorsionInfo torsion_info;

namespace core {
namespace pose {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////
bool
mutate_position( pose::Pose & pose, Size const i, char const new_seq ) {

	using namespace core::conformation;
	using namespace core::chemical;

	if ( new_seq == pose.sequence()[i-1] ) return false;

	ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose( pose.residue_type( i ).mode() );

	ResidueProperty base_property = ( pose.residue_type( i ).is_RNA() ) ? RNA : NO_PROPERTY;
	ResidueTypeCOP new_rsd_type( ResidueTypeFinder( *rsd_set ).name1( new_seq ).variants( pose.residue_type(i).variant_types() ).base_property( base_property ).get_representative_type() );

	ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );

	Real const save_chi = pose.chi(i);
	pose.replace_residue( i, *new_rsd, false );
	pose.set_chi( i, save_chi );

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////
bool
mutate_position( pose::Pose & pose, Size const i, chemical::ResidueType const & rt ) {

	using namespace core::conformation;
	using namespace core::chemical;

	if ( rt.name() == pose.residue_type( i ).name() ) return false;

	ResidueOP new_rsd( ResidueFactory::create_residue( rt, pose.residue( i ), pose.conformation() ) );

	Real const save_chi = pose.chi(i);
	pose.replace_residue( i, *new_rsd, false );
	pose.set_chi( i, save_chi );

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////
bool
mutate_position( pose::Pose & pose, Size const i, std::string const & name3 ) {

	using namespace core::conformation;
	using namespace core::chemical;

	std::string clean_current;
	std::map< Size, std::string > special_res;
	remove_and_store_bracketed( pose.annotated_sequence(), clean_current, special_res );

	if ( name3 == special_res[i-1] ) return false;

	ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose( pose.residue_type( i ).mode() );

	ResidueProperty base_property = ( pose.residue_type( i ).is_RNA() ) ? RNA : NO_PROPERTY;
	ResidueTypeCOP new_rsd_type( ResidueTypeFinder( *rsd_set ).name3( name3 ).variants( pose.residue_type(i).variant_types() ).base_property( base_property ).get_representative_type() );

	ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );

	Real const save_chi = pose.chi(i);
	pose.replace_residue( i, *new_rsd, false );
	pose.set_chi( i, save_chi );

	return true;
}

//////////////////////////////////////////////////////
void
figure_out_reasonable_rna_fold_tree(
	pose::Pose & pose,
	bool force_cut_at_rna_chainbreak /* = false */  )
{
	using namespace core::conformation;
	using namespace core::chemical;

	//Look for chainbreaks in PDB.
	Size const nres = pose.size();
	kinematics::FoldTree f( nres );

	Size m( 0 );

	for ( Size i = 1; i < nres; ++i ) {
		bool new_jump( false );
		if ( (  pose.residue_type(i).is_RNA() != pose.residue_type(i+1).is_RNA() ) ||
				(  pose.residue_type(i).is_protein() != pose.residue_type(i+1).is_protein() ) ||
				(  pose.pdb_info() && ( pose.pdb_info()->chain( i ) != pose.pdb_info()->chain( i+1 ) ) ) ||
				(  pose.residue_type(i).has_variant_type( CUTPOINT_LOWER ) && pose.residue_type(i+1).has_variant_type( CUTPOINT_UPPER ) ) ) {
			new_jump = true;
		}

		// probably would be OK to also check for protein or general polymer chainbreak based on UPPER/LOWER (util function may exist somewhere?).
		bool rna_chainbreak( false );
		if ( pose.residue_type(i).is_RNA() && pose.residue_type(i+1).is_RNA() ) rna_chainbreak =  pose::rna::is_rna_chainbreak( pose, i );

		if ( !new_jump &&
				pose.pdb_info() && ( pose.pdb_info()->number(i+1) != pose.pdb_info()->number(i)+1 ) &&
				!basic::options::option[ basic::options::OptionKeys::full_model::allow_jump_in_numbering ]() ) {
			if ( rna_chainbreak ) {
				new_jump = true;
			} else {
				// user may not realize that jump in numbering is confusing.
				//    TR << "There appears to be a break in numbering but not chain geometry at " << i << " to " << i+1;
				//    if ( pose.pdb_info() ) TR << " [in PDBinfo numbering: "  << pose.pdb_info()->chain(i) << ":" << pose.pdb_info()->number(i) << " to " << pose.pdb_info()->chain(i+1) << ":" << pose.pdb_info()->number(i+1) << "]" << std::endl;
				//    TR << " so adding to cutpoint list. If that is wrong, use flag -allow_jump_in_numbering, and be careful with other parts of Rosetta." << TR.Reset << std::endl;
			}
		}

		// Following is a legacy condition -- back when we used to make RNA poses that had all the same chain.
		if ( !new_jump && rna_chainbreak ) {
			if ( force_cut_at_rna_chainbreak ||  basic::options::option[ basic::options::OptionKeys::rna::cut_at_rna_chainbreak ]() ) {
				new_jump = true;
			} else {
				TR << TR.Red << "Found possible chainbreak at: " << i << " to " << i+1;
				if ( pose.pdb_info() ) TR << " [in PDBinfo numbering: "  << pose.pdb_info()->chain(i) << ":" << pose.pdb_info()->number(i) << " to " << pose.pdb_info()->chain(i+1) << ":" << pose.pdb_info()->number(i+1) << "]" << std::endl;
				TR << " but not enforcing cutpoint. If that is wrong, use flag -cut_at_rna_chainbreak; or, preferably, change PDB number or chains to reflect cutpoint." << std::endl;
			}
		}

		if ( new_jump ) {
			f.new_jump( i, i+1, i );
			m++;

			if ( pose.residue_type(i).is_RNA() && pose.residue_type(i+1).is_RNA() ) {
				ResidueType const & current_rsd( pose.residue_type( i   ) ) ;
				ResidueType const &    next_rsd( pose.residue_type( i+1 ) ) ;
				f.set_jump_atoms( m,
					chemical::rna::chi1_torsion_atom( current_rsd ),
					chemical::rna::chi1_torsion_atom( next_rsd )   );
			}
		}
	}

	pose.fold_tree( f );
}

////////////////////////////////////////////////////////
void
virtualize_5prime_phosphates( pose::Pose & pose ){

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( i==1 || ( pose.fold_tree().is_cutpoint( i-1 ) &&
				!pose.residue_type( i-1 ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
				!pose.residue_type( i   ).has_variant_type( chemical::CUTPOINT_UPPER )
				) ) {
			if ( pose.residue_type(i).is_RNA() && !pose.residue_type(i).has_variant_type( chemical::FIVEPRIME_CAP ) ) {
				pose::add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_PHOSPHATE, i );
			}
		}
	}
}


//////////////////////////////////////////////////////
bool
is_cutpoint_open( Pose const & pose, Size const i ) {

	if ( i < 1 ) return true; // user may pass zero -- sometimes checking if we are at a chain terminus, and this would corresponde to n-1 with n=1.
	if ( i >= pose.size() ) return true;

	if ( ! pose.fold_tree().is_cutpoint(i) ) return false;
	if ( pose.residue_type( i   ).has_variant_type( chemical::CUTPOINT_LOWER ) ||
			pose.residue_type( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) return false;

	return true;
}

//////////////////////////////////////////////////////
bool
is_rna_chainbreak( Pose const & pose, Size const i ) {

	static Real const CHAINBREAK_CUTOFF2 ( 2.5 * 2.5 );

	if ( i >= pose.size() ) return true;
	if ( i < 1 ) return true;

	conformation::Residue const & current_rsd( pose.residue( i   ) ) ;
	conformation::Residue const &    next_rsd( pose.residue( i+1 ) ) ;
	runtime_assert ( current_rsd.is_RNA() );
	runtime_assert ( next_rsd.is_RNA() );

	Size atom_O3prime = current_rsd.type().RNA_info().o3prime_atom_index();
	Size atom_P       =    next_rsd.type().RNA_info().p_atom_index();
	Real const dist2 =
		( current_rsd.atom( atom_O3prime ).xyz() - next_rsd.atom( atom_P ).xyz() ).length_squared();

	if ( dist2 > CHAINBREAK_CUTOFF2 ) {
		TR.Debug << "Found chainbreak at residue "<< i << " .  O3'-P distance: " << sqrt( dist2 ) << std::endl;
		return true;
	}

	if ( pose.pdb_info() ) {
		if ( pose.pdb_info()->number( i ) + 1 != pose.pdb_info()->number( i+1 ) ) return true;
		if ( pose.pdb_info()->chain( i ) != pose.pdb_info()->chain( i+1 ) ) return true;
	}

	return false;
}

/////////////////////////////////////////////////////////
void
prepare_scratch_residue(
	core::conformation::ResidueOP & scratch_rsd,
	core::conformation::Residue const & start_rsd,
	utility::vector1< Vector > const & non_main_chain_sugar_coords,
	Pose const & pose)
{
	// Can't do anything like this if this isn't an RNA residue.
	if ( ! scratch_rsd->type().is_RNA() ) return;

	for ( Size j = 1; j < scratch_rsd->first_sidechain_atom(); j++ ) {
		scratch_rsd->set_xyz( j , start_rsd.xyz( j ) );
	}

	kinematics::Stub const input_stub( scratch_rsd->xyz( scratch_rsd->type().RNA_info().c3prime_atom_index() ),
		scratch_rsd->xyz( scratch_rsd->type().RNA_info().c3prime_atom_index() ),
		scratch_rsd->xyz( scratch_rsd->type().RNA_info().c4prime_atom_index() ),
		scratch_rsd->xyz( scratch_rsd->type().RNA_info().c5prime_atom_index() ) );

	for ( Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
		//Desired location
		Size const j = scratch_rsd->atom_index( non_main_chain_sugar_atoms[ n ] );
		Vector v2 = input_stub.local2global( non_main_chain_sugar_coords[ n ] );
		scratch_rsd->set_xyz( j, v2 );
	}

	Size const o2prime_index( scratch_rsd->type().RNA_info().o2prime_index() );
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

	for ( auto const & atom_name : atoms_for_which_we_need_new_dofs ) {
		Size const j = scratch_rsd->atom_index( atom_name );

		AtomID const aid( j, i );
		// "Don't do update" --> my hack to prevent lots of refolds. I just want information about whether the
		// atom is a jump_atom, what its stub atoms are, etc... in principle could try to use input_stub_atom1_id(), etc.
		kinematics::tree::AtomCOP current_atom ( reference_pose.atom_tree().atom_dont_do_update( aid ).get_self_ptr() );
		kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );

		if ( !input_stub_atom1 ) continue;
		if ( (input_stub_atom1->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom1->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

		Real const d = ( scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ) -
			scratch_rsd->xyz( (current_atom->id()).atomno() ) ).length();
		pose.set_dof( DOF_ID( aid, D), d );

		if ( input_stub_atom1->is_jump() ) continue;

		kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
		if ( !input_stub_atom2 ) continue;
		if ( (input_stub_atom2->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom2->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

		Real const theta = numeric::angle_radians(
			scratch_rsd->xyz( (current_atom->id()).atomno() ) ,
			scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom2->id()).atomno() ) );

		pose.set_dof( DOF_ID( aid, THETA), numeric::constants::d::pi - theta );

		// I commented out the following because otherwise, O4' at the 5' end of a pose did not get set properly. (RD, Nov. 2010)
		//  but there may be fallout.
		// if ( input_stub_atom2->is_jump() ) continue; //HEY NEED TO BE CAREFUL HERE.

		kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );

		if ( !input_stub_atom3 ) continue;
		if ( (input_stub_atom3->id()).rsd() != (current_atom->id()).rsd() ) continue;
		if ( (input_stub_atom3->id()).atomno() > scratch_rsd->first_sidechain_atom() ) continue;

		Real const phi = numeric::dihedral_radians(
			scratch_rsd->xyz( (current_atom->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom1->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom2->id()).atomno() ),
			scratch_rsd->xyz( (input_stub_atom3->id()).atomno() ) );

		pose.set_dof( DOF_ID( aid, PHI), phi );
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
	using namespace chemical;

	// Which way does atom_tree connectivity flow, i.e. is sugar drawn after base,
	// or after backbone?
	// This is admittedly very ugly, and very RNA specific.
	// ... perhaps this will be figured out in the Cartesian Fragment class?
	//
	ResidueType const & rsd( pose.residue_type( i ) );
	RNA_Info const & rna_type( rsd.RNA_info() );

	// AMW: We don't need to ask ALL of these by RNA_info -- note
	// that they're using strings anyway. The key is since some residues don't
	// have O4' but S4'.

	kinematics::tree::AtomCOP c1prime_atom ( pose.atom_tree().atom( AtomID( rna_type.c1prime_atom_index(), i ) ).get_self_ptr() );
	kinematics::tree::AtomCOP o2prime_atom ( pose.atom_tree().atom( AtomID( rna_type.o2prime_index(), i ) ).get_self_ptr() );
	kinematics::tree::AtomCOP c2prime_atom ( pose.atom_tree().atom( AtomID( rna_type.c2prime_atom_index(), i ) ).get_self_ptr() );

	if ( (c1prime_atom->parent()->id()).atomno() == first_base_atom_index( rsd ) ) {
		// There's a jump to this residue.
		//std::cout << "RESIDUE WITH JUMP CONNECTIVITY : " <<  i << std::endl;
		atoms_for_which_we_need_new_dofs.push_back( " C2'" );
		atoms_for_which_we_need_new_dofs.push_back( " C3'" );
		atoms_for_which_we_need_new_dofs.push_back( rsd.atom_name( rna_type.o4prime_atom_index() ) );
		atoms_for_which_we_need_new_dofs.push_back( " C4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C5'" );
		atoms_for_which_we_need_new_dofs.push_back( " O3'" );
	} else if ( (c2prime_atom->parent()->id()).atomno() ==  (o2prime_atom->id()).atomno() ) {
		atoms_for_which_we_need_new_dofs.push_back( " C1'" );
		atoms_for_which_we_need_new_dofs.push_back( " C3'" );
		atoms_for_which_we_need_new_dofs.push_back( rsd.atom_name( rna_type.o4prime_atom_index() ) );
		atoms_for_which_we_need_new_dofs.push_back( " C4'" );
		atoms_for_which_we_need_new_dofs.push_back( " C5'" );
		atoms_for_which_we_need_new_dofs.push_back( " O3'" );
	} else {
		atoms_for_which_we_need_new_dofs.push_back( " C1'" );
		atoms_for_which_we_need_new_dofs.push_back( " C2'" );
		atoms_for_which_we_need_new_dofs.push_back( rsd.atom_name( rna_type.o4prime_atom_index() ) );
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
) {
	using namespace id;

	/////////////////////////////////////////////
	// Save desired torsion values.
	utility::vector1< Real > start_torsions;
	for ( Size j = 1; j <= NUM_RNA_TORSIONS; j++ ) {
		id::TorsionID rna_torsion_id( i, id::BB, j );
		if ( j > NUM_RNA_MAINCHAIN_TORSIONS ) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
		start_torsions.push_back( reference_pose.torsion( rna_torsion_id ) );
	}

	/////////////////////////////////////////////
	//What DOFS do I need to get the ring atoms where I want them?
	utility::vector1< std::string > atoms_for_which_we_need_new_dofs;

	initialize_atoms_for_which_we_need_new_dofs( atoms_for_which_we_need_new_dofs,  pose, i );

	fix_sugar_coords( atoms_for_which_we_need_new_dofs, non_main_chain_sugar_coords, pose, reference_pose, i );

	/////////////////////////////////////////////
	// Reapply desired torsion values.
	for ( Size j = 1; j <= NUM_RNA_TORSIONS; j++ ) {
		id::TorsionID rna_torsion_id( i, id::BB, j );
		if ( j > NUM_RNA_MAINCHAIN_TORSIONS ) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
		pose.set_torsion( rna_torsion_id, start_torsions[ j ] );
	}
}

// AMW TODO: find invocations of this function, replace with RNA_IdealCoord
// Used in ConstrainToIdealMover
////////////////////////////////////////////////////////////////////
//FANG: All these sugar coord stuffs should be deprecated in favor of
//RNA_IdealCoord class? Are there any performance concern for using
//copy_dof_match_atom_name there?
void
apply_ideal_c2endo_sugar_coords(
	Pose & pose,
	Size const i
) {
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
	debug_assert( pucker_state <= 2 );

	static const RNA_IdealCoord ideal_coord;
	Real delta, nu1, nu2;

	PuckerState const curr_pucker = assign_pucker( pose, i );
	if ( skip_same_state && pucker_state == curr_pucker ) return;

	if ( pucker_state == core::chemical::rna::ANY_PUCKER ) pucker_state = curr_pucker;

	if ( idealize_coord ) {
		ideal_coord.apply_pucker(pose, i, pucker_state);
	} else {
		if ( pucker_state == NORTH ) {
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


//When a CUTPOINT_UPPER is added to 3' chain_break residue, the EXISTENCE of the CUTPOINT_UPPER atoms means that the alpha torsion which previously DOES NOT exist due to the chain_break now exist. The alpha value is automatically defined to the A-form value by Rosetta. However Rosetta does not automatically adjust the OP2 and OP1 atom position to account for this fact. So it is important that the OP2 and OP1 atoms position are correctly set to be consistent with A-form alpha torsion before the CUTPOINT_UPPER IS ADDED Parin Jan 2, 2009
void
position_cutpoint_phosphate_torsions( pose::Pose & current_pose,
	Size const five_prime_chainbreak,
	Size three_prime_chainbreak /* = 0 */ )
{

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::io::pdb;

	if ( three_prime_chainbreak == 0 ) three_prime_chainbreak = five_prime_chainbreak + 1;

	static const ResidueTypeSetCOP rsd_set = current_pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t );

	chemical::AA res_aa = chemical::na_rad; //aa_from_name( "RAD" );
	ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set->get_representative_type_aa( res_aa ) ) );

	current_pose.prepend_polymer_residue_before_seqpos( *new_rsd, three_prime_chainbreak, true );
	chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

	//Actually just by prepending the residue causes the alpha torsion to automatically be set to -64.0274,
	//so the manual setting below is actually not needed, May 24, 2010.. Parin S.
	//These are the initial value of virtual upper and lower cutpoint atom.
	//Actually only the alpha (id::BB, 1) is important here since it set the position of O3' (LOWER) atom which in turn determines  OP2 and OP1 atom
	current_pose.set_torsion( TorsionID( three_prime_chainbreak + 1, id::BB, 1 ), -64.027359 );
	current_pose.delete_polymer_residue( three_prime_chainbreak );
}


///////
bool
is_cutpoint_closed_by_atom_name(
	chemical::ResidueType const & rsd_1,
	chemical::ResidueType const & rsd_2,
	chemical::ResidueType const & rsd_3,
	chemical::ResidueType const & rsd_4,
	id::AtomID const & id1,
	id::AtomID const & id2,
	id::AtomID const & id3,
	id::AtomID const & id4 ){
	return(
		is_cutpoint_closed_atom( rsd_1, id1 ) ||
		is_cutpoint_closed_atom( rsd_2, id2 ) ||
		is_cutpoint_closed_atom( rsd_3, id3 ) ||
		is_cutpoint_closed_atom( rsd_4, id4 ) );
}

// Useful functions to torsional potential
bool
is_torsion_valid(
	pose::Pose const & pose,
	id::TorsionID const & torsion_id,
	bool verbose,
	bool skip_chainbreak_torsions
) {
	id::AtomID id1, id2, id3, id4;

	bool is_fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
	if ( verbose ) TR << "torsion_id: " << torsion_id << std::endl;
	if ( is_fail ) {
		if ( verbose ) TR << "fail to get torsion!, perhap this torsion is located at a chain_break " << std::endl;
		return false;
	}

	chemical::ResidueType const & rsd_1 = pose.residue_type( id1.rsd() );
	chemical::ResidueType const & rsd_2 = pose.residue_type( id2.rsd() );
	chemical::ResidueType const & rsd_3 = pose.residue_type( id3.rsd() );
	chemical::ResidueType const & rsd_4 = pose.residue_type( id4.rsd() );

	if ( rsd_1.has_variant_type( chemical::FIVEPRIME_CAP ) ) return false;
	if ( rsd_2.has_variant_type( chemical::FIVEPRIME_CAP ) ) return false;
	if ( rsd_3.has_variant_type( chemical::FIVEPRIME_CAP ) ) return false;
	if ( rsd_4.has_variant_type( chemical::FIVEPRIME_CAP ) ) return false;

	if ( !rsd_1.is_RNA() || !rsd_2.is_RNA() ||
			!rsd_3.is_RNA() || !rsd_4.is_RNA() ) return false;

	bool is_virtual_torsion = (
		rsd_1.is_virtual( id1.atomno() ) ||
		rsd_2.is_virtual( id2.atomno() ) ||
		rsd_3.is_virtual( id3.atomno() ) ||
		rsd_4.is_virtual( id4.atomno() ) );

	if ( is_virtual_torsion && verbose ) {
		print_torsion_info( pose, torsion_id );
	}

	// Check for cutpoint_closed
	// (Since these torsions will contain virtual atom(s),
	// but want to score these torsions
	bool const is_cutpoint_closed1 = is_cutpoint_closed_torsion( pose, torsion_id );

	// Oh, either this OR one of em has 5prime_cap
	// Note that having to make this exception is awful.... at least make it earlier.
	debug_assert( is_cutpoint_closed1 == is_cutpoint_closed_by_atom_name( rsd_1, rsd_2, rsd_3, rsd_4, id1, id2, id3, id4) ); // takes time

	// AMW: If this isn't a cutpoint in the FT, the below is not a problem (fingers crossed -- hack)
	if ( pose.fold_tree().is_cutpoint( id1.rsd() ) && is_cutpoint_closed1 && !is_virtual_torsion ) {
		print_torsion_info( pose, torsion_id );
		utility_exit_with_message( "is_cutpoint_closed1 == true && is_virtual_torsion == false!!" );
	}

	runtime_assert( id1.rsd() <= id2.rsd() );
	runtime_assert( id2.rsd() <= id3.rsd() );
	runtime_assert( id3.rsd() <= id4.rsd() );
	runtime_assert( ( id1.rsd() == id4.rsd() ) || ( id1.rsd() == ( id4.rsd() - 1 ) ) );

	bool const inter_residue_torsion = ( id1.rsd() != id4.rsd() );

	bool is_chain_break_torsion( false );

	if ( inter_residue_torsion ) {
		// Note that chain_break_torsion does not neccessarily have to be located
		// at a cutpoint_open. For example, in RNA might contain multiple strands,
		// but the user not have specified them as cutpoint_open
		// This happen frequently, for example when modeling single-stranded RNA
		// loop (PNAS December 20, 2011 vol. 108 no. 51 20573-20578).
		// Actually if chain_break is cutpoint_open,
		// pose.conformation().get_torsion_angle_atom_ids() should fail, which
		// leads to the EARLY RETURN FALSE statement at the beginning of
		// this function.
		bool const violate_max_O3_prime_to_P_bond_dist =
			pose::rna::is_rna_chainbreak( pose, id1.rsd() );

		// Note that cutpoint_closed_torsions are NOT considered as
		// chain_break_torsion since we want to score them EVEN when
		// skip_chainbreak_torsions_=true!
		// Necessary since for cutpoint_closed_torsions, the max
		// O3_prime_to_P_bond_dist might be violated during stages of the
		// Fragment Assembly and Stepwise Assembly where the chain is
		// not yet closed.
		is_chain_break_torsion = ( violate_max_O3_prime_to_P_bond_dist
			&& !is_cutpoint_closed1 );
	}

	// Warning before Jan 20, 2012, this used to be
	// "Size should_score_this_torsion= true;"
	bool should_score_this_torsion = true;

	if ( is_virtual_torsion && !is_cutpoint_closed1 ) should_score_this_torsion = false;
	if ( skip_chainbreak_torsions && is_chain_break_torsion ) should_score_this_torsion = false;

	if ( verbose ) {
		output_boolean( " should_score_torsion = ", should_score_this_torsion );
		output_boolean( " | is_cutpoint_closed = ", is_cutpoint_closed1 );
		output_boolean( " | is_virtual_torsion = ", is_virtual_torsion );
		output_boolean( " | skip_chainbreak_torsions = ", skip_chainbreak_torsions );
		output_boolean( " | is_chain_break_torsion   = ", is_chain_break_torsion ) ;
		TR << std::endl;
	}
	return should_score_this_torsion;
}

void
print_torsion_info(
	pose::Pose const & pose,
	id::TorsionID const & torsion_id
) {
	id::AtomID id1, id2, id3, id4;
	pose.conformation().get_torsion_angle_atom_ids(
		torsion_id, id1, id2, id3, id4 );

	TR << "torsion_id: " << torsion_id << std::endl;

	bool is_fail = pose.conformation().get_torsion_angle_atom_ids(
		torsion_id, id1, id2, id3, id4 );
	if ( is_fail ) {
		TR << "fail to get torsion!, perhap this torsion is located at a chain_break " << std::endl;
		return;
	}

	chemical::ResidueType const & rsd_1 = pose.residue_type( id1.rsd() );
	chemical::ResidueType const & rsd_2 = pose.residue_type( id2.rsd() );
	chemical::ResidueType const & rsd_3 = pose.residue_type( id3.rsd() );
	chemical::ResidueType const & rsd_4 = pose.residue_type( id4.rsd() );

	TR << " Torsion containing one or more virtual atom( s )" << std::endl;
	TR << "  torsion_id: " << torsion_id;
	TR << "  atom_id: " << id1 << " " << id2 << " " << id3 << " " << id4 << std::endl;
	TR << "  name: " << rsd_1.atom_name( id1.atomno() ) << " " <<
		rsd_2.atom_name( id2.atomno() ) << " " <<
		rsd_3.atom_name( id3.atomno() ) << " " <<
		rsd_4.atom_name( id4.atomno() ) << std::endl;
	TR << "  type: " << rsd_1.atom_type( id1.atomno() ).name() << " " <<
		rsd_2.atom_type( id2.atomno() ).name() << " " <<
		rsd_3.atom_type( id3.atomno() ).name() << " " <<
		rsd_4.atom_type( id4.atomno() ).name() << std::endl;
	TR << "\t\tatom_type_index: " << rsd_1.atom( id1.atomno() ).atom_type_index() <<
		" " << rsd_2.atom( id2.atomno() ).atom_type_index() << " " <<
		rsd_3.atom( id3.atomno() ).atom_type_index() << " " <<
		rsd_4.atom( id4.atomno() ).atom_type_index() << std::endl;
	TR << "\t\tatomic_charge: " << rsd_1.atom( id1.atomno() ).charge() <<
		" " << rsd_2.atom( id2.atomno() ).charge() << " " <<
		rsd_3.atom( id3.atomno() ).charge() << " " <<
		rsd_4.atom( id4.atomno() ).charge() << std::endl;
}

void
output_boolean( std::string const & tag, bool boolean ) {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	TR << tag;
	if ( boolean ) {
		TR << A( 4, "T" );
	} else {
		TR << A( 4, "F" );
	}
}

bool
is_cutpoint_closed_atom(
	core::chemical::ResidueType const & rsd,
	id::AtomID const & id
) {
	std::string const & atom_name = rsd.atom_name( id.atomno() );
	if ( atom_name == "OVU1" || atom_name == "OVL1" || atom_name == "OVL2" ) {
		return true;
	} else {
		return false;
	}
}

bool
is_cutpoint_closed_torsion(
	pose::Pose const & pose,
	id::TorsionID const & torsion_id
) {
	using namespace ObjexxFCL;
	using namespace core::scoring::methods;
	Size torsion_seq_num = torsion_id.rsd();
	Size lower_seq_num = 0;
	Size upper_seq_num = 0;

	if ( torsion_id.type() != id::BB ) return false;

	if ( torsion_id.torsion() == ALPHA ) { //COULD BE A UPPER RESIDUE OF A CHAIN_BREAK_CLOSE

		if ( !pose.residue_type( torsion_seq_num ).has_variant_type( chemical::CUTPOINT_UPPER ) ) return false;
		lower_seq_num = get_lower_cutpoint_partner_for_upper( pose, torsion_seq_num );
		upper_seq_num = torsion_seq_num;

	} else if ( torsion_id.torsion() == EPSILON || torsion_id.torsion() == ZETA ) {

		if ( !pose.residue_type( torsion_seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ) return false;
		lower_seq_num = torsion_seq_num;
		upper_seq_num = get_upper_cutpoint_partner_for_lower( pose, torsion_seq_num );

	} else {
		if ( torsion_id.torsion() != DELTA && torsion_id.torsion() != BETA && torsion_id.torsion() != GAMMA ) {
			utility_exit_with_message( "The torsion should be DELTA( lower ), BETA( upper ) or GAMMA( upper ) !!" );
		}
		return false;
	}

	if ( upper_seq_num == 0 ) return false;

	if ( lower_seq_num == 0 ) return false;

	if ( pose.residue_type( lower_seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ) {
		if ( !pose.residue_type( upper_seq_num ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			utility_exit_with_message( "seq_num " + string_of( lower_seq_num ) + " is a CUTPOINT_LOWER but seq_num " + string_of( upper_seq_num ) + " is not a cutpoint CUTPOINT_UPPER??" );
		}
		return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
apply_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num, bool const apply_check ){

	utility::vector1< Size > working_cutpoint_closed_list;
	working_cutpoint_closed_list.clear(); //empty list
	apply_virtual_rna_residue_variant_type( pose, seq_num, working_cutpoint_closed_list, apply_check );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
apply_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num, utility::vector1< core::Size > const & working_cutpoint_closed_list, bool const apply_check ){

	using namespace core::chemical;
	using namespace ObjexxFCL;

	runtime_assert( pose.size() >= seq_num );
	if ( pose.residue_type( seq_num ).has_variant_type( VIRTUAL_RNA_RESIDUE ) ) return;
	// AMW TODO: no reason these should actually be incompatible!
	//Basically the two variant type are not compatible, VIRTUAL_RNA_RESIDUE variant type currently does not virtualize the protonated H1 atom.
	if ( pose.residue( seq_num ).has_variant_type( PROTONATED_N1_ADENOSINE ) ) {
		TR << "Removing PROTONATED_N1_ADENOSINE variant_type from seq_num = " << seq_num <<
			" before adding VIRTUAL_RNA_RESIDUE variant_type since the two variant_types are not compatible!" <<
			std::endl;
		pose::remove_variant_type_from_pose_residue( pose, PROTONATED_N1_ADENOSINE, seq_num );
	}

	//OK PROTONATED_N1_ADENOSINE variant type should also be removed when adding VIRTUAL_RNA_RESIDUE_EXCLUDE_PHOSPHATE variant type or BULGE variant type.
	//However these two variant type are not currently used in standard SWA run (May 04, 2011)
	bool is_cutpoint_closed = false;
	if ( pose.residue_type( seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ) {
		runtime_assert( pose.residue_type( seq_num + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) );
		is_cutpoint_closed = true;
	}

	//Ok another possibility is that the CUTPOINT_LOWER and CUTPOINT_UPPER variant type had not been applied yet..so check the working_cutpoint_closed_list
	for ( auto const cutpoint_closed : working_cutpoint_closed_list ) {
		if ( seq_num == cutpoint_closed ) {
			is_cutpoint_closed = true;
			break;
		}
	}

	bool const is_cutpoint_open = ( pose.fold_tree().is_cutpoint( seq_num ) && !is_cutpoint_closed );
	if ( apply_check ) {
		if ( is_cutpoint_open ) {
			std::cerr << pose.annotated_sequence() << std::endl;
			std::cerr << pose.fold_tree();
			utility_exit_with_message( "Cannot apply VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". The residue is 5' of a OPEN cutpoint" );
		}
		if ( pose.size() == seq_num ) {
			utility_exit_with_message( "Cannot apply VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". pose.size() == seq_num" );
		}
	}

	pose::add_variant_type_to_pose_residue( pose, VIRTUAL_RNA_RESIDUE, seq_num );

	if ( ( pose.size() != seq_num ) &&  ( !is_cutpoint_open ) ) { //April 6, 2011
		pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, seq_num + 1 );
	}
}

//////////////////////////////////////////////////////////////////////////////////////
void
remove_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num ){

	using namespace core::chemical;
	using namespace ObjexxFCL;

	runtime_assert( seq_num < pose.size() );
	pose::remove_variant_type_from_pose_residue( pose, VIRTUAL_RNA_RESIDUE, seq_num );
	pose::remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, seq_num + 1 );
}

//////////////////////////////////////////////////////////////////////////////////////
bool
has_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num ){

	using namespace ObjexxFCL;

	if ( ! pose.residue_type( seq_num ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return false;

	if ( ( seq_num + 1 ) > pose.size() ) { //Check in range
		TR << "( seq_num + 1 ) = " << ( seq_num + 1 )  << std::endl;
		utility_exit_with_message( "( seq_num + 1 ) > pose.size()!" );
	}

	if ( ! pose.residue_type( seq_num + 1 ).has_variant_type( chemical::VIRTUAL_PHOSPHATE ) ) {
		TR << "Problem seq_num = " << seq_num << std::endl;
		utility_exit_with_message( "res ( " + string_of( seq_num ) +
			" ) has_variant_type VIRTUAL_RNA_RESIDUE but res seq_num + 1 ( " + string_of( seq_num + 1 ) +
			" )does not have variant_type VIRTUAL_PHOSPHATE" );
	}

	return true;
}

/////////////////////////////////////////////////////////////////////
bool
check_in_base_pair_list( pose::rna::BasePair const & base_pair /*from native*/,
	utility::vector1< core::pose::rna::BasePair > const & base_pair_list /*for decoy*/)
{
	using namespace pose::rna;
	return std::any_of( base_pair_list.begin(), base_pair_list.end(),
		[&]( BasePair const & bp2 ) {
			return bp2 == base_pair || bp2 == base_pair.flipped();
		} );
}

void
get_number_base_pairs( pose::Pose const & pose, Size & nwc, Size & nnwc, Size & bs ) {
	using namespace pose::rna;
	using namespace conformation;
	using namespace core::chemical::rna;

	utility::vector1< core::pose::rna::BasePair > base_pair_list;
	utility::vector1< bool > is_bulged;
	core::pose::rna::classify_base_pairs( pose, base_pair_list, is_bulged );

	Size N_WC( 0 ), N_NWC( 0 );

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		BasePair const base_pair = base_pair_list[ n ];

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		BaseEdge const k = base_pair.edge1();
		BaseEdge const m = base_pair.edge2();

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		// Virtual residues (perhaps added by build_full_model) don't count.
		if ( rsd_i.type().has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) continue;
		if ( rsd_j.type().has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) continue;

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK && base_pair.orientation() == ANTIPARALLEL )  &&
				core::chemical::rna::possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )  {
			N_WC++;
		} else {
			N_NWC++;
		}
	}
	nwc = N_WC;
	nnwc = N_NWC;
	bs = core::pose::rna::get_number_base_stacks( pose );
}

/////////////////////////////////////////////////////////////////////
void
add_number_base_pairs( pose::Pose & pose )
{
	Size N_WC, N_NWC, N_BS;
	get_number_base_pairs( pose, N_WC, N_NWC, N_BS );

	setPoseExtraScore( pose, "N_WC",  N_WC );
	setPoseExtraScore( pose, "N_NWC", N_NWC );
	setPoseExtraScore( pose, "N_BS",  N_BS );
}

/////////////////////////////////////////////////////////////////////
void
add_number_base_pairs( pose::Pose const & pose, io::silent::SilentStruct & s )
{
	Size N_WC, N_NWC, N_BS;
	get_number_base_pairs( pose, N_WC, N_NWC, N_BS );

	s.add_string_value( "N_WC",  ObjexxFCL::format::I( 9, N_WC ) );
	s.add_string_value( "N_NWC", ObjexxFCL::format::I( 9, N_NWC ) );
	s.add_string_value( "N_BS",  ObjexxFCL::format::I( 9, N_BS ) );
}

/////////////////////////////////////////////////////////////////////
void
get_number_native_base_pairs(pose::Pose & pose, pose::Pose const & native_pose,
	Size & pN_WC,
	Size & pN_NWC,
	Size & pN_BP,
	Size & pnatWC,
	Size & pnatNWC,
	Size & pnatBP,
	Real & pf_natWC,
	Real & pf_natNWC,
	Real & pf_natBP
) {
	// This is fundamentally comparative. When it was developed for FARFAR,
	// there was a guarantee that pose and native_pose would have the same
	// residues. Not so now.
	//
	// Rather than messing with FullModelInfo (which may not be set up for the
	// native, especially since FARFAR still needs this functionality) we'll use
	// the PDBInfo.

	using namespace chemical::rna;
	using namespace conformation;

	utility::vector1< core::pose::rna::BasePair > base_pair_list;
	utility::vector1< bool > is_bulged;
	core::pose::rna::classify_base_pairs( pose, base_pair_list, is_bulged );
	utility::vector1< core::pose::rna::BasePair > base_pair_list_native;
	utility::vector1< bool > is_bulged_native;
	core::pose::rna::classify_base_pairs( native_pose, base_pair_list_native, is_bulged_native );

	Size N_WC_NATIVE( 0 ), N_NWC_NATIVE( 0 );
	Size N_WC( 0 ), N_NWC( 0 );
	for ( BasePair const & native_base_pair : base_pair_list_native ) {
		Size const i = native_base_pair.res1();
		Size const j = native_base_pair.res2();

		BaseEdge const k = native_base_pair.edge1();
		BaseEdge const m = native_base_pair.edge2();

		// Does the pose-of-interest have these residues, too?
		// Map to new seqpos numbering, then create a base pair.
		Size i_decoy = 0;
		Size j_decoy = 0;
		if ( pose.size() == native_pose.size() ) {
			i_decoy = i;
			j_decoy = j;
		} else {
			// handle partially built stepwise decoys
			// chains might be unspecified for single-chain.
			for ( Size ii = 1; ii <= pose.size(); ++ii ) {
				if ( native_pose.pdb_info()->chain( i ) == pose.pdb_info()->chain( ii )
						|| pose.pdb_info()->chain( ii ) == ' '
						|| pose.pdb_info()->chain( ii ) == '^' ) {
					if ( native_pose.pdb_info()->number( i ) == pose.pdb_info()->number( ii ) ) {
						i_decoy = ii;
					}
				}
				if ( native_pose.pdb_info()->chain( j ) == pose.pdb_info()->chain( ii )
						|| pose.pdb_info()->chain( ii ) == ' '
						|| pose.pdb_info()->chain( ii ) == '^' ) {
					if ( native_pose.pdb_info()->number( j ) == pose.pdb_info()->number( ii ) ) {
						j_decoy = ii;
					}
				}
				if ( i_decoy != 0 && j_decoy != 0 ) break;
			}
		}
		// Do not proceed unless we have found both residues of interest.
		if ( !( i_decoy && j_decoy ) ) continue;

		Residue const & rsd_i( pose.residue( i_decoy ) );
		Residue const & rsd_j( pose.residue( j_decoy ) );

		BasePair decoy_base_pair( native_base_pair );
		decoy_base_pair.set_res1( i_decoy );
		decoy_base_pair.set_res2( j_decoy );

		//std::cout << " NATIVE BASE PAIR " << i << " " << j << " " << k << " " << m << " " << it->first << std::endl;

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK && decoy_base_pair.orientation() == ANTIPARALLEL )  &&
				possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )  {
			N_WC_NATIVE++;
			if ( check_in_base_pair_list( decoy_base_pair /*from native*/, base_pair_list /*for decoy*/) ) N_WC++;
		} else {
			N_NWC_NATIVE++;
			if ( check_in_base_pair_list( decoy_base_pair /*from native*/, base_pair_list /*for decoy*/) ) {
				N_NWC++;
			} else {
				std::cout << "Missing native base pair " << pose.residue_type( i ).name1() << i << " " << pose.residue_type(j).name1() << j << "  " << get_edge_from_num( k ) << " " << get_edge_from_num( m ) << " " << std::endl;
			}
		}
	}

	Real f_natWC( 0.0 ), f_natNWC( 0.0 ), f_natBP( 0.0 );
	if ( N_WC_NATIVE > 0 ) f_natWC = ( N_WC / (1.0 * N_WC_NATIVE) );
	if ( N_NWC_NATIVE > 0 ) f_natNWC = ( N_NWC / (1.0 * N_NWC_NATIVE) );
	if ( (N_WC_NATIVE + N_NWC_NATIVE) > 0 ) f_natBP = ( (N_WC+N_NWC) / (1.0 * (N_WC_NATIVE + N_NWC_NATIVE) ));

	// Adding these here helps control for false positives found in the other function
	// (Which might count as 'recovered' BPs not found in the native model.
	pN_WC = N_WC;
	pN_NWC = N_NWC;
	pN_BP = N_WC + N_NWC;
	pnatWC = N_WC_NATIVE;
	pnatNWC = N_NWC_NATIVE;
	pnatBP = N_WC_NATIVE + N_NWC_NATIVE;
	pf_natWC = f_natWC;
	pf_natNWC = f_natNWC;
	pf_natBP = f_natBP;
}

/////////////////////////////////////////////////////////////////////
void
add_number_native_base_pairs(pose::Pose & pose, pose::Pose const & native_pose, io::silent::SilentStruct & s )
{
	Size N_WC, N_NWC, N_BP, natWC, natNWC, natBP;
	Real f_natWC, f_natNWC, f_natBP;

	get_number_native_base_pairs( pose, native_pose, N_WC, N_NWC, N_BP, natWC, natNWC, natBP, f_natWC, f_natNWC, f_natBP );

	// Adding these here helps control for false positives found in the other function
	// (Which might count as 'recovered' BPs not found in the native model.
	s.add_string_value( "N_WC",  ObjexxFCL::format::I( 9, N_WC) );
	s.add_string_value( "N_NWC", ObjexxFCL::format::I( 9, N_NWC ) );
	s.add_string_value( "N_BP",  ObjexxFCL::format::I( 9, N_BP ) );

	// I would personally like to see the native number esp. of NWC because
	// otherwise the fraction could be deflated by false 0s.
	s.add_string_value( "natWC" , ObjexxFCL::format::I( 9, natWC ) );
	s.add_string_value( "natNWC", ObjexxFCL::format::I( 9, natNWC ) );
	s.add_string_value( "natBP" , ObjexxFCL::format::I( 9, natBP ) );

	s.add_energy( "f_natWC" , f_natWC );
	s.add_energy( "f_natNWC", f_natNWC );
	s.add_energy( "f_natBP" , f_natBP );
}


/////////////////////////////////////////////////////////////////////
void
add_number_native_base_pairs(pose::Pose & pose, pose::Pose const & native_pose )
{
	Size N_WC, N_NWC, N_BP, natWC, natNWC, natBP;
	Real f_natWC, f_natNWC, f_natBP;

	get_number_native_base_pairs( pose, native_pose, N_WC, N_NWC, N_BP, natWC, natNWC, natBP, f_natWC, f_natNWC, f_natBP );

	// Adding these here helps control for false positives found in the other function
	// (Which might count as 'recovered' BPs not found in the native model.
	setPoseExtraScore( pose, "N_WC",  N_WC );
	setPoseExtraScore( pose, "N_NWC", N_NWC );
	setPoseExtraScore( pose, "N_BP",  N_BP );

	// I would personally like to see the native number esp. of NWC because
	// otherwise the fraction could be deflated by false 0s.
	setPoseExtraScore( pose, "natWC" , natWC );
	setPoseExtraScore( pose, "natNWC", natNWC );
	setPoseExtraScore( pose, "natBP" , natBP );

	setPoseExtraScore( pose, "f_natWC" , f_natWC );
	setPoseExtraScore( pose, "f_natNWC", f_natNWC );
	setPoseExtraScore( pose, "f_natBP" , f_natBP );
}

//////////////////////////////////////////////////////////////////////////////////////
void
apply_Aform_torsions( pose::Pose & pose, Size const n ){
	using namespace core::id;
	// following does not really matter, since sugar & phosphate will be erased anyway.
	core::chemical::rna::RNA_FittedTorsionInfo const torsion_info_;
	pose.set_torsion( TorsionID( n, BB, 1), torsion_info_.alpha_aform() );
	pose.set_torsion( TorsionID( n, BB, 2), torsion_info_.beta_aform() );
	pose.set_torsion( TorsionID( n, BB, 3), torsion_info_.gamma_aform() );
	pose.set_torsion( TorsionID( n, BB, 4), torsion_info_.delta_north() );
	pose.set_torsion( TorsionID( n, BB, 5), torsion_info_.epsilon_aform() );
	pose.set_torsion( TorsionID( n, BB, 6), torsion_info_.zeta_aform() );
	apply_pucker(pose, n, NORTH, false /*skip_same_state*/, true /*use_phenix_geo_*/);
}


//////////////////////////////////////////////////////////////////////////////////////
ChiState
get_residue_base_state( core::pose::Pose const & pose, Size const seq_num ){
	return core::chemical::rna::get_residue_base_state( pose.residue( seq_num ) );
}

//////////////////////////////////////////////////////////////////////////////////////
PuckerState
get_residue_pucker_state( core::pose::Pose const & pose, Size const seq_num ){
	return core::chemical::rna::get_residue_pucker_state( pose.residue( seq_num ) );
}

/////////////////////////////////////////////////////////////////////////////////////
/// Could this be made more general? Figured out through knowledge of where side
///  chain atoms connect to polymeric backbone?
utility::vector1< std::pair< core::id::TorsionID, core::Real > >
get_suite_torsion_info( core::pose::Pose const & pose, Size const i )
{
	using namespace utility;
	using namespace core::id;
	vector1< TorsionID > torsion_ids;
	if ( pose.residue_type( i ).is_NA() ) {
		torsion_ids.push_back( TorsionID( i  , BB, 5 ) ); // epsilon
		torsion_ids.push_back( TorsionID( i  , BB, 6 ) ); // zeta
		torsion_ids.push_back( TorsionID( i+1, BB, 1 ) ); // alpha
		torsion_ids.push_back( TorsionID( i+1, BB, 2 ) ); // beta
		torsion_ids.push_back( TorsionID( i+1, BB, 3 ) ); // gamma
	} else if ( pose.residue_type( i ).is_protein() ) {
		torsion_ids.push_back( TorsionID( i  , BB, 2 ) ); // psi
		torsion_ids.push_back( TorsionID( i  , BB, 3 ) ); // omega
		torsion_ids.push_back( TorsionID( i+1, BB, 1 ) ); // phi
	}
	vector1< std::pair< TorsionID, Real > > suite_torsion_info;
	for ( Size n = 1; n <= torsion_ids.size(); n++ ) {
		suite_torsion_info.push_back( std::make_pair( torsion_ids[ n ], pose.torsion( torsion_ids[ n ] ) ) );
	}
	return suite_torsion_info;
}

/////////////////////////////////////////////////////////////////////////////////////
void
apply_suite_torsion_info( core::pose::Pose & pose,
	utility::vector1< std::pair< core::id::TorsionID, core::Real > > const & suite_torsion_info ) {
	for ( auto const & elem : suite_torsion_info ) {
		pose.set_torsion( elem.first, elem.second );
	}
}

////////////////////////////////////////////////////////////////////////////////////
Real
get_op2_op1_sign( pose::Pose const & pose ) {

	Real sign= 0;
	bool found_valid_sign=false;

	for ( Size i = 2; i <= pose.size(); i++ ) {

		conformation::Residue const & rsd( pose.residue( i )  );
		if ( !rsd.is_RNA() ) continue;

		sign = dot( rsd.xyz( " O5'" ) - rsd.xyz( " P  " ), cross( rsd.xyz( " OP1" ) - rsd.xyz( " P  " ), rsd.xyz( " OP2" ) - rsd.xyz( " P  " ) ) );
		found_valid_sign = true;
		break;
	}

	if ( found_valid_sign==false ) utility_exit_with_message("found_valid_sign==false");

	return sign;
}

////////////////////////////////////////////////////////////////////////////////////
//This version used to be called get_op2_op1_sign_parin()
Real
get_op2_op1_sign( pose::Pose const & pose , Size res_num) {

	if ( res_num > pose.size() ) utility_exit_with_message("res_num > pose.size()");

	conformation::Residue const & rsd( pose.residue(res_num)  );

	//SML PHENIX conference cleanup
	if ( basic::options::option[basic::options::OptionKeys::rna::erraser::rna_prot_erraser].value() ) {
		if ( !rsd.is_RNA() ) return 0.0;
	} else {
		if ( rsd.is_RNA()==false ) utility_exit_with_message("rsd.is_RNA()==false!");
	}

	Real const sign = dot( rsd.xyz( " O5'" ) - rsd.xyz( " P  " ), cross( rsd.xyz( " OP1" ) - rsd.xyz( " P  " ), rsd.xyz( " OP2" ) - rsd.xyz( " P  " ) ) );

	return sign;
}

////////////////////////////////////////////////////////////////
void
make_phosphate_nomenclature_matches_mini( pose::Pose & pose)
{
	using namespace ObjexxFCL;
	pose::Pose mini_pose;
	make_pose_from_sequence( mini_pose, "aa", pose.residue_type_set_for_pose() );
	Real const sign2 = get_op2_op1_sign( mini_pose);

	for ( Size res_num=1; res_num<=pose.size(); res_num++ ) {

		if ( !pose.residue( res_num ).is_RNA() ) continue;
		Real sign1 = get_op2_op1_sign( pose,  res_num);

		if ( sign1 * sign2 < 0 ) {
			//std::cout << " Flipping OP2 <--> OP1 " << "res_num " << res_num << " | sign1: " << sign1 << " | sign2: " << sign2 << std::endl;
			conformation::Residue const & rsd( pose.residue(res_num) );

			if ( rsd.is_RNA()==false ) { //Consistency check!
				std::cout << "residue # " << res_num << " should be a RNA nucleotide!" << std::endl;
				utility_exit_with_message("residue # " + string_of(res_num)+ " should be a RNA nucleotide!");
			};

			Vector const temp1 = rsd.xyz( " OP2" );
			Vector const temp2 = rsd.xyz( " OP1" );
			pose.set_xyz( id::AtomID( rsd.atom_index( " OP2" ), res_num ), temp2 );
			pose.set_xyz( id::AtomID( rsd.atom_index( " OP1" ), res_num ), temp1 );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_virtual_O2Prime_hydrogen( core::pose::Pose & pose ){
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		if ( !pose.residue_type( i ).is_RNA() ) continue;
		pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_O2PRIME_HYDROGEN, i );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Took these functions out of the class so that they are accessible in the VDWGridEnergy for scoring
// without actually needing to construct an RNA_VDW_BinChecker object
Atom_Bin
get_atom_bin( numeric::xyzVector< core::Real > const & atom_pos, numeric::xyzVector< core::Real > const & ref_xyz,
	core::Real const atom_bin_size, int const bin_offset ) {

	numeric::xyzVector< core::Real > const atom_pos_ref_frame = atom_pos - ref_xyz;

	Atom_Bin atom_bin;
	atom_bin.x = int( atom_pos_ref_frame[0]/atom_bin_size );
	atom_bin.y = int( atom_pos_ref_frame[1]/atom_bin_size );
	atom_bin.z = int( atom_pos_ref_frame[2]/atom_bin_size );


	if ( atom_pos_ref_frame[0] < 0 ) atom_bin.x--;
	if ( atom_pos_ref_frame[1] < 0 ) atom_bin.y--;
	if ( atom_pos_ref_frame[2] < 0 ) atom_bin.z--;

	//////////////////////////////////////////////////////////
	atom_bin.x += bin_offset; //Want min bin to be at one.
	atom_bin.y += bin_offset; //Want min bin to be at one.
	atom_bin.z += bin_offset; //Want min bin to be at one.


	//////////////////////////////////////////////////////////
	return atom_bin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_atom_bin_in_range( Atom_Bin const & atom_pos_bin, int const bin_max ) {

	if ( atom_pos_bin.x < 1 || ( atom_pos_bin.x > ( bin_max*2 ) ) ||
			atom_pos_bin.y < 1 || ( atom_pos_bin.y > ( bin_max*2 ) ) ||
			atom_pos_bin.z < 1 || ( atom_pos_bin.z > ( bin_max*2 ) ) ) {
		return false;
	} else {
		return true;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< std::string >
tokenize( std::string const & str, std::string const & delimiters ){
	using namespace std;

	utility::vector1< std::string > tokens;

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of( delimiters, lastPos );

	while ( string::npos != pos || string::npos != lastPos ) {
		// Found a token, add it to the vector.
		tokens.push_back( str.substr( lastPos, pos - lastPos ) );
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of( delimiters, pos );
		// Find next "non-delimiter"
		pos = str.find_first_of( delimiters, lastPos );
	}
	return tokens;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
core::Size
string_to_int( std::string const & input_string ){

	Size int_of_string; //misnomer
	std::stringstream ss ( std::stringstream::in | std::stringstream::out );

	ss << input_string;

	if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

	ss >> int_of_string;

	if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> int_of_string | string ( " + input_string + " )" );

	return int_of_string;
}

std::string
remove_bracketed( std::string const & sequence ) {
	std::string return_val = sequence;

	std::string::size_type i = return_val.find( '[' );
	while ( i != std::string::npos ) {

		std::string::size_type j = return_val.find( ']' );
		if ( j == std::string::npos ) utility_exit_with_message( "String with imbalanced brackets passed to remove_bracketed! (" + sequence + ")" );

		return_val.erase( i, j - i + 1 );

		i = return_val.find( '[' );
	}

	return return_val;
}

void
remove_and_store_bracketed( std::string const & working_sequence, std::string & working_sequence_clean, std::map< Size, std::string > & special_res ) {
	working_sequence_clean = working_sequence;

	std::string::size_type i = working_sequence_clean.find( '[' );
	while ( i != std::string::npos ) {

		std::string::size_type j = working_sequence_clean.find( ']' );
		if ( j == std::string::npos ) utility_exit_with_message( "String with imbalanced brackets passed to remove_bracketed! (" + working_sequence + ")" );

		std::string temp = working_sequence_clean.substr( i + 1, j - i - 1 );
		special_res[ i - 1 ] = temp;
		working_sequence_clean.erase( i, j - i + 1 );

		i = working_sequence_clean.find( '[' );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
add_chi_constraints( pose::Pose & pose,
	core::scoring::func::FuncOP chi_potential_restraint,
	utility::vector1< Size > const & rna_chi_res )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::pose::full_model_info;

	utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
	for ( Size const n : rna_chi_res ) {
		//TR << "Thinking about " << n << std::endl;
		if ( ! res_list.has_value( n ) ) continue;

		//TR << "Trying " << n << std::endl;

		Size const j = res_list.index( n );
		runtime_assert( pose.residue_type( j ).is_RNA() );

		AtomID id1,id2,id3,id4;
		pose.conformation().get_torsion_angle_atom_ids( TorsionID( j, core::id::CHI, 1 ), id1, id2, id3, id4 );
		pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( id1, id2, id3, id4, chi_potential_restraint, rna_torsion ) ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
add_syn_chi_constraints( core::pose::Pose & pose ) {
	using namespace core::scoring::func;
	using namespace core::pose::full_model_info;
	if ( !full_model_info_defined( pose ) ) return;
	// define func
	utility::vector1< Real > energies( 36, 20.0 ); // spaced at 5, 10, ... 355 degrees. Set penalty to be +20.0
	// zero out function at permissible vals. (syn should be -120.0 to 0.0 --> 240.0 to 360.0)
	for ( Size q = 24; q <= 36; q++ ) energies[ q ] = 0.0;
	FuncOP chi_potential_syn_restraint( new CircularSplineFunc( 1.0, energies, true /*convert_to_degrees*/ ) );
	add_chi_constraints( pose, chi_potential_syn_restraint, const_full_model_info( pose ).rna_syn_chi_res() );
}

////////////////////////////////////////////////////////////////////////////////////////
void
add_anti_chi_constraints( core::pose::Pose & pose ) {
	using namespace core::scoring::func;
	using namespace core::pose::full_model_info;
	if ( !full_model_info_defined( pose ) ) return;
	// define func
	utility::vector1< Real > energies( 36, 20.0 ); // spaced at 5, 10, ... 355 degrees. Set penalty to be +20.0
	// zero out function at permissible vals. (anti should be 0.0 to 120.0)
	for ( Size q = 1; q <= 12; q++ ) energies[ q ] = 0.0;
	FuncOP chi_potential_anti_restraint( new CircularSplineFunc( 1.0, energies, true /*convert_to_degrees*/ ) );
	add_chi_constraints( pose, chi_potential_anti_restraint, const_full_model_info( pose ).rna_anti_chi_res() );
}

////////////////////////////////////////////////////////////////////////////////////////
void
add_syn_anti_chi_constraints( core::pose::Pose & pose ) {
	using namespace core::pose::full_model_info;
	if ( !full_model_info_defined( pose ) ) return;
	add_syn_chi_constraints( pose );
	add_anti_chi_constraints( pose );
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::id::TorsionID >
get_suite_torsion_ids( Size const i )
{
	using namespace core::id;
	using namespace utility::tools;
	return make_vector1(
		TorsionID( i    , BB, EPSILON ),
		TorsionID( i    , BB, ZETA ),
		TorsionID( i + 1, BB, ALPHA ),
		TorsionID( i + 1, BB, BETA ),
		TorsionID( i + 1, BB, GAMMA ) );
}

void
get_stub_stub( conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	kinematics::Stub & stub1,
	kinematics::Stub & stub2,
	StubStubType const & stub_stub_type )
{
	using namespace core::kinematics;
	typedef numeric::xyzMatrix< Real > Matrix;
	switch ( stub_stub_type ) {
	case O3P_TO_O5P :
		// takeoff
		stub1 = Stub( rsd1.xyz( " O3'") /* center */,
			rsd1.xyz( " O3'") /* a */,
			rsd1.xyz( " C3'") /* b  [b->a defines x] */,
			rsd1.xyz( " C4'") /* c  [c->b defines y] */ );
		stub1.M = Matrix::cols( stub1.M.col_y(), stub1.M.col_z(), stub1.M.col_x() ); // Prefer to have C3'->O3' (takeoff vector) along z

		// landing
		stub2 = Stub( rsd2.xyz( " O5'") /* center */,
			rsd2.xyz( " C5'") /* a */,
			rsd2.xyz( " O5'") /* b  [b->a defines x] */,
			rsd2.xyz( " C4'") /* c  [c->b defines y] */ );
		stub2.M = Matrix::cols( stub2.M.col_y(), stub2.M.col_z(), stub2.M.col_x() ); // Prefer to have O5'->C5' (landing vector) along z
		return;
	case CHAINBREAK :
		if ( !rsd1.has_variant_type( chemical::CUTPOINT_LOWER ) ) return;
		if ( !rsd2.has_variant_type( chemical::CUTPOINT_UPPER ) ) return;
		// takeoff
		stub1 = Stub(
			rsd1.xyz( "OVL1" ) /* center */,
			rsd1.xyz( "OVL2" ) /* a */,
			rsd1.xyz( "OVL1" ) /* b  [b->a defines x] */,
			rsd1.xyz( rsd1.mainchain_atoms()[ rsd1.mainchain_atoms().size() ]  ) /* c  [c->b defines y] */ );
		stub1.M = Matrix::cols( stub1.M.col_y(), stub1.M.col_z(), stub1.M.col_x() ); // Prefer to have C3'->O3' (takeoff vector) along z

		// landing
		stub2 = Stub(
			rsd2.xyz( rsd2.mainchain_atoms()[ 1 ] ) /* center */,
			rsd2.xyz( rsd2.mainchain_atoms()[ 2 ] ) /* a */,
			rsd2.xyz( rsd2.mainchain_atoms()[ 1 ] ) /* b  [b->a defines x] */,
			rsd2.xyz( "OVU1" ) /* c  [c->b defines y] */ );
		stub2.M = Matrix::cols( stub2.M.col_y(), stub2.M.col_z(), stub2.M.col_x() ); // Prefer to have O5'->C5' (landing vector) along z
		return;
	case BASE_CENTROID :
		stub1 = get_rna_base_coordinate_system_stub( rsd1 );
		stub2 = get_rna_base_coordinate_system_stub( rsd2 );
		return;
	default :
		utility_exit_with_message( "Unrecognized stub_stub_type" );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
output_base_pairs( std::ostream & out, core::pose::rna::RNA_BasePairList const & base_pair_list, pose::Pose const & pose  )
{
	using namespace core::pose::rna;
	using namespace core::chemical::rna;
	RNA_BasePairList base_pair_list_sort = base_pair_list;
	std::sort( base_pair_list_sort.end(), base_pair_list_sort.end() );
	for ( auto const & base_pair : base_pair_list_sort ) {
		out <<
			pose.pdb_info()->tag( base_pair.res1() ) << " " <<
			pose.pdb_info()->tag( base_pair.res2() ) << " " <<
			get_edge_from_num( base_pair.edge1() ) << " " <<
			get_edge_from_num( base_pair.edge2() ) << " " <<
			get_LW_orientation_from_num( base_pair.LW_orientation() ) << " " <<
			std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
output_base_stacks( std::ostream & out, core::pose::rna::RNA_BaseStackList const & base_stack_list, pose::Pose const & pose  )
{
	using namespace core::pose::rna;
	using namespace core::chemical::rna;
	RNA_BaseStackList base_stack_list_sort = base_stack_list;
	std::sort( base_stack_list_sort.end(), base_stack_list_sort.end() );
	for ( auto const & base_stack : base_stack_list_sort ) {
		out <<
			pose.pdb_info()->tag( base_stack.res1() ) << " " <<
			pose.pdb_info()->tag( base_stack.res2() ) << " " <<
			get_side_from_num( base_stack.which_side() ) << " " <<
			get_orientation_from_num( base_stack.orientation() ) << " " <<
			std::endl;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////
// @details
// used by rna_motif app.
//
// Really should deprecate this in favor of figure_out_stems() function in RNA_SecStruct()!
//  I had forgotten about that when I wrote this function in late 2017, oops. -- rhiju.
//
/////////////////////////////////////////////////////////////////////////////////////////////
void
output_stems( std::ostream & out, RNA_Motifs const & rna_motifs, pose::Pose const & pose )
{
	using utility::vector1;

	// Order so that stacked pairs will be contiguous.
	// Otherwise risking that a stack seeded at one end will be 'different' from
	// stack at seeded at other end, even if they're part of one long stack.

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Should actually separate out the function that gets STACKED PAIRS from pose  -- don't need to run all of RNA_Motif...
	// Could also figure out stacked_pairs just from an ordered list of WC base pairs...
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector1< RNA_Motif > stacked_pairs = rna_motifs.get_motifs( WC_STACKED_PAIR );
	std::sort( stacked_pairs.begin(), stacked_pairs.end() );

	// group stacked_pairs into long stems
	vector1< Size > in_stem( pose.size(), 0 );
	vector1< std::set< Size > > stem_residues;
	for ( auto const & stacked_pair : stacked_pairs ) {
		Size stem = 0;
		if ( stacked_pair[1] > stacked_pair[3] ) continue; // would be better to force uniqueness above. anyway.
		for ( auto const & m : stacked_pair.residues() ) {
			if ( in_stem[ m ] > 0 ) {
				if ( stem > 0 ) assert( stem == in_stem[ m ] );
				else stem = in_stem[ m ];
			}
		}
		if ( stem == 0 ) {
			// no residue in this stacked pair is shared with prior stem. So make a new stem.
			stem_residues.push_back( std::set< Size >() );
			stem = stem_residues.size();
		}
		for ( auto const & m : stacked_pair.residues() ) {
			in_stem[ m ] = stem;
			stem_residues[ stem ].insert( m );
		}
	}

	// Each long stem should have two strands, recognizable as two consecutive chunks of residues.
	vector1< vector1< vector1< Size > > > strand_sets;
	for ( auto const & stem_reslist : stem_residues ) {
		vector1< vector1< Size > > strand_set;
		vector1< Size > strand;
		Size prev_m = 0;
		for ( auto const m : stem_reslist ) {
			if ( strand.size() > 0 && ( m != prev_m + 1 ||
					pose.pdb_info()->chain(m) != pose.pdb_info()->chain(prev_m ) ||
					pose.pdb_info()->segmentID(m) != pose.pdb_info()->segmentID(prev_m ) ||
					pose.pdb_info()->number(m) != pose.pdb_info()->number(prev_m)+1  ||
					pose.fold_tree().is_cutpoint( prev_m ) ) ) {
				strand_set.push_back( strand );
				strand.clear();
			}
			strand.push_back( m );
			prev_m = m;
		}
		strand_set.push_back( strand );
		if ( strand_set.size() != 2 ) {
			std::cout << "Problem with strand_set: " << strand_set << std::endl;
			utility_exit_with_message( "strand_set.size() != 2" );
		}
		strand_sets.push_back( strand_set );
	}

	// output stacks -- could be used for secondary structure determination.
	for ( auto const & strand_set : strand_sets ) {
		vector1< std::string > tags;
		for ( auto const & strand : strand_set ) {
			vector1< int > res;
			vector1< char > chain;
			vector1< std::string > segid;
			for ( auto const & m : strand ) {
				res.push_back( pose.pdb_info()->number( m ) );
				chain.push_back( pose.pdb_info()->chain( m ) );
				segid.push_back( pose.pdb_info()->segmentID( m ) );
			}
			tags.push_back( make_tag_with_dashes( res, chain, segid ) );
		}
		out << tags[1] << " " << tags[2] << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< std::pair< char, std::string > >
figure_out_rna_chains(
	pose::Pose const & pose,
	utility::vector1< std::string > const & chains /* = utility::vector1< std::string >*/
)
{
	utility::vector1< std::pair< char, std::string > > all_rna_chain_segids;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue_type( n ).is_RNA() ) {
			ChainSegID chain_segid( std::make_pair( pose.pdb_info()->chain( n ), strip( pose.pdb_info()->segmentID( n ) ) ) );
			if ( !all_rna_chain_segids.has_value( chain_segid ) ) all_rna_chain_segids.push_back( chain_segid );
		}
	}

	if ( chains.size() == 0 ) return all_rna_chain_segids;

	utility::vector1< ChainSegID > chain_segids;
	for ( auto const & chain_string : chains ) {
		char const chain( chain_string[ 0 ] );
		std::string segid;
		if ( chain_string.size() > 1 ) {
			if ( chain_string[1] == ':' ) segid = chain_string.substr( 2 );
			else segid = chain_string.substr(1);
		}
		ChainSegID chain_segid( chain, segid );
		if ( ! all_rna_chain_segids.has_value( chain_segid ) ) {
			std::cout << "all_rna_chain_segids: " << all_rna_chain_segids << std::endl;
			utility_exit_with_message( "Did not find desired chain "+std::string(1,chain)+" segid "+segid+" in all_rna_chain_segids above." );
		}
		chain_segids.push_back( std::make_pair( chain, segid ) );
	}

	return chain_segids;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
pose::Pose
extract_rna_chains( pose::Pose const & full_pose, utility::vector1< ChainSegID > const & chain_segids )
{
	using namespace core::pose::rna;
	utility::vector1< Size > slice_res;
	for ( Size n = 1; n <= full_pose.size(); n++ ) {
		ChainSegID chain_segid( std::make_pair( full_pose.pdb_info()->chain( n ), strip( full_pose.pdb_info()->segmentID( n ) ) ) );
		if ( chain_segids.has_value( chain_segid ) ) {
			if ( !full_pose.residue_type( n ).is_RNA() ) continue;
			slice_res.push_back( n );
		}
	}
	pose::Pose pose;
	pdbslice( pose, full_pose, slice_res );
	remove_upper_lower_variants_from_RNA( pose );
	figure_out_reasonable_rna_fold_tree( pose, true /*put cutpoints at chainbreaks*/  );
	return pose;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_map( std::map< ChainSegID, std::set< core::Size > > & ligand_map,
	std::map< ChainSegID, std::string > & ligand_tag,
	Size const & i /* partner of ligand in RNA chain*/,
	ChainSegID const & chain_segid_j /*ligand*/,
	core::conformation::Residue const & rsd_j /*ligand*/)
{
	ligand_map[ chain_segid_j ].insert( i );

	if ( rsd_j.type().is_polymer() ) {
		if ( rsd_j.type().is_RNA() ) {
			ligand_tag[ chain_segid_j ] = "RNA";
		} else if ( rsd_j.type().is_protein() ) {
			ligand_tag[ chain_segid_j ] = "protein";
		} else if ( rsd_j.type().is_DNA() ) {
			ligand_tag[ chain_segid_j ] = "DNA";
		} else {
			ligand_tag[ chain_segid_j ] = "other";
		}
	} else {
		ligand_tag[ chain_segid_j ] = rsd_j.name();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details TODO: We should recognize non-RNA bits that are supplied in the PDB with the same chains as the RNA polymers.
///            Happens, e.g., in the glycine riboswitch 3P49 where the glycines are in chain A along with the RNA.
void
output_ligands( std::ostream & out, pose::Pose const & pose,  utility::vector1< ChainSegID > const & chain_segids )
{
	using namespace core::scoring;
	using namespace core::conformation;
	using namespace utility::graph;

	TenANeighborGraph const & tenA(pose.energies().tenA_neighbor_graph());

	std::map< ChainSegID, std::set< core::Size > > ligand_map;
	std::map< ChainSegID, std::string > ligand_tag;


	// For each neighboring residue to the donor
	for ( EdgeListConstIterator ni = tenA.const_edge_list_begin(),
			ni_end = tenA.const_edge_list_end();
			ni != ni_end; ++ni ) {
		Size const & i = (*ni)->get_first_node_ind();
		Size const & j = (*ni)->get_second_node_ind();
		ChainSegID const & chain_segid_i = std::make_pair( pose.pdb_info()->chain( i ), strip(pose.pdb_info()->segmentID( i ) ) );
		ChainSegID const & chain_segid_j = std::make_pair( pose.pdb_info()->chain( j ), strip(pose.pdb_info()->segmentID( j ) ) );

		// will be covered by RNA base pair/motif check -- so not considered a 'ligand' interaction.
		if ( chain_segids.has_value( chain_segid_i ) && chain_segids.has_value( chain_segid_j ) ) continue;

		static Distance const CONTACT_DIST_CUTOFF( 4.0 );

		// now look for a direct contact between these residues
		// might be faster to look up atom-level neighbor graph
		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );
		bool found_contact( false );
		for ( Size m = 1; m <= rsd_i.nheavyatoms(); m++ ) {
			for ( Size n = 1; n <= rsd_j.nheavyatoms(); n++ ) {
				if ( ( rsd_i.xyz( m ) - rsd_j.xyz( n ) ).length() <  CONTACT_DIST_CUTOFF ) {
					found_contact = true;
					break;
				}
			}
		}
		if ( !found_contact ) continue;

		if ( chain_segids.has_value( chain_segid_i ) && !chain_segids.has_value( chain_segid_j ) ) {
			update_map( ligand_map, ligand_tag, i, chain_segid_j, rsd_j );
		} else if ( !chain_segids.has_value( chain_segid_i ) && chain_segids.has_value( chain_segid_j ) ) {
			update_map( ligand_map, ligand_tag, j, chain_segid_i, rsd_i );
		} else {
			//TR << chain_segid_i << " " << chain_segid_j << std::endl;
		}

	}

	for ( auto const & it : ligand_map ) {
		ChainSegID ligand_chain_segid( it.first );
		out << ligand_chain_segid.first;
		if ( ligand_chain_segid.second.size() > 0 ) out << ":" << ligand_chain_segid.second;
		out << "     ";

		out << ObjexxFCL::right_string_of( ligand_tag[ ligand_chain_segid ], 7);
		out << "     ";

		utility::vector1< int > resnum;
		utility::vector1< char > chain;
		utility::vector1< std::string > segid;
		for ( auto const & pos: it.second ) {
			resnum.push_back( pose.pdb_info()->number( pos ) );
			chain.push_back( pose.pdb_info()->chain( pos ) );
			segid.push_back( pose.pdb_info()->segmentID( pos ) );
		}
		out << make_tag_with_dashes( resnum, chain, segid ) << std::endl;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_other_contacts( std::ostream & out, pose::Pose const & pose )
{
	using namespace core::scoring;
	using namespace core::conformation;
	using namespace utility::graph;
	TenANeighborGraph const & tenA( pose.energies().tenA_neighbor_graph() );

	static Distance const CONTACT_DIST_CUTOFF( 3.0 );

	std::map< Size, std::map< Size, bool > > base_base_interacting;
	for ( auto const & base_pair : rna_scoring_info_from_pose( pose ).rna_filtered_base_base_info().base_pair_list() ) {
		base_base_interacting[ base_pair.res1() ][ base_pair.res2() ] = true;
		base_base_interacting[ base_pair.res2() ][ base_pair.res1() ] = true;
	}
	for ( auto const & base_stack : rna_scoring_info_from_pose( pose ).rna_filtered_base_base_info().base_stack_list() ) {
		base_base_interacting[ base_stack.res1() ][ base_stack.res2() ] = true;
		base_base_interacting[ base_stack.res2() ][ base_stack.res1() ] = true;
	}


	// For each neighboring residue to the donor
	for ( EdgeListConstIterator ni = tenA.const_edge_list_begin(),
			ni_end = tenA.const_edge_list_end();
			ni != ni_end; ++ni ) {
		Size const & i = (*ni)->get_first_node_ind();
		Size const & j = (*ni)->get_second_node_ind();
		if ( abs( int(i) - int(j) ) <= 1 ) continue;
		if ( base_base_interacting[i][j] ) continue;

		// now look for a direct contact between these residues
		// might be faster to look up atom-level neighbor graph
		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );
		bool found_contact( false );
		Size m( 0 ), n( 0 );
		for ( m = 1; m <= rsd_i.nheavyatoms(); m++ ) {
			for ( n = 1; n <= rsd_j.nheavyatoms(); n++ ) {
				if ( ( rsd_i.xyz( m ) - rsd_j.xyz( n ) ).length() <  CONTACT_DIST_CUTOFF ) {
					out << pose.pdb_info()->tag( i ) << " " << pose.pdb_info()->tag( j ) <<  " " << rsd_i.atom_name(m) << " " << rsd_j.atom_name(n) << std::endl;
					found_contact = true;
					break;
				}
			}
			if ( found_contact ) break;
		}
		if ( !found_contact ) continue;

	}
}


//////////////////////////////////////////////////////////////////////////////////////
// @details remove Rosetta variants (preparing for output of FASTA)
void
remove_upper_lower_variants_from_RNA( pose::Pose & pose )
{
	using namespace core::chemical;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue_type( n ).is_RNA() ) {
			remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS_VARIANT, n );
			remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS_VARIANT, n );
		}
	}

}

} //ns rna
} //ns pose
} //ns core
