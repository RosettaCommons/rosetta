// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/AtomLevelDomainMap.cc
/// @brief
/// @details
/// @author Rhiju Das

// Unit Headers
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/toolbox/AtomID_Mapper.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/chemical/rna/util.hh> // for information on phosphate atoms.

#include <utility/exit.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using ObjexxFCL::format::I;

// C++ headers
#include <map>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.AtomLevelDomainMap" );

using namespace core;
using core::id::AtomID;
using core::id::NamedAtomID;

///////////////////////////////////////////////////////////////////////////
//
// This object stores at the atom level information on domains in the pose,
//  with a code for parts that are definitely moving and parts that are
//  definitely fixed.
//
// Very useful for defining move-maps [including variable bond geometry.]
//
// Currently tied to a specific pose and atom numbering -- there's a mapping function for
//  when pose variants change, based on NamedAtomID, but could bear some improvement.
//
// NOTE: If pose changes any residue variants, must call renumber_after_variant_changes( pose )!
//
// Note on domain code:
//
//     0      = moving
//    1,2,... = part of input pose 1, 2, ...
//   999      = FIXED
//
//    -- Rhiju Das, 2014
//
///////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace toolbox {

//Size const FIXED_DOMAIN( 999 ); // now defined in CopyDofs.cc

AtomLevelDomainMap::AtomLevelDomainMap( core::pose::Pose const & pose,
	bool const map_to_vanilla_pose /* = false */,
	utility::vector1< Size > const & allow_insert_res /* = blank */ )
{
	initialize( pose, map_to_vanilla_pose, allow_insert_res );
}

//////////////////////////////////////////////////////////////////
AtomLevelDomainMap::~AtomLevelDomainMap() {}

//////////////////////////////////////////////////////////////////
AtomLevelDomainMapOP
AtomLevelDomainMap::clone() const
{
	AtomLevelDomainMapOP new_atom_level_domain_map( new AtomLevelDomainMap( *this ) );
	return new_atom_level_domain_map;
}

///////////////////////////////////////////////////////////////////////////////
AtomLevelDomainMap::AtomLevelDomainMap( AtomLevelDomainMap const & src ):
	ReferenceCount()
{
	*this = src;
}

//////////////////////////////////////////////////////////////////
bool
AtomLevelDomainMap::get( Size const & i ) const{
	return ( get_domain( i ) == 0 );
}

//////////////////////////////////////////////////////////////////
bool
AtomLevelDomainMap::get( core::id::AtomID const & atom_id  ) const{
	return ( get_domain( atom_id ) == 0 );
}

//////////////////////////////////////////////////////////////////
bool
AtomLevelDomainMap::get( core::id::TorsionID const & torsion_id, core::conformation::Conformation const & conformation ) const{

	// Check allow insert -- get atoms associated with these torsions.
	// are any allowed to move?
	id::AtomID id1,id2,id3,id4;

	bool const fail = conformation.get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

	if ( fail ) return false;

	if  ( !get( id1 ) && !get( id4 ) && ( get_domain( id1 ) == get_domain( id4 ) ) ) return false;

	return true;

}

//////////////////////////////////////////////////////////////////
bool
AtomLevelDomainMap::get_jump( Size const & jump_number, core::conformation::Conformation const & conformation ) const {

	id::AtomID id1, id2;
	bool const fail = conformation.get_jump_atom_ids( jump_number, id1, id2 );
	if ( fail ) return false;

	bool fixed_jump = ( !get( id1 ) && !get( id2 ) && ( get_domain( id1 ) == get_domain( id2 ) ) );

	return ( !fixed_jump );
}

//////////////////////////////////////////////////////////////////
Size
AtomLevelDomainMap::get_domain( Size const & i ) const{
	if ( i > atom_id_mapper_->nres() ) {
		utility_exit_with_message( "Out of bounds for atom_level_domain_map" );
	}
	Size domain( 0 );
	utility::vector1< AtomID > atom_ids_in_res( atom_id_mapper_->atom_ids_in_res( i ) );
	for ( Size j = 1; j <= atom_ids_in_res.size(); j++ ) {
		domain = get_domain( atom_ids_in_res[ j ] );
		if ( domain == 0 ) return 0;
	}
	return domain;
}

//////////////////////////////////////////////////////////////////
bool
AtomLevelDomainMap::has_domain( core::id::AtomID const & atom_id ) const{
	return atom_id_mapper_->has_atom_id( atom_id );
}

//////////////////////////////////////////////////////////////////
Size
AtomLevelDomainMap::get_domain( core::id::AtomID const & atom_id  ) const {
	if ( !atom_id_mapper_->has_atom_id( atom_id ) ) return FIXED_DOMAIN;
	std::map< core::id::AtomID, Size >::const_iterator it = domain_map_.find( atom_id_mapper_->map_to_reference( atom_id ) );
	if ( it == domain_map_.end() )  utility_exit_with_message( "Asked atom_level_domain_map for an atom_id it does not know about!" );
	return it->second;
}

//////////////////////////////////////////////////////////////////
Size
AtomLevelDomainMap::get_domain( core::id::NamedAtomID const & named_atom_id, pose::Pose const & pose  ) const{
	AtomID atom_id = named_atom_id_to_atom_id( named_atom_id, pose );
	return get_domain( atom_id );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_domain( Size const & i, Size const & setting  ){
	utility::vector1< AtomID > atom_ids_in_res( atom_id_mapper_->atom_ids_in_res( i ) );
	for ( Size j = 1; j <= atom_ids_in_res.size(); j++ ) {
		set_domain( atom_ids_in_res[ j ], setting );
	}
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_domain( core::id::AtomID const & atom_id, Size const & setting  ){
	if ( !atom_id_mapper_->has_atom_id( atom_id ) ) {
		std::cerr << "Problem ID: " << atom_id << std::endl;
		utility_exit_with_message( "Asked atom_level_domain_map to set atom_id that cannot be mapped to original pose!" );
	}
	AtomID original_atom_id = atom_id_mapper_->map_to_reference( atom_id );
	if ( domain_map_.find( original_atom_id ) == domain_map_.end() )  utility_exit_with_message( "Asked atom_level_domain_map to set atom_id it does not know about!" );
	domain_map_[ original_atom_id ] = setting;
}


//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_domain( core::id::NamedAtomID const & named_atom_id, pose::Pose const & pose, Size const & setting  )
{
	AtomID atom_id = named_atom_id_to_atom_id( named_atom_id, pose );
	set_domain( atom_id, setting );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_domain( Size const & setting  )
{
	for ( std::map< AtomID, Size >::iterator it = domain_map_.begin(),
			end = domain_map_.end(); it != end; ++it ) {
		it->second = setting;
	}
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_phosphate_domain( Size const & i,
	pose::Pose const & pose,
	Size const & setting )
{
	if ( pose.residue(i).is_coarse() ) {
		set_domain( AtomID( named_atom_id_to_atom_id( NamedAtomID( " P  ", i ), pose ) ), setting );
	} else {

		utility::vector1< std::string > const & atoms_involved = core::chemical::rna::atoms_involved_in_phosphate_torsion;
		for ( Size n = 1; n <= atoms_involved.size(); n++ ) {
			set_domain( AtomID( named_atom_id_to_atom_id( NamedAtomID( atoms_involved[ n ], i ), pose ) ), setting );
		}

	}
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_sugar_domain( Size const & i,
	pose::Pose const & pose,
	Size const & setting ){

	runtime_assert( !pose.residue(i).is_coarse() );
	utility::vector1< std::string > const & atoms_involved = core::chemical::rna::sugar_atoms;
	for ( Size n = 1; n <= atoms_involved.size(); n++ ) {
		set_domain( AtomID( named_atom_id_to_atom_id( NamedAtomID( atoms_involved[ n ], i ), pose ) ), setting );
	}

}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set( bool const & setting  ){
	set_domain( ( setting ) ? 0 : FIXED_DOMAIN  );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set( Size const & i, bool const & setting  ){
	set_domain( i, ( setting ) ? 0 : FIXED_DOMAIN  );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_fixed_if_moving( AtomID const & atom_id ){
	if ( get( atom_id ) ) set_domain( atom_id, FIXED_DOMAIN );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_fixed_if_moving( Size const & i ){
	utility::vector1< AtomID > atom_ids_in_res( atom_id_mapper_->atom_ids_in_res( i ) );
	for ( Size j = 1; j <= atom_ids_in_res.size(); j++ ) {
		set_fixed_if_moving( atom_ids_in_res[ j ]);
	}
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set( AtomID const & atom_id, bool const & setting  ){
	set_domain( atom_id, ( setting ) ? 0 : FIXED_DOMAIN  );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set( NamedAtomID const & named_atom_id, core::pose::Pose const & pose, bool const & setting  ){
	AtomID atom_id = named_atom_id_to_atom_id( named_atom_id, pose );
	if ( setting ) {
		set_domain( atom_id, 0 );
	} else {
		set_fixed_if_moving( atom_id );
	}
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_phosphate( Size const & i,
	pose::Pose const & pose,
	bool const & setting ){

	set_phosphate_domain( i, pose, ( ( setting ) ? 0 : FIXED_DOMAIN ) );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::set_sugar( Size const & i,
	pose::Pose const & pose,
	bool const & setting ){

	set_sugar_domain( i, pose, ( ( setting ) ? 0 : FIXED_DOMAIN ) );
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::show() const {
	for ( Size i = 1; i <= atom_id_mapper_->nres(); i++ ) {
		std::cout << "RES" << i;
		utility::vector1< core::id::AtomID > const & atom_ids_in_res( atom_id_mapper_->atom_ids_in_res( i ) );
		for ( Size j = 1; j <= atom_ids_in_res.size(); j++ ) {
			std::cout << ' ' << I( 3, get_domain( atom_ids_in_res[ j ] ) );
		}
		std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////
std::map< AtomID, Size >
AtomLevelDomainMap::calculate_atom_id_domain_map( core::pose::Pose const & pose ) const
{
	std::map< core::id::AtomID, Size > calculated_atom_id_domain_map;

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {

			AtomID const atom_id( j, i );

			calculated_atom_id_domain_map[ atom_id ] = get_domain( atom_id );

		}
	}

	return calculated_atom_id_domain_map;
}

////////////////////////////////////////////////////////////////////////////
// whoa, this is easy now!
////////////////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::setup_movemap( core::kinematics::MoveMap & mm,
	core::pose::Pose const & pose ){

	using namespace core::id;

	runtime_assert( pose.total_residue() == atom_id_mapper_->nres() );

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		utility::vector1< TorsionID > torsion_ids;
		for ( Size torsion_number = 1; torsion_number <= pose.residue( i ).mainchain_torsions().size(); torsion_number++ ) {
			torsion_ids.push_back( TorsionID( i, id::BB, torsion_number ) );
		}
		for ( Size torsion_number = 1; torsion_number <= pose.residue_type( i ).nchi(); torsion_number++ ) {
			torsion_ids.push_back( TorsionID( i, id::CHI, torsion_number ) );
		}

		for ( Size n = 1; n <= torsion_ids.size(); n++ ) {
			TorsionID const & torsion_id  = torsion_ids[ n ];
			if ( get( torsion_id, pose.conformation() ) ) mm.set( torsion_id, true );
		}
	}

	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
		if ( get_jump( n, pose.conformation() ) ) {
			mm.set_jump( n, true );
		}
	}
}

void
AtomLevelDomainMap::renumber_after_variant_changes( core::pose::Pose const & pose )
{
	AtomID_MapperOP atom_id_mapper_new( atom_id_mapper_->clone() );
	atom_id_mapper_new->renumber_after_variant_changes( pose );
	atom_id_mapper_ = atom_id_mapper_new; // save as COP
}

//////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::initialize( core::pose::Pose const & pose,
	bool const map_to_vanilla_pose,
	utility::vector1< Size > const & allow_insert_res )
{
	using namespace core::pose;
	using namespace core::pose::full_model_info;

	// need a map to a reference pose without variants.
	atom_id_mapper_ = AtomID_MapperCOP( new AtomID_Mapper( pose, map_to_vanilla_pose ) );

	utility::vector1< Size > fixed_domain( pose.total_residue(), 0 );
	if ( full_model_info_defined( pose ) ) {
		// Pose can store some information on separate domains... check inside.
		// Any domains that are not claimed as fixed_domains will be assigned domain 0 (i.e., free)
		fixed_domain = get_fixed_domain_from_full_model_info_const( pose );
		runtime_assert( allow_insert_res.size() == 0 ); // allow_insert_res is a legacy of rna_denovo with params files.
	}

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {
			AtomID const atom_id( j, i );
			if ( atom_id_mapper_->has_atom_id( atom_id ) ) {
				domain_map_[ atom_id_mapper_->map_to_reference( atom_id ) ] = fixed_domain[ i ];
			}
		}
	}

	if ( allow_insert_res.size() > 0 ) apply_allow_insert_res( allow_insert_res );

	// following is totally hacky, but I've got to put it somewhere -- rhiju.
	update_to_not_move_virtual_phosphates( pose );
	update_to_not_move_last_virtual_residue( pose );
	update_to_move_internal_phosphates( pose );

}

////////////////////////////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::update_to_not_move_last_virtual_residue( pose::Pose const & pose )
{
	Size n = pose.total_residue();
	if ( pose.residue_type( n ).aa() == core::chemical::aa_vrt && get( n ) /* moving */ ) set( n, false );
}


////////////////////////////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::update_to_not_move_virtual_phosphates( pose::Pose const & pose )
{
	using namespace chemical;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( pose.residue_type( n ).has_variant_type( VIRTUAL_PHOSPHATE ) ) {
			set_phosphate( n, pose, false );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::update_to_move_internal_phosphates( pose::Pose const & pose ) {
	// new -- make sure loops are moveable (& closeable!) at the 3'-endpoint
	for ( Size i = 1; i < pose.total_residue(); i++ ) {
		if ( pose.residue_type( i   ).is_RNA() &&
				pose.residue_type( i+1 ).is_RNA() &&
				get_domain( id::NamedAtomID( " C1'", i ), pose ) == 0 &&
				!pose.residue_type( i+1 ).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) ) {
			set_phosphate( i+1, pose, true );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
AtomLevelDomainMap::apply_allow_insert_res( utility::vector1< Size > const & allow_insert_res )
{
	set( false );
	for ( Size n = 1; n <= allow_insert_res.size(); n++ ) {
		Size const i = allow_insert_res[ n ];
		set( i, true );
	}
}


}
}


