// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/AllowInsert.cc
/// @brief
/// @detailed
/// @author Rhiju Das

// Unit Headers
#include <protocols/toolbox/AllowInsert.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
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

static thread_local basic::Tracer TR( "protocols.toolbox.AllowInsert" );

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
// Note on domain code:
//
//     0      = moving
//    1,2,... = part of input pose 1, 2, ...
//   999      = FIXED
//
//    -- Rhiju Das, 2014
//
///////////////////////////////////////////////////////////////////////////

namespace protocols{
namespace toolbox{

	Size const FIXED_DOMAIN( 999 );

	AllowInsert::AllowInsert( core::pose::Pose const & pose ):
		nres_( pose.total_residue() ),
		force_ideal_chainbreak_( false )
	{
		initialize( pose );
	}

	//////////////////////////////////////////////////////////////////
	AllowInsert::~AllowInsert() {}

	//////////////////////////////////////////////////////////////////
	AllowInsertOP
	AllowInsert::clone() const
	{
		AllowInsertOP new_allow_insert( new AllowInsert( *this ) );
		return new_allow_insert;
	}

///////////////////////////////////////////////////////////////////////////////
	AllowInsert &
	AllowInsert::operator=( AllowInsert const & src )
	{
		if ( this == &src ) return *this;

		allow_insert_ = src.allow_insert_;
		named_atom_id_map_ = src.named_atom_id_map_;
		atom_ids_in_res_ = src.atom_ids_in_res_;

		calculated_atom_id_domain_map_ = src.calculated_atom_id_domain_map_;

		map_to_original_ = src.map_to_original_;

		nres_ = src.nres_;

		force_ideal_chainbreak_ = src.force_ideal_chainbreak_;

		return *this;
	}

///////////////////////////////////////////////////////////////////////////////
	AllowInsert::AllowInsert( AllowInsert const & src ):
		ReferenceCount()
	{
		*this = src;
	}



	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::initialize( core::pose::Pose const & pose )
	{

		using namespace core::pose::full_model_info;
		utility::vector1< Size > fixed_domain( pose.total_residue(), 0 );

		if ( full_model_info_defined( pose ) ) {
			// pose can store some information on separate domains... check inside.
			fixed_domain = get_fixed_domain_from_full_model_info_const( pose );
		}

		allow_insert_.clear();
		for (Size i = 1; i <= pose.total_residue(); i++ ) {

			utility::vector1< core::id::AtomID > atom_ids;

			for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {

				AtomID const atom_id( j, i );

				allow_insert_[ atom_id ] = fixed_domain[ i ];

				// Needed in case of changes in atom names/indices
				// The main allow_insert map is keyed on number but not names for speed.
				NamedAtomID const named_atom_id( pose.residue_type(i).atom_name( j ), i );
				named_atom_id_map_[ named_atom_id ] = atom_id;

				map_to_original_[ atom_id ] = atom_id;

				atom_ids.push_back( atom_id );

			}

			atom_ids_in_res_.push_back( atom_ids );

		}

		calculate_atom_id_domain_map( pose );
	}

	//////////////////////////////////////////////////////////////////
	bool
	AllowInsert::get( Size const & i ) const{
		return ( get_domain( i ) == 0 );
	}

	//////////////////////////////////////////////////////////////////
	bool
	AllowInsert::get( core::id::AtomID const & atom_id  ) const{
		return ( get_domain( atom_id ) == 0 );
	}

	//////////////////////////////////////////////////////////////////
	bool
	AllowInsert::get( core::id::TorsionID const & torsion_id, core::conformation::Conformation const & conformation ) const{

		// Check allow insert -- get atoms associated with these torsions.
		// are any allowed to move?
		id::AtomID id1,id2,id3,id4;

		bool const fail = conformation.get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

		if (fail) return false;

		//		TR << torsion_id << ": " << id1 << " " << get_domain( id1 ) << " -- " << id4 << " " << get_domain( id4 ) << std::endl;
		if  ( !get( id1 ) && !get( id4 ) && ( get_domain( id1 ) == get_domain( id4 ) ) ) return false;

		return true;

	}

	//////////////////////////////////////////////////////////////////
	bool
	AllowInsert::get_jump( Size const & jump_number, core::conformation::Conformation const & conformation ) const {

		id::AtomID id1, id2;
		bool const fail = conformation.get_jump_atom_ids( jump_number, id1, id2 );
		if ( fail ) return false;

		bool fixed_jump = ( !get( id1 ) && !get( id2 ) && ( get_domain( id1 ) == get_domain( id2 ) ) );
		//TR << TR.Cyan << jump_number << ": " << id1 << " " << get_domain( id1 ) << " -- " << id2 << " " << get_domain( id2 ) << "  MOVE? " << (!fixed_jump) << TR.Reset << std::endl;

		return ( !fixed_jump );
	}

	//////////////////////////////////////////////////////////////////
	Size
	AllowInsert::get_domain( Size const & i ) const{
		if ( i > nres_ ) {
			utility_exit_with_message( "Out of bounds for allow_insert" );
		}
		Size domain( 0 );
		for ( Size j = 1; j <= atom_ids_in_res_[ i ].size(); j++ ) {
			domain = get_domain( atom_ids_in_res_[ i ][ j ] );
			if ( domain == 0 ) return 0;
		}
		return domain;
	}

	//////////////////////////////////////////////////////////////////
	bool
	AllowInsert::has_domain( core::id::AtomID const & atom_id  ) const{
		std::map< AtomID, AtomID >::const_iterator it_original =	map_to_original_.find( atom_id );
		return !( it_original == map_to_original_.end() );
	}

	//////////////////////////////////////////////////////////////////
	Size
	AllowInsert::get_domain( core::id::AtomID const & atom_id  ) const{

		std::map< AtomID, AtomID >::const_iterator it_original =	map_to_original_.find( atom_id );
		if ( it_original == map_to_original_.end() ) return FIXED_DOMAIN;
		AtomID original_atom_id = it_original->second;

		std::map< core::id::AtomID, Size >::const_iterator it =	allow_insert_.find( original_atom_id );
		if ( it == allow_insert_.end() ){
			utility_exit_with_message( "Asked allow_insert for an atom_id it does not know about!" );
		}
		return it->second;
	}

	//////////////////////////////////////////////////////////////////
	Size
	AllowInsert::get_domain( core::id::NamedAtomID const & named_atom_id, pose::Pose const & pose  ) const{
		AtomID atom_id = named_atom_id_to_atom_id( named_atom_id, pose );
		return get_domain( atom_id );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set_domain( Size const & i, Size const & setting  ){
		for ( Size j = 1; j <= atom_ids_in_res_[ i ].size(); j++ ) {
			set_domain( atom_ids_in_res_[ i ][ j ], setting );
		}
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set_domain( core::id::AtomID const & atom_id, Size const & setting  ){
		std::map< AtomID, AtomID >::const_iterator it_original =	map_to_original_.find( atom_id );
		if ( it_original == map_to_original_.end() ) {
			std::cerr << "Problem ID: " << atom_id << std::endl;
			utility_exit_with_message( "Asked allow_insert to set atom_id that cannot be mapped to original pose!" );
		}

		AtomID original_atom_id = it_original->second;

		if ( allow_insert_.find( original_atom_id ) == allow_insert_.end() ){
			utility_exit_with_message( "Asked allow_insert to set atom_id it does not know about!" );
		}
		allow_insert_[ original_atom_id ] = setting;
	}


	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set_domain( core::id::NamedAtomID const & named_atom_id, pose::Pose const & pose, Size const & setting  ){
		AtomID atom_id = named_atom_id_to_atom_id( named_atom_id, pose );
		set_domain( atom_id, setting );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set_domain( Size const & setting  ){
		for( std::map< AtomID, Size >::iterator it = allow_insert_.begin();
				 it != allow_insert_.end(); it++ ){
			it->second = setting;
		}
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set_phosphate_domain( Size const & i,
																		 pose::Pose const & pose,
																		 Size const & setting ){

		if ( pose.residue(i).is_coarse() ){
			set_domain( AtomID( named_atom_id_to_atom_id( NamedAtomID( " P  ", i ), pose ) ), setting );
		} else {

			utility::vector1< std::string > const & atoms_involved = core::chemical::rna::atoms_involved_in_phosphate_torsion;
			for ( Size n = 1; n <= atoms_involved.size(); n++ ){
				set_domain( AtomID( named_atom_id_to_atom_id( NamedAtomID( atoms_involved[ n ], i ), pose ) ), setting );
			}

		}
	}


	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set( bool const & setting  ){
		set_domain( ( setting ) ? 0 : FIXED_DOMAIN  );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set( Size const & i, bool const & setting  ){
		set_domain( i, ( setting ) ? 0 : FIXED_DOMAIN  );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set( AtomID const & atom_id, bool const & setting  ){
		set_domain( atom_id, ( setting ) ? 0 : FIXED_DOMAIN  );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set( NamedAtomID const & named_atom_id, core::pose::Pose const & pose, bool const & setting  ){
		AtomID atom_id = named_atom_id_to_atom_id( named_atom_id, pose );
		set_domain( atom_id, ( setting ) ? 0 : FIXED_DOMAIN  );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::set_phosphate( Size const & i,
															pose::Pose const & pose,
															bool const & setting ){

		set_phosphate_domain( i, pose, ( ( setting ) ? 0 : FIXED_DOMAIN ) );
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::show() {
		for( Size i = 1; i <= nres_; i++ ) {
			std::cout << "RES" << i;
			for ( Size j = 1; j <= atom_ids_in_res_[ i ].size(); j++ ) {
				std::cout << ' ' << I( 3, allow_insert_[ atom_ids_in_res_[ i ][ j ] ] );
			}
			std::cout << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::and_allow_insert(AllowInsertOP allow_insert_in ){

		for( std::map< AtomID, Size >::iterator it = allow_insert_.begin();
				 it != allow_insert_.end(); it++ ){

			Size const & current_setting = it->second;
			Size const & other_setting = allow_insert_in->get_domain( it->first );

			if ( other_setting   > current_setting ) it->second = other_setting;

		}

	}

	//////////////////////////////////////////////////////////////////
	std::map< AtomID, Size > const &
	AllowInsert::calculated_atom_id_domain_map() {
		return calculated_atom_id_domain_map_;
	}

	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::calculate_atom_id_map( core::pose::Pose const & pose,
																			std::map< core::Size, core::Size > const & res_map,
																			core::kinematics::FoldTree const & scratch_fold_tree,
																			std::map< AtomID, AtomID > & atom_id_map  ){

		atom_id_map.clear();

		std::map< core::Size, core::Size > in_source_res; //basically reverse of res_map.
		for( std::map< Size, Size >::const_iterator  it  = res_map.begin(); it != res_map.end(); it++ ){
			Size const & insert_pos = it->first;
			Size const & source_pos = it->second;
			in_source_res[ source_pos ] = insert_pos;
		}

		for( std::map< Size, Size >::const_iterator  it  = res_map.begin(); it != res_map.end(); it++ ){

			Size const & insert_pos = it->first;
			Size const & source_pos = it->second;

			for ( Size j = 1; j <= pose.residue_type( insert_pos ).natoms(); j++ ){

				AtomID const atom_id( j, insert_pos );

				std::map< AtomID, AtomID >::const_iterator it_original =	map_to_original_.find( atom_id );
				if ( it_original == map_to_original_.end() )  continue;
				//				{
				//					utility_exit_with_message( "Asked allow_insert to set atom_id that cannot be mapped to original pose!" );
				//				}

				AtomID original_atom_id = it_original->second;

				int const rsd_offset = int( original_atom_id.rsd() ) - int( atom_id.rsd() ); //this is nonzero for OVL1, OVU1, etc. (chainbreak atoms)
				Size const source_atomno = original_atom_id.atomno();

				Size const source_pos_offset = source_pos + rsd_offset;

				if ( in_source_res.find( source_pos_offset ) == in_source_res.end() ) continue;

				if ( rsd_offset == +1 && scratch_fold_tree.is_cutpoint( source_pos   ) ) continue;
				if ( rsd_offset == -1 && scratch_fold_tree.is_cutpoint( source_pos-1 ) ) continue;

				// PUTTING IN MATT'S CRAZY OPTION FOR HOMOLOGY MODELING -- make this an option though!
				//if ( true ){
				//					if ( rsd_offset != 0) continue;
				//				}

				//				if ( source_pos + rsd_offset == 0 ) {
				//					std::cout << pose.annotated_sequence( true ) << std::endl;
				//					std::cout << "FAIL on Atom " << atom_id << "; ''original'' atom: " << original_atom_id << std::endl;
				//					utility_exit_with_message( "mapping atom_id to rsd 0?" );
				//				}

				AtomID source_atom_id( source_atomno, source_pos + rsd_offset );

				atom_id_map[ atom_id ] = source_atom_id;

				// change this to an assert later?
				//				if ( AtomID( named_atom_id_map_[ atom_id ], pose ) != atom_id  ) utility_exit_with_message( "allow insert atom_id mismatch" );

			}

		}


	}


	//////////////////////////////////////////////////////////////////
	void
	AllowInsert::renumber_after_variant_changes( core::pose::Pose const & pose ){


		if ( pose.total_residue() != nres_ ) {
			utility_exit_with_message( "AllowInsert cannot currenty handle changes in no. residues!" );
		}

		map_to_original_.clear();

		for( Size i = 1; i <= nres_; i++ ) {

			for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {

				AtomID const new_atom_id( j, i );

				std::string const & atom_name = pose.residue_type(i).atom_name( j );
				NamedAtomID const new_named_atom_id( atom_name, i );

				std::map< NamedAtomID, AtomID >::const_iterator it = named_atom_id_map_.find( new_named_atom_id );

				if ( it != named_atom_id_map_.end() ) { //awesome, this atom is recognizable by its name.

					AtomID const & original_atom_id = it->second;
					map_to_original_[ new_atom_id ] = original_atom_id;

				} else { // there are some special cases....

					NamedAtomID alternative_named_atom_id;

					// if we prevent mapping of chainbreak virtual atoms OVL1, OVL2, OVU1 to their
					// corresponding real atoms, they won't move during fragment insertions.
					if ( force_ideal_chainbreak_ ) continue;

					//note that this is hardcoded, but it is also hardcoded in Conformation.cc so I don't feel so bad.
					// later generalize based on mainchain[ ... ], so can handle any polymer with cutpoint variants.
					// later can add H1, H for protein terminus variants.
					core::chemical::ResidueType const & rsd_type = pose.residue_type( i );
					if ( rsd_type.is_RNA() ) {
						if ( rsd_type.is_coarse() ){
							if ( atom_name == "OVL1" ){
								alternative_named_atom_id = NamedAtomID( " P  ", i+1 );
							} else if ( atom_name == "OVL2" ) {
								alternative_named_atom_id = NamedAtomID( " S  ", i+1 );
							} else if ( atom_name == "OVU1" ) {
								alternative_named_atom_id = NamedAtomID( " S  ", i-1 );
							}	else {
								continue;
							}
						} else {
							if ( atom_name == "OVL1" ){
								alternative_named_atom_id = NamedAtomID( " P  ", i+1 );
							} else if ( atom_name == "OVL2" ) {
								alternative_named_atom_id = NamedAtomID( " O5'", i+1 );
							} else if ( atom_name == "OVU1" ) {
								alternative_named_atom_id = NamedAtomID( " O3'", i-1 );
							}	else {
								continue;
							}
						}
					}

					std::map< NamedAtomID, AtomID >::const_iterator it2 = named_atom_id_map_.find( alternative_named_atom_id );

					if ( it2 != named_atom_id_map_.end() ) { //awesome, this atom is recognizable by its name.
						AtomID const & original_atom_id = it2->second;
						map_to_original_[ new_atom_id ] = original_atom_id;
					}

				}
			}
		}

		// std::cout << "NEW ALLOW INSERT FOR CHANGED POSE " << std::endl;
		// for( Size i = 1; i <= nres_; i++ ) {
		// 	std::cout << "RES " << i << "  ";
		// 	for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {
		// 		if ( get_domain( AtomID(j,i) ) == 0 ) {
		// 			std::cout << ' ' <<  pose.residue_type( i ).atom_name( j ); //I(3,get_domain( AtomID( j,i ) ) );
		// 		}
		// 	}
		// 	std::cout << std::endl;
		// }

		calculate_atom_id_domain_map( pose );
	}

	//////////////////////////////////////////////////////////////////////////

	std::map< AtomID, Size > const &
	AllowInsert::calculate_atom_id_domain_map( core::pose::Pose const & pose ){

		calculated_atom_id_domain_map_.clear();

		for( Size i = 1; i <= pose.total_residue(); i++ ) {

			for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {

				AtomID const atom_id( j, i );

				calculated_atom_id_domain_map_[ atom_id ] = get_domain( atom_id );

			}

		}

		return calculated_atom_id_domain_map_;
	}

	//////////////////////////////////////////////////////////////////
	// Used in appending virtual residue to pose.
	void
	AllowInsert::append_residue( core::pose::Pose const & pose,
															 Size const & i,
															 bool const & setting ){

		using namespace core::id;
		using namespace core::chemical;

		Size const domain = setting ? 0 : FIXED_DOMAIN;

		utility::vector1< AtomID > atom_ids;
		for ( Size j = 1; j <= pose.residue( i ).natoms(); j++ ) {
			AtomID const atom_id( j, i );
			allow_insert_[ atom_id ] = domain;
			named_atom_id_map_[ atom_id_to_named_atom_id( atom_id, pose ) ] = atom_id;
			atom_ids.push_back( atom_id );
			map_to_original_[ atom_id ] = atom_id;
		}
		atom_ids_in_res_.push_back( atom_ids );
		nres_++;

	}

	////////////////////////////////////////////////////////////////////////////
	// whoa, this is easy now!
	////////////////////////////////////////////////////////////////////////////
	void
	AllowInsert::setup_movemap( core::kinematics::MoveMap & mm,
															core::pose::Pose const & pose ){

		using namespace core::id;

		runtime_assert( pose.total_residue() == nres() );

		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );

		for ( Size i = 1; i <= nres(); i++ ){

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

		for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			if ( get_jump( n, pose.conformation() ) ) {
				mm.set_jump( n, true );
			}
		}
	}

}
}



