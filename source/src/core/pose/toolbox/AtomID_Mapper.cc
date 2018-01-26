// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/toolbox/AtomID_Mapper.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/toolbox/AtomID_Mapper.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/methods/chainbreak_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.toolbox.AtomID_Mapper" );

using namespace core::id;
using namespace core::pose;

///////////////////////////////////////////////////////////////////////////
/// @details
///
/// Developed for RNA_DeNovo (farna) applications; also in use (as part
///  of AtomLevelDomainMap) in ConstrainToIdealMover.
///
/// A lot of times we need to map dofs from fragments or from 'chunks' of
///  library PDBs into our target pose. Its usually pretty easy
///  to make this map by matching atom_names, but that takes time due
///  to string lookups. This object helps cache those lookups.
///
/// In addition, the object formally holds some non-trivial mappings of
///  chainbreak atoms (OVL1, OVL2, OVU1) to their cognate regular atoms
///  in the pose.
///
/// Typical use case is the following:
///
///   AtomID_Mapper atom_id_mapper( reference_pose );
///   atom_id_mapper.renumber_after_variant_changes( new_pose );
///    ...
///   desired_atom_id_in_new_pose = atom_id_mapper.map_to_reference( some_atom_id_in_new_pose );
///
///            -- Rhiju, 2015
///
///////////////////////////////////////////////////////////////////////////

namespace core {
namespace pose {
namespace toolbox {

//Constructor
AtomID_Mapper::AtomID_Mapper( core::pose::Pose const & pose,
	bool const map_to_vanilla_pose /* = false */  )
{
	initialize( pose, map_to_vanilla_pose );
}

//constructor
AtomID_Mapper::AtomID_Mapper( AtomID_Mapper const & src ) {
	*this = src;
}

//constructor
AtomID_MapperOP
AtomID_Mapper::clone() const {
	return AtomID_MapperOP( new AtomID_Mapper( *this ) );
}


//Destructor
AtomID_Mapper::~AtomID_Mapper() = default;


//////////////////////////////////////////////////////////////////
void
AtomID_Mapper::initialize( core::pose::Pose const & pose, bool const map_to_vanilla_pose ) {
	if ( map_to_vanilla_pose ) {
		// following is pretty awful -- atom_level_domain_map makes use of atom indices for speed, and
		//  some downstream applications (RNA_ChunkLibrary) assume that it was set up in a pose
		//  without weird variants.
		//  Pose pose_without_variants;
		//
		//  // Want to keep any variants for protein stuff though
		//  std::string seq_without_RNA_variants;
		//  for ( core::Size i = 1; i <= pose.total_residue(); ++i ){
		//   char c = pose.residue( i ).name1();
		//   seq_without_RNA_variants += c;
		//   if (!pose.residue( i ).is_RNA()) {
		//    // get rid of the variants for RNA
		//    // keep the variants for anything other than RNA
		//    if (( !core::chemical::oneletter_code_specifies_aa(c) || core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code(c) ) != pose.residue(i).name() )) {
		//     seq_without_RNA_variants += '[';
		//     seq_without_RNA_variants += pose.residue(i).name();
		//     seq_without_RNA_variants += ']';
		//    }
		//   }
		//  }
		//  make_pose_from_sequence( pose_without_variants, seq_without_RNA_variants,
		//   pose.residue_type( 1 ).residue_type_set(), false /*auto_termini*/ );
		// try Andy's fix instead

		// AHA! Here is the issue. We NEED to have some info from annotated_sequence()
		// because we need the base name of everything.
		Pose pose_without_variants = pose;
		//make_pose_from_sequence( pose_without_variants, pose.sequence() /* note this is not annotated_sequence(), which would include variants*/,
		// pose.residue_type( 1 ).residue_type_set(), false /*auto_termini*/ );
		using namespace core::conformation;
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {

			if ( pose.residue( ii ).is_polymer() && !pose.residue( ii ).is_protein() ) {
				// If this residue type does not havea base name, just use the name.
				// With any luck, these will mostly be CCD-generated residues... which
				// may have few variants anyway!
				if ( pose.residue_type( ii ).base_name() != "" ) {
					ResidueOP new_rsd = ResidueFactory::create_residue( pose.residue_type_set_for_pose( pose.residue_type(ii).mode() )->name_map( pose.residue_type( ii ).base_name() ) );
					pose_without_variants.replace_residue( ii, *new_rsd, true );
				}
			}
		}
		initialize_from_pose( pose_without_variants );
		renumber_after_variant_changes( pose );
	} else {
		initialize_from_pose( pose );
	}
}

//////////////////////////////////////////////////////////////////
void
AtomID_Mapper::initialize_from_pose( core::pose::Pose const & pose ) {

	for ( Size i = 1; i <= pose.size(); i++ ) {

		utility::vector1< AtomID > atom_ids;

		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {

			AtomID atom_id( j, i );

			// Needed in case of changes in atom names/indices
			// The main atom_level_domain_map map is keyed on number but not names for speed.
			NamedAtomID const named_atom_id( pose.residue_type(i).atom_name( j ), i );
			named_atom_id_map_[ named_atom_id ] = atom_id;

			map_to_reference_[ atom_id ]   = atom_id; /* primary map in use */
			map_from_reference_[ atom_id ] = atom_id; /* reverse map */

			atom_ids.push_back( atom_id );
		}

		atom_ids_in_res_.push_back( atom_ids );
	}

	sequence_ = pose.sequence();

}


//////////////////////////////////////////////////////////////////
bool
AtomID_Mapper::has_atom_id( AtomID const & atom_id ) const
{
	auto it_reference = map_to_reference_.find( atom_id );
	return ( it_reference != map_to_reference_.end() );
}

//////////////////////////////////////////////////////////////////
AtomID const &
AtomID_Mapper::map_to_reference( AtomID const & atom_id ) const
{
	auto it_reference = map_to_reference_.find( atom_id );
	runtime_assert( it_reference != map_to_reference_.end() );
	return ( it_reference->second );
}

//////////////////////////////////////////////////////////////////
AtomID const &
AtomID_Mapper::map_from_reference( AtomID const & atom_id ) const
{
	auto it_reference = map_from_reference_.find( atom_id );
	runtime_assert( it_reference != map_from_reference_.end() );
	return ( it_reference->second );
}

//////////////////////////////////////////////////////////////////
std::map< AtomID, AtomID >
AtomID_Mapper::calculate_atom_id_map(
	core::pose::Pose const & target_pose,
	std::map< core::Size, core::Size > const & res_map /* from target to source */,
	core::kinematics::FoldTree const & source_fold_tree,
	AtomID_MapperCOP source_mapper_to_vanilla /* = 0 */ ) const
{
	std::map< AtomID, AtomID > atom_id_map;

	std::map< core::Size, core::Size > in_source_res; //basically reverse of res_map.
	for ( auto const & it : res_map ) {
		Size const & target_pos = it.first;
		Size const & source_pos = it.second;
		in_source_res[ source_pos ] = target_pos;
	}

	for ( auto const & it : res_map ) {

		Size const & target_pos = it.first;
		Size const & source_pos = it.second;

		for ( Size j = 1; j <= target_pose.residue_type( target_pos ).natoms(); j++ ) {

			AtomID const target_atom_id( j, target_pos );

			auto it_reference = map_to_reference_.find( target_atom_id );
			if ( it_reference == map_to_reference_.end() )  continue;

			AtomID reference_atom_id = it_reference->second;

			int const rsd_offset = int( reference_atom_id.rsd() ) - int( target_atom_id.rsd() ); //this is nonzero for OVL1, OVU1, etc. (chainbreak atoms)
			Size const source_atomno = reference_atom_id.atomno();

			Size const source_pos_offset = source_pos + rsd_offset; // will not be correct for cyclization

			if ( in_source_res.find( source_pos_offset ) == in_source_res.end() ) continue;

			// will not be correct for cyclization:
			if ( rsd_offset == +1 && source_fold_tree.is_cutpoint( source_pos   ) ) continue;
			if ( rsd_offset == -1 && source_fold_tree.is_cutpoint( source_pos-1 ) ) continue;

			AtomID source_atom_id( source_atomno, source_pos + rsd_offset );
			if ( source_mapper_to_vanilla != nullptr ) {
				// no mapping. happens with 2'-OH in DNA source poses.
				if ( !source_mapper_to_vanilla->in_map_from_reference( source_atom_id ) ) continue;
				source_atom_id = source_mapper_to_vanilla->map_from_reference( source_atom_id ); /* note 'reverse': from a vanilla pose to source pose. */
			}

			atom_id_map[ target_atom_id ] = source_atom_id;
		}

	}

	return atom_id_map;

}


//////////////////////////////////////////////////////////////////
// This is *really* specific to RNA. This basically assumes that the initial pose has no
//  chainbreak variants and has pretty standard numbering scheme.
void
AtomID_Mapper::renumber_after_variant_changes( core::pose::Pose const & pose ){

	using namespace core::scoring::methods;

	if ( pose.sequence() != sequence_ &&
			utility::replace_in( sequence_, "t", "X" ) !=
			utility::replace_in( pose.sequence(), "t", "X" ) /*horrible hack*/ ) {
		utility_exit_with_message( "AtomID_Mapper cannot currenty handle changes in sequence, just changes in variants! Your old sequence was " + sequence_ + " and the new sequence is " + pose.sequence() );
	}

	map_to_reference_.clear();
	map_from_reference_.clear();
	atom_ids_in_res_.clear();

	for ( Size i = 1; i <= pose.size(); i++ ) {

		utility::vector1< AtomID > atom_ids;

		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {

			AtomID const new_atom_id( j, i );

			std::string const & atom_name = pose.residue_type(i).atom_name( j );
			NamedAtomID const new_named_atom_id( atom_name, i );

			std::map< NamedAtomID, AtomID >::const_iterator it = named_atom_id_map_.find( new_named_atom_id );

			if ( it != named_atom_id_map_.end() ) { //awesome, this atom is recognizable by its name.

				AtomID const & reference_atom_id = it->second;
				map_to_reference_  [ new_atom_id      ] = reference_atom_id;
				map_from_reference_[ reference_atom_id ] = new_atom_id;

			} else { // there are some special cases....

				NamedAtomID alternative_named_atom_id;

				//note that this is hardcoded, but it is also hardcoded in Conformation.cc so I don't feel so bad.
				// later generalize based on mainchain[ ... ], so can handle any polymer with cutpoint variants.
				// later can add H1, H for protein terminus variants.
				core::chemical::ResidueType const & rsd_type = pose.residue_type( i );
				if ( rsd_type.is_RNA() ) {
					if ( rsd_type.is_coarse() ) {
						if ( atom_name == "OVL1" ) {
							alternative_named_atom_id = NamedAtomID( " P  ", get_upper_cutpoint_partner_for_lower( pose, i ) );
						} else if ( atom_name == "OVL2" ) {
							alternative_named_atom_id = NamedAtomID( " S  ", get_upper_cutpoint_partner_for_lower( pose, i ) );
						} else if ( atom_name == "OVU1" ) {
							alternative_named_atom_id = NamedAtomID( " S  ", get_lower_cutpoint_partner_for_upper( pose, i ) );
						} else {
							continue;
						}
					} else {
						if ( atom_name == "OVL1" ) {
							alternative_named_atom_id = NamedAtomID( " P  ", get_upper_cutpoint_partner_for_lower( pose, i ) );
						} else if ( atom_name == "OVL2" ) {
							alternative_named_atom_id = NamedAtomID( " O5'", get_upper_cutpoint_partner_for_lower( pose, i ) );
						} else if ( atom_name == "OVU1" ) {
							alternative_named_atom_id = NamedAtomID( " O3'", get_lower_cutpoint_partner_for_upper( pose, i ) );
						} else {
							continue;
						}
					}
				}

				std::map< NamedAtomID, AtomID >::const_iterator it2 = named_atom_id_map_.find( alternative_named_atom_id );

				if ( it2 != named_atom_id_map_.end() ) { //awesome, this atom is recognizable by its name.
					AtomID const & reference_atom_id = it2->second;
					map_to_reference_[ new_atom_id ] = reference_atom_id;
					// we will not update map_from_reference -- assume those reference atoms will go to their namesakes, not these alternatives.
				}

			}

			if ( map_to_reference_.find( new_atom_id ) != map_to_reference_.end() ) atom_ids.push_back( new_atom_id );
		}

		atom_ids_in_res_.push_back( atom_ids );
	}

}

} //toolbox
} //toolbox
} //protocols
