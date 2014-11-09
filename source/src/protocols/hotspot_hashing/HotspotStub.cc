// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pose/hotspot_hashing/HotspotStub.cc
/// @brief  HotspotStub class
/// @author John Karanicolas, Jacob Corn (jecorn@u.washington.edu), Sarel Fleishman

// Unit headers
#include <protocols/hotspot_hashing/HotspotStub.hh>

// Package headers
// AUTO-REMOVED #include <protocols/hotspot_hashing/HotspotStubSet.hh>

//Project Headers
// AUTO-REMOVED #include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <core/types.hh>
//

// C++ Headers
#include <map>


#include <core/pose/util.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace hotspot_hashing {

HotspotStub::HotspotStub(
	core::conformation::ResidueCOP const & residue,
	core::Real const bonus_value,
	core::pose::PoseOP pose,
	core::Size chain_to_design,
	protocols::filters::FilterCOP filter
) :
	utility::pointer::ReferenceCount(),
	residue_(core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue(*residue) ) )), // jk note: deep copy
	bonus_value_(bonus_value),
	pose_( pose ),
	filter_( filter ),
	chain_to_design_( chain_to_design )
//	my_set_( parent_set )
{ }

HotspotStub::HotspotStub( HotspotStub const & src ) :
	utility::pointer::ReferenceCount( src ),
	residue_( src.residue_ ),  // jk note: this is NOT a deep copy, since residue_ is effectively const
                            // (no member functions to move it). This must be a deep copy if residue_ can change.
	bonus_value_( src.bonus_value_ ),
	pose_( src.pose_ ),
	filter_( src.filter_ ),	// How should we deal with the reference to the containing set, since the new stub might be associ with a diff set?
	chain_to_design_( src.chain_to_design_ ),
//	my_set_( src.my_set_ ),
	scaffold_status_( src.scaffold_status_ )
{}

HotspotStub::~HotspotStub() {}

#ifndef BOINC // gives build error
#ifndef WIN32 // gives build error

HotspotStub & HotspotStub::operator=( HotspotStub const & src )
{
	residue_ = src.residue_;  // jk note: this is NOT a deep copy, since residue_ is effectively const
                            // (no member functions to move it). This must be a deep copy if residue_ can change.
	bonus_value_ = src.bonus_value_;

	// How should we deal with the reference to the containing set, since the new stub might be associ with a diff set?	my_set_ = src.my_set_ ;
//	my_set_ = src.my_set_;
	scaffold_status_ = src.scaffold_status_;
	pose_ = src.pose_;
	chain_to_design_ = src.chain_to_design_;
	filter_ = src.filter_;
	return *this;
}

bool HotspotStub::operator< ( HotspotStub const & right ) const {
	return ( bonus_value_ < right.bonus_value_ );
}

#endif
#endif

core::Real HotspotStub::bonus_value() const { return bonus_value_; }
core::conformation::ResidueCOP HotspotStub::residue() const { return residue_; }

/// @brief create scaffold_status_ appropriate for given scaffold and set it all to unchecked
void
HotspotStub::pair_with_scaffold( ) {
//	runtime_assert( my_set_ );
	core::Size const design_begin( pose_->conformation().chain_begin( chain_to_design_ ) );
	core::Size const design_end( pose_->conformation().chain_end( chain_to_design_ ) );
	scaffold_status_.clear();
	for (core::Size i = design_begin; i <= design_end; ++i ) {
		scaffold_status_.push_back( unchecked );
	}
}

void
HotspotStub::pair_with_scaffold( core::pose::PoseOP pose, protocols::filters::FilterCOP filter, core::Size chain_to_design )
{
	pose_ = pose;
	filter_ = filter;
	chain_to_design_ = chain_to_design;

	core::Size const design_begin( pose_->conformation().chain_begin( chain_to_design_ ) );
	core::Size const design_end( pose_->conformation().chain_end( chain_to_design_ ) );
	scaffold_status_.clear();
	for (core::Size i = design_begin; i <= design_end; ++i ) {
		scaffold_status_.push_back( unchecked );
	}
}

/// @brief Check whether stub matches a given position on the scaffold (must be within max_distance)
bool
HotspotStub::scaffold_match( core::Size const seqpos )
{
	Size const host_chain_begin( pose_->conformation().chain_begin( chain_to_design_) );
	Size const host_chain_end( pose_->conformation().chain_end( chain_to_design_) );
	runtime_assert( seqpos <= host_chain_end );
	runtime_assert( seqpos >= host_chain_begin );
	core::Size const scaffold_position = seqpos - host_chain_begin;
	if ( scaffold_status_[ scaffold_position ] == unchecked ) {
		using namespace core::conformation;
		ResidueCOP saved_res( ResidueOP( new Residue( pose_->residue( seqpos ) ) ) );
		ResidueCOP stub( residue_ );
		using namespace core::chemical;

		pose_->replace_residue( seqpos, *stub, true );
		if( seqpos > host_chain_begin )
			core::pose::remove_upper_terminus_type_from_pose_residue( *pose_, seqpos );
		if( seqpos < host_chain_end )
			core::pose::remove_lower_terminus_type_from_pose_residue( *pose_, seqpos );
		pose_->conformation().update_polymeric_connection( seqpos ); //o/w residue connections mess up
		if( scaffold_position > 0 )
			pose_->conformation().update_polymeric_connection( seqpos - 1 );
		pose_->update_residue_neighbors();
		if( filter_->apply( *pose_ ) )
			scaffold_status_[ scaffold_position ] = accept;
		else
			scaffold_status_[ scaffold_position ] = reject;
// SJF IMPORTANT! Next step is crucial b/c in current setup ALL stubs feed off the SAME POSE. Woe to those who forget
// to revert the pose to its original state!
		pose_->replace_residue( seqpos, *saved_res, true ); // SJF to revert the pose_ to its original state
		pose_->update_residue_neighbors();
	}
	return( scaffold_status_[ scaffold_position ] == accept );
}

/// @brief Get status, setting if unchecked
bool
HotspotStub::get_scaffold_status( core::Size const seqpos )
{
	runtime_assert( seqpos <= pose_->conformation().chain_end( chain_to_design_) );
	runtime_assert( seqpos >= pose_->conformation().chain_begin( chain_to_design_) );
	core::Size const scaffold_position = seqpos - pose_->conformation().chain_begin( chain_to_design_ );
	if( scaffold_status_[ scaffold_position ] != unchecked )
		return( scaffold_status_[ scaffold_position ] == accept );
	else
	{
		bool const status = scaffold_match( seqpos );
		return status;
	}
}

/// @brief manually set scaffold status
void
HotspotStub::set_scaffold_status( core::Size const seqpos, StubStatus const status)
{
	runtime_assert( seqpos <= pose_->conformation().chain_end( chain_to_design_) );
	runtime_assert( seqpos >= pose_->conformation().chain_begin( chain_to_design_) );
	core::Size const scaffold_position = seqpos - pose_->conformation().chain_begin( chain_to_design_ );
	scaffold_status_[ scaffold_position ] = status;
}

void HotspotStub::set_filter( protocols::filters::FilterCOP filter ) {
	filter_ = filter;
}

core::Size
HotspotStub::get_nearest_residue( core::pose::Pose const & pose ) const
{
	core::Size const chain_begin( pose.conformation().chain_begin( chain_to_design_ ) );
	core::Size const chain_end  ( pose.conformation().chain_end  ( chain_to_design_ ) );

	core::Real nearest_distance( 100000 );
	core::Size nearest_residue( 0 );
	for( core::Size resi( chain_begin ); resi<=chain_end; ++resi ){
    core::Real const distance( pose.residue( resi ).xyz( "CB" ).distance( residue()->xyz( "CB" ) ) );
		if( distance<=nearest_distance ){
			nearest_distance = distance;
			nearest_residue = resi;
		}
	}
	runtime_assert( nearest_residue );
	return( nearest_residue );
}

} // hotspot_hashing
} // protocols
