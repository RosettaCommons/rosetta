// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/scoring/InterfaceInfo.cc
/// @brief  Statistically derived rotamer pair potentials
/// @details For docking (or between chains) only those residues at the interface
///      and between the two interfaces need to be evaluated
/// @author Monica Berrondo


// Unit headers
#include <protocols/scoring/InterfaceInfo.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Basic headers

// Utility headers
#include <utility/vector1.hh>

#include <protocols/scoring/Interface.hh> // AUTO IWYU For Interface


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace scoring {

/// @brief Default constructor with rb jump = 1
InterfaceInfo::InterfaceInfo()
{
	rb_jump_.push_back( 1 );
	initialize();
}

/// @brief Constructor with arguments for non-default rb jump
InterfaceInfo::InterfaceInfo( core::Size rb_jump_in )
{
	rb_jump_.push_back( rb_jump_in );
	initialize();
}

/// @brief Constructor with arguments for multiple jumps
InterfaceInfo::InterfaceInfo( utility::vector1_size const & rb_jump_in )
{
	rb_jump_ = rb_jump_in;
	initialize();
}

/// @brief Copy constructor
InterfaceInfo::InterfaceInfo( InterfaceInfo const & src ) :
	CacheableData(src)
{
	rb_jump_ = src.rb_jump_;
	num_jump_ = src.num_jump_;
	distance_ = src.distance_;

	interface_list_ = src.interface_list_;
	initialize();
}

void
InterfaceInfo::setup_interfaces()
{
	interface_list_.resize(num_jump_);

	//initialize interface objects for each interface in pose
	for ( core::Size i = 1; i <= num_jump_; i++ ) {
		interface_list_[i] = utility::pointer::make_shared< protocols::scoring::Interface >( rb_jump_[i] );
		interface_list_[i]->distance(distance_);
	}
}

void
InterfaceInfo::initialize()
{
	//get the number of jumps in the fold tree
	num_jump_ = rb_jump_.size();

	// Initialize the interface objects
	setup_interfaces();
}

bool
InterfaceInfo::is_interface(
	core::conformation::Residue rsd
) const
{

	bool is_interface(false);
	for ( core::Size i = 1; i <= num_jump_; i++ ) {
		// If the seqpos is within the size of the pose (as understood by the Interface), check if the rsd is_interface
		if ( rsd.seqpos() <= interface_list_[i]->contact_list().size() ) {
			if ( interface_list_[i]->is_interface( rsd ) ) is_interface = true;
		}
	}

	return is_interface;
}

bool
InterfaceInfo::is_pair(
	core::conformation::Residue rsd1,
	core::conformation::Residue rsd2
) const
{

	bool is_pair(false);

	for ( core::Size i = 1; i<= num_jump_; i++ ) {
		if ( interface_list_[i]->is_pair(rsd1, rsd2) ) is_pair = true;
	}

	return is_pair;
}

core::Size
InterfaceInfo::interface_nres( core::Size jump_num ) const
{
	return interface_list_[jump_num]->interface_nres();
}

InterfaceCOP
InterfaceInfo::interface( core::Size interface_num ) const {
	return interface_list_[ interface_num ];
}

/// @brief Removes all jumps from the interface calculation
void
InterfaceInfo::clear_jumps()
{
	rb_jump_.clear();
	initialize();
}

/// @brief Adds another jump to the interface calculation, for example
///for multi-body docking
void
InterfaceInfo::add_jump( core::Size jump_in )
{
	rb_jump_.push_back( jump_in );
	initialize();
}

/// @brief Sets the distance cutoff for interface calculations, and updates interfaces accordingly
void
InterfaceInfo::distance( core::Real distance_in )
{
	if ( distance_in != distance_ ) {
		distance_ = distance_in;

		// Update the interfaces using the new distance
		setup_interfaces();
	}
}

void
InterfaceInfo::calculate( core::pose::Pose const & pose )
{
	for ( core::Size i = 1; i <= num_jump_; i++ ) {
		// Only calculate the interface if we have enough jumps in the pose
		if ( pose.num_jump() >= i ) {
			interface_list_[i]->calculate( pose );
		}
	}
}

core::Size
InterfaceInfo::closest_interface_residue( core::pose::Pose const & pose, core::Size src_rsd, core::Real & distance ) const
{
	core::Size ret_rsd (0);
	core::Real min_distance (100000.0), temp_dist (0.0);
	for ( core::Size i = 1; i <= num_jump_; ++i ) {
		core::Size temp_rsd = interface_list_[i]->closest_interface_residue( pose, src_rsd, temp_dist );
		if ( temp_dist < min_distance ) {
			ret_rsd = temp_rsd;
			min_distance = temp_dist;
		}
	}
	distance = temp_dist;
	return ret_rsd;
}

}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::scoring::InterfaceInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( calculated_ ) ); // _Bool
	arc( CEREAL_NVP( num_jump_ ) ); // core::Size
	arc( CEREAL_NVP( distance_ ) ); // core::Real
	arc( CEREAL_NVP( rb_jump_ ) ); // utility::vector1_size
	arc( CEREAL_NVP( interface_list_ ) ); // utility::vector1<InterfaceOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::scoring::InterfaceInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( calculated_ ); // _Bool
	arc( num_jump_ ); // core::Size
	arc( distance_ ); // core::Real
	arc( rb_jump_ ); // utility::vector1_size
	arc( interface_list_ ); // utility::vector1<InterfaceOP>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::scoring::InterfaceInfo );
CEREAL_REGISTER_TYPE( protocols::scoring::InterfaceInfo )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_scoring_InterfaceInfo )
#endif // SERIALIZATION
