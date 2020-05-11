// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/DEERDataCache.cc
/// @brief  Contains data for DEER energy method, stored as DEERData objects
/// @details To prevent the energy method from reading from the command line every scoring
///      round and parsing the input file, this method stores a list of DEER decay data.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/DEERIO.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <set>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Destructor
DEERDataCache::~DEERDataCache() {}

/// @brief Copy function, overrides parent class
basic::datacache::CacheableDataOP
DEERDataCache::clone() const {
	return basic::datacache::CacheableDataOP( new DEERDataCache( *this ) );
}

/// @brief Returns non-const data for a given ID
DEERDataOP &
DEERDataCache::operator[](
	Size const & i
) {
	return data_[ i ];
}

/// @brief Returns const data for a given ID
DEERDataOP const &
DEERDataCache::at(
	Size const & i
) const {
	return data_.at( i );
}

void
DEERDataCache::append(
	DEERDataOP const & data
) {
	data_[ size() + 1 ] = data;
}

/// @brief Returns number of DEERData objects stored here
Size
DEERDataCache::size() const {
	if ( data_.size() == 0 ) {
		return 0;
	} else {
		return data_.rbegin()->first;
	}
}

/// @brief Fetches are parses data from command line using DEERIO object
void
DEERDataCache::fetch_and_organize_data(
	pose::Pose & pose
) {
	// First get the data from the command line via DEERIO object
	data_ = DEERIO().generate_data( pose );

	// Go through the residues and store that info
	for ( auto const & edge : data_ ) {
		auto const & deer_res = edge.second->residues();
		for ( auto const & res : deer_res ) {
			labeled_residues_.insert( res );
		}
		// Store all combinations of residue pairs
		for ( Size i = 1; i < deer_res.size(); ++i ) {
			Real weight_per_pair = ( deer_res.size() * ( deer_res.size() - 1.0 ) ) / 2.0;
			for ( Size j = i + 1; j <= deer_res.size(); ++j ) {
				auto const & min_res = std::min( deer_res.at( i ).first, deer_res.at( j ).first );
				auto const & max_res = std::max( deer_res.at( i ).first, deer_res.at( j ).first );
				auto res_pair = std::make_pair( min_res, max_res );
				if ( pair_to_data_map_.find( res_pair ) == pair_to_data_map_.end() ) {
					pair_to_data_map_[ res_pair ] = std::map< Size, Real >();
				}
				pair_to_data_map_[ res_pair ][ edge.first ] = weight_per_pair;
			}
		}
	}

	// Next pull custom coordinates from command line
	auto sl_data = DEERIO().pull_coords();
	if ( sl_data.size() == 0 ) return;
	for ( auto const & sl_weight : sl_data ) {
		spin_labels_.push_back( sl_weight.first );
		sl_weights_.push_back( sl_weight.second );
	}
}

/// @brief Returns the number of edges corresponding to the pair of residues
std::map< Size, Real >
DEERDataCache::edge(
	Size const & rsd1,
	Size const & rsd2
) const {
	std::pair< Size, Size > respair = std::make_pair( std::min( rsd1, rsd2 ), std::max( rsd1, rsd2 ) );
	if ( pair_to_data_map_.find( respair ) == pair_to_data_map_.end() ) {
		return std::map< Size, Real >();
	} else return pair_to_data_map_.at( respair );
}

/// @brief Returns list of residues and appropriate spin labels to compute
std::set< std::pair< Size, std::string > >
DEERDataCache::labeled_residues() const {
	return labeled_residues_;
}

/// @brief Stores the normalized coordinates for residues with CUSTOM spin labels
utility::vector1< EPRSpinLabel >
DEERDataCache::labels() const {
	return spin_labels_;
}

/// @brief Stores the normalized coordinates for residues with CUSTOM spin labels
void
DEERDataCache::set_labels(
	utility::vector1< EPRSpinLabel > const & labels
) {
	spin_labels_ = labels;
}

/// @brief Returns the weights assigned to CUSTOM spin label coordinates
utility::vector1< Real > const &
DEERDataCache::sl_weights() const {
	return sl_weights_;
}

/// @brief Stores the weights assigned to CUSTOM spin label coordinates
void
DEERDataCache::set_sl_weights(
	utility::vector1< Real > const & weights
) {
	sl_weights_ = weights;
}

bool
DEERDataCache::has_f1_force(
	Size const & res
) const {
	return f1_forces_.find( res ) != f1_forces_.end();
}

bool
DEERDataCache::has_f2_force(
	Size const & res
) const {
	return f2_forces_.find( res ) != f2_forces_.end();
}

/// @brief Stores the F1 force applied to a residue, used for gradient minimzation
void
DEERDataCache::f1_force(
	Size const & res,
	numeric::xyzVector< Real > const & vec
) {
	f1_forces_[ res ] = vec;
}

/// @brief Stores the F2 force applied to a residue, used for gradient minimzation
void
DEERDataCache::f2_force(
	Size const & res,
	numeric::xyzVector< Real > const & vec
) {
	f2_forces_[ res ] = vec;
}

/// @brief Returns the F1 force applied to a residue, used for gradient minimzation
numeric::xyzVector< Real > const &
DEERDataCache::f1_force(
	Size const & res
) const {
	return f1_forces_.at( res );
}

/// @brief Returns the F2 force applied to a residue, used for gradient minimzation
numeric::xyzVector< Real > const &
DEERDataCache::f2_force(
	Size const & res
) const {
	return f2_forces_.at( res );
}

/// @brief Sets the relative weight assigned to the pose storing this object. Used for multi-pose fitting
void
DEERDataCache::set_ensemble_weight(
	Real const & input
) {
	ensemble_weight_ = input;
}

/// @brief Returns the relative weight assigned to the pose storing this object. Used for multi-pose fitting
Real
DEERDataCache::ensemble_weight() const {
	return ensemble_weight_;
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
