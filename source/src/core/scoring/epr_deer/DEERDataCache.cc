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
#include <core/scoring/epr_deer/DEERIO.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/metrics/DEERData.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/types.hh>

// Basic headers
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <set>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.metrics.DEERDataCache" );

/// @brief Destructor
DEERDataCache::~DEERDataCache() {}

/// @brief Copy function, overrides parent class
basic::datacache::CacheableDataOP
DEERDataCache::clone() const {
	return basic::datacache::CacheableDataOP( new DEERDataCache( *this ) );
}

/// @brief Returns non-const data for a given ID (or nullptr if empty)
/// @param  i: Index of data object ()
/// @return Requested DEERDataOP object (or nullptr)
metrics::DEERDataOP &
DEERDataCache::operator[](
	Size const & i
) {
	if ( data_.find( i ) == data_.end() ) {
		throw CREATE_EXCEPTION( utility::excn::KeyError, "Data with key "
			+ std::to_string( i ) + " not found in DEERDataCache!" );
	} else {
		return data_[ i ];
	}
}

/// @brief Returns const data for a given ID (or nullptr if empty)
/// @param  i: Index of data object ()
/// @return Requested DEERDataOP object (or nullptr)
metrics::DEERDataOP const &
DEERDataCache::at(
	Size const & i
) const {
	if ( data_.find( i ) == data_.end() ) {
		throw CREATE_EXCEPTION( utility::excn::KeyError, "Data with key "
			+ std::to_string( i ) + " not found in DEERDataCache!" );
	} else {
		return data_.at( i );
	}
}

/// @brief Adds DEERDataOP object to DEERDataCache at earliest spot
/// @param data: metrics::DEERDataOP object to insert
/// @detail Not threadsafe
void
DEERDataCache::append(
	metrics::DEERDataOP const & data
) {

	// Failsafe in case there are no objects in vector
	if ( data_.size() == 0 ) {
		data_[ 1 ] = data;
	} else {
		for ( Size i = 1; i <= size() + 1; ++i ) {
			if ( data_.find( i ) == data_.end() ) {
				data_[ i ] = data;
				return;
			}
		}
	}
}

/// @brief Adds DEERDataOP object to DEERDataCache at predefined spot
/// @param data: Object to insert
/// @param i: Position to insert object
/// @detail Not threadsafe
void
DEERDataCache::append(
	metrics::DEERDataOP const & data,
	Size const & i
) {
	if ( data_.find( i ) != data_.end() ) {
		TR.Warning << "Overwriting data at pos " << i << std::endl;
		TR.Warning << "Existing data has been discarded!" << std::endl;
	}
	data_[ i ] = data;
}

/// @brief Adds DEERDataOP object to DEERDataCache at last spot
/// @param data: metrics::DEERDataOP object to insert
/// @detail Not threadsafe
void
DEERDataCache::push_back(
	metrics::DEERDataOP const & data
) {

	// Failsafe in case there are no objects in vector
	if ( data_.size() == 0 ) {
		data_[ 1 ] = data;
	} else {
		data_[ data_.rbegin()->first + 1 ] = data;
	}
}

/// @brief Returns number of DEERDataOPs stored here
/// @return Size object ofnumber of DEERDataOPs in map
/// @detail Note: output could be less than data_.rbegin()->first)
Size
DEERDataCache::size() const {
	return data_.size();
}

/// @brief Returns set of keys in the data map
/// @return Set object containing all keys in map
std::set< Size >
DEERDataCache::indices() const {
	std::set< Size > output;
	for ( auto const & data : data_ ) {
		output.insert( data.first );
	}
	return output;
}

/// @brief Fetches are parses data from command line using DEERIO object
/// @param pose: Pose object from which data will be fetched
void
DEERDataCache::fetch_and_organize_data(
	pose::Pose const & pose
) {

	// Data IO from file
	data_ = DEERIO().generate_data( pose );
	set_up_residues();

	// In case epr_deer:custom_coords is used
	auto const sl_data = DEERIO().pull_coords( pose );
	if ( sl_data.size() ) {
		for ( auto const & sl_weight : sl_data ) {
			spin_labels_.push_back( sl_weight.first );
			sl_weights_.push_back( sl_weight.second );
		}
	}
}

/// @brief Initializer function following reading of data
void
DEERDataCache::set_up_residues() {

	// Clear existing data
	labeled_residues_.clear();
	pair_data_.clear();

	// Iterate across sets for redistribution of data
	for ( auto const & edge : data_ ) {
		auto const & ids = edge.second->residues();
		for ( auto const & res : ids ) {
			labeled_residues_.insert( res );
		}

		// Iterate across pairs of residues in that dataset
		for ( Size i = 1; i < ids.size(); ++i ) {
			for ( Size j = i + 1; j <= ids.size(); ++j ) {

				// Pair ID is ( first, last )
				auto const pair = std::make_pair(
					std::min( ids.at( i ).first, ids.at( j ).first ),
					std::max( ids.at( i ).first, ids.at( j ).first ) );

				// Add pair (and weight) to appropriate map
				if ( pair_data_.find( pair ) == pair_data_.end() ) {
					pair_data_[ pair ] = std::map< Size, Real >();
				}
				auto const w = ids.size() * ( ids.size() - 1 ) / 2.0;
				pair_data_[ pair ][ edge.first ] = w;
			}
		}
	}
}

/// @brief Returns number of edges corresponding to the pair of residues
/// @param rsd1: identity of first residue (AA sequence)
/// @param rsd2: identity of second residue (AA sequence)
/// @return Map with edges and weights for residue pair
std::map< Size, Real >
DEERDataCache::edge(
	Size const & rsd1,
	Size const & rsd2
) const {

	// Pair ID is ( first, last )
	std::pair< Size, Size > const respair = std::make_pair(
		std::min( rsd1, rsd2 ),
		std::max( rsd1, rsd2 ) );

	// Return the map (if found) or an empty map (if not found)
	if ( pair_data_.find( respair ) == pair_data_.end() ) {
		return std::map< Size, Real >();
	} else {
		return pair_data_.at( respair );
	}
}

/// @brief Returns set of residues and appropriate spin labels to compute
/// @return Residues in data
//std::set< std::pair< Size, std::string > >
//DEERDataCache::labeled_residues() const {
// return labeled_residues_;
//}

/// @brief Returns vector of residues and appropriate spin labels to compute
/// @return Residues in data
utility::vector1< std::pair< Size, std::string > >
DEERDataCache::labeled_residues() const {
	utility::vector1< std::pair< Size, std::string > > output;
	for ( auto const & res_sl : labeled_residues_ ) {
		output.push_back( res_sl );
	}
	return output;
}

/// @brief Returns list of residues spin labels to compute
/// @return Labels required by data
std::set< std::string >
DEERDataCache::label_types() const {
	std::set< std::string > output;
	for ( auto const & res_sl : labeled_residues_ ) {
		output.insert( res_sl.second );
	}
	return output;
}

/// @brief Stores the coordinates for residues with CUSTOM spin labels
/// @return utility::vector1< EPRSpinLabel > of custom SL rotamers in data
utility::vector1< EPRSpinLabel >
DEERDataCache::labels() const {
	return spin_labels_;
}

/// @brief Stores the normalized coordinates for residues with CUSTOM spin labels
/// @param labels: Custom rotamers across model
void
DEERDataCache::set_labels(
	utility::vector1< EPRSpinLabel > const & labels
) {
	spin_labels_ = labels;
}

/// @brief Returns the weights assigned to CUSTOM spin label coordinates
/// @return Vector of weights
utility::vector1< Real > const &
DEERDataCache::sl_weights() const {
	return sl_weights_;
}

/// @brief Stores the weights assigned to CUSTOM spin label coordinates
/// @param weights: Vector of weights for custom SLs
void
DEERDataCache::set_sl_weights(
	utility::vector1< Real > const & weights
) {
	sl_weights_ = weights;
}

/// @brief Returns if F1 force is stored for residue
/// @param res: Residue ID of interest
/// @return True or false if residue is stored in F1 force map
bool
DEERDataCache::has_f1_force(
	Size const & res
) const {
	return f1_forces_.find( res ) != f1_forces_.end();
}

/// @brief Returns if F2 force is stored for residue
/// @param res: Residue ID of interest
/// @return True or false if residue is stored in F2 force map
bool
DEERDataCache::has_f2_force(
	Size const & res
) const {
	return f2_forces_.find( res ) != f2_forces_.end();
}

/// @brief Sets F1 force stored for residue
/// @param res: Residue ID of interest
/// @param vec: Vector of F1 force
void
DEERDataCache::f1_force(
	Size const & res,
	numeric::xyzVector< Real > const & vec
) {
	f1_forces_[ res ] = vec;
}

/// @brief Sets F2 force stored for residue
/// @param res: Residue ID of interest
/// @param vec: Vector of F2 force
void
DEERDataCache::f2_force(
	Size const & res,
	numeric::xyzVector< Real > const & vec
) {
	f2_forces_[ res ] = vec;
}

/// @brief Returns the F1 force applied to a residue
/// @param res: Residue ID of interest
/// @return F1 force
numeric::xyzVector< Real > const &
DEERDataCache::f1_force(
	Size const & res
) const {
	return f1_forces_.at( res );
}

/// @brief Returns the F2 force applied to a residue
/// @param res: Residue ID of interest
/// @return F2 force
numeric::xyzVector< Real > const &
DEERDataCache::f2_force(
	Size const & res
) const {
	return f2_forces_.at( res );
}

/// @brief Removes data at specific spot and reparameterizes
/// @param i: Key of data in map
void
DEERDataCache::delete_data(
	core::Size const & i
) {
	data_.erase( i );
	set_up_residues();
}

/// @brief Check and return if data in spot is occupied
/// @param i: Key to check
/// @return Whether key is being used
bool
DEERDataCache::has_data(
	core::Size const & i
) {
	return data_.find( i ) != data_.end();
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
