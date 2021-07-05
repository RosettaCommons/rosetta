// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/DEERDataCache.hh
/// @brief  Contains data for DEER energy method, stored as DEERData objects
/// @details To prevent the energy method from reading from the command line every scoring
///      round and parsing the input file, this method stores a list of DEER decay data.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_DEERDataCache_hh
#define INCLUDED_core_scoring_epr_deer_DEERDataCache_hh

// Unit headers
#include <core/scoring/epr_deer/DEERDataCache.fwd.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/metrics/DEERData.fwd.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>
#include <map>

#include <set> // AUTO IWYU For set

namespace core {
namespace scoring {
namespace epr_deer {

class DEERDataCache : public basic::datacache::CacheableData {
public:

	/// @brief Constructor
	DEERDataCache() = default;

	/// @brief Copy constructor
	DEERDataCache( DEERDataCache const & ) = default;

	/// @brief Destructor
	~DEERDataCache() override;

	/// @brief Copy function, overrides parent class
	basic::datacache::CacheableDataOP
	clone() const override;

	/// @brief Returns non-const data for a given ID (or nullptr if empty)
	/// @param  i: Index of data object ()
	/// @return Requested DEERDataOP object (or nullptr)
	metrics::DEERDataOP &
	operator[](
		Size const & i
	);

	/// @brief Returns const data for a given ID (or nullptr if empty)
	/// @param  i: Index of data object ()
	/// @return Requested DEERDataOP object (or nullptr)
	metrics::DEERDataOP const &
	at(
		Size const & i
	) const;

	/// @brief Adds DEERDataOP object to DEERDataCache at earliest spot
	/// @param data: metrics::DEERDataOP object to insert
	void
	append(
		metrics::DEERDataOP const & data
	);

	/// @brief Adds DEERDataOP object to DEERDataCache at predefined spot
	/// @param data: Object to insert
	/// @param i: Position to insert object
	void
	append(
		metrics::DEERDataOP const & data,
		Size const & i
	);

	/// @brief Adds DEERDataOP object to DEERDataCache at last spot
	/// @param data: metrics::DEERDataOP object to insert
	void
	push_back(
		metrics::DEERDataOP const & data
	);

	/// @brief Returns number of DEERDataOPs stored here
	/// @return Size object ofnumber of DEERDataOPs in map
	/// @detail Note: output could be less than data_.rbegin()->first)
	Size
	size() const;

	/// @brief Returns set of keys in the data map
	/// @return Set object containing all keys in map
	std::set< Size >
	indices() const ;

	/// @brief Fetches are parses data from command line using DEERIO object
	/// @param pose: Pose object from which data will be fetched
	void
	fetch_and_organize_data(
		pose::Pose const & pose
	);

	/// @brief Initializer function following reading of data
	void
	set_up_residues();

	/// @brief Returns number of edges corresponding to the pair of residues
	/// @param rsd1: identity of first residue (AA sequence)
	/// @param rsd2: identity of second residue (AA sequence)
	/// @return Map with edges and weights for residue pair
	std::map< Size, Real >
	edge(
		Size const & rsd1,
		Size const & rsd2
	) const;

	/// @brief Returns set of residues and appropriate spin labels to compute
	/// @return Residues in data
	//std::set< std::pair< Size, std::string > >
	//labeled_residues() const;

	/// @brief Returns vector of residues and appropriate spin labels to compute
	/// @return Residues in data
	utility::vector1< std::pair< Size, std::string > >
	labeled_residues() const;

	/// @brief Returns list of residues spin labels to compute
	/// @return Labels required by data
	std::set< std::string >
	label_types() const;

	/// @brief Stores the coordinates for residues with CUSTOM spin labels
	/// @return utility::vector1< EPRSpinLabel > of custom SL rotamers in data
	utility::vector1< EPRSpinLabel >
	labels() const;

	/// @brief Stores the normalized coordinates for residues with CUSTOM spin labels
	/// @param labels: Custom rotamers across model
	void
	set_labels(
		utility::vector1< EPRSpinLabel > const & labels
	);

	/// @brief Returns the weights assigned to CUSTOM spin label coordinates
	/// @return Vector of weights
	utility::vector1< Real > const &
	sl_weights() const;

	/// @brief Stores the weights assigned to CUSTOM spin label coordinates
	/// @param weights: Vector of weights for custom SLs
	void
	set_sl_weights(
		utility::vector1< Real > const & weights
	);

	/// @brief Returns if F1 force is stored for residue
	/// @param res: Residue ID of interest
	/// @return True or false if residue is stored in F1 force map
	bool
	has_f1_force(
		Size const & res
	) const;

	/// @brief Checks if residues has an F1 force, used for gradient minimzation
	bool
	has_f2_force(
		Size const & res
	) const;

	/// @brief Sets F1 force stored for residue
	/// @param res: Residue ID of interest
	/// @param vec: Vector of F1 force
	void
	f1_force(
		Size const & res,
		numeric::xyzVector< Real > const & vec
	);

	/// @brief Sets F2 force stored for residue
	/// @param res: Residue ID of interest
	/// @param vec: Vector of F2 force
	void
	f2_force(
		Size const & res,
		numeric::xyzVector< Real > const & vec
	);

	/// @brief Returns the F1 force applied to a residue
	/// @param res: Residue ID of interest
	/// @return F1 force
	numeric::xyzVector< Real > const &
	f1_force(
		Size const & res
	) const;

	/// @brief Returns the F2 force applied to a residue
	/// @param res: Residue ID of interest
	/// @return F2 force
	numeric::xyzVector< Real > const &
	f2_force(
		Size const & res
	) const;

	/// @brief Removes data at specific spot and reparameterizes
	/// @param i: Key of data in map
	void
	delete_data(
		core::Size const & i
	);

	/// @brief Check and return if data in spot is occupied
	/// @param i: Key to check
	/// @return Whether key is being used
	bool
	has_data(
		core::Size const & i
	);

private:

	/// @brief Data container
	std::map< Size, metrics::DEERDataOP > data_;

	/// @brief Container for F1 forces
	std::map< Size, numeric::xyzVector< Real > > f1_forces_;

	/// @brief Container for F2 forces
	std::map< Size, numeric::xyzVector< Real > > f2_forces_;

	/// @brief Set of all labeled residues (string is SL type)
	std::set< std::pair< Size, std::string > > labeled_residues_;

	/// @brief Map of IDs and weights corresponding to residue pairs
	std::map< std::pair< Size, Size >, std::map< Size, Real > > pair_data_;

	/// @brief List of CUSTOM EPRSpin labels (obtained by multilateration)
	utility::vector1< EPRSpinLabel > spin_labels_;

	/// @brief Weights for CUSTOM labels (obtained by multilateration)
	utility::vector1< Real > sl_weights_;

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
