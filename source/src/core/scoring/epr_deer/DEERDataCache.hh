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

#ifndef INCLUDED_core_scoring_epr_deer_DEERDataCache_hh
#define INCLUDED_core_scoring_epr_deer_DEERDataCache_hh

// Unit headers
#include <core/scoring/epr_deer/DEERDataCache.fwd.hh>
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>

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

	/// @brief Returns non-const data for a given ID
	DEERDataOP &
	operator[](
		Size const & i
	);

	/// @brief Returns const data for a given ID
	DEERDataOP const &
	at(
		Size const & i
	) const;

	/// @brief Adds data to vector
	void
	append(
		DEERDataOP const & data
	);

	/// @brief Returns number of DEERData objects stored here
	Size
	size() const;

	/// @brief Fetches are parses data from command line using DEERIO object
	void
	fetch_and_organize_data(
		pose::Pose & pose
	);

	/// @brief Returns the number of edges corresponding to the pair of residues
	std::map< Size, Real >
	edge(
		Size const & rsd1,
		Size const & rsd2
	) const;

	/// @brief Returns list of residues and appropriate spin labels to compute
	std::set< std::pair< Size, std::string > >
	labeled_residues() const;

	/// @brief Stores the normalized coordinates for residues with CUSTOM spin labels
	void
	set_labels(
		utility::vector1< EPRSpinLabel > const & labels
	);

	/// @brief Returns the weights assigned to CUSTOM spin label coordinates
	utility::vector1< Real > const &
	sl_weights() const;

	/// @brief Stores the weights assigned to CUSTOM spin label coordinates
	void
	set_sl_weights(
		utility::vector1< Real > const & weights
	);

	/// @brief Returns the normalized coordinates for residues with CUSTOM spin labels
	utility::vector1< EPRSpinLabel >
	labels() const;

	/// @brief Stores the F1 force applied to a residue, used for gradient minimzation
	void
	f1_force(
		Size const & res,
		numeric::xyzVector< Real > const & vec
	);

	/// @brief Stores the F2 force applied to a residue, used for gradient minimzation
	void
	f2_force(
		Size const & res,
		numeric::xyzVector< Real > const & vec
	);

	/// @brief Checks if residues has an F1 force, used for gradient minimzation
	bool
	has_f1_force(
		Size const & res
	) const;

	/// @brief Checks if residues has an F1 force, used for gradient minimzation
	bool
	has_f2_force(
		Size const & res
	) const;

	/// @brief Returns the F1 force applied to a residue, used for gradient minimzation
	numeric::xyzVector< Real > const &
	f1_force(
		Size const & res
	) const;

	/// @brief Returns the F2 force applied to a residue, used for gradient minimzation
	numeric::xyzVector< Real > const &
	f2_force(
		Size const & res
	) const;

	/// @brief Sets the relative weight assigned to the pose storing this object. Used for multi-pose fitting
	void
	set_ensemble_weight(
		Real const & input
	);

	/// @brief Returns the relative weight assigned to the pose storing this object. Used for multi-pose fitting
	Real
	ensemble_weight() const;

private:

	std::map< Size, DEERDataOP > data_;
	std::map< Size, numeric::xyzVector< Real > > f1_forces_;
	std::map< Size, numeric::xyzVector< Real > > f2_forces_;
	//std::map< Size, Real > rmsd_values_;
	std::set< std::pair< Size, std::string > > labeled_residues_;
	std::map< std::pair< Size, Size >, std::map< Size, Real > > pair_to_data_map_;
	utility::vector1< EPRSpinLabel > spin_labels_;
	utility::vector1< Real > sl_weights_;
	Real ensemble_weight_ = 0.0;

	//std::map< Size, utility::vector1< PseudoElectron > > electrons_;
	//bool using_custom_electrons_;
};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
