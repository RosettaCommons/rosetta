// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCSingle.hh
/// @brief   class that stores and handles data related to one single RDC
/// @details last Modified: 06/21/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_rdc_RDCSingle_HH
#define INCLUDED_core_scoring_nmr_rdc_RDCSingle_HH

// Unit headers
#include <core/scoring/nmr/rdc/RDCSingle.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Package headers
#include <core/io/nmr/AtomSelection.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/json_spirit/json_spirit_value.h>

// C++ headers
#include <iostream>
#include <string>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

class RDCSingle {

public: // Methods

	/// @brief default constructor
	RDCSingle();

	/// @brief constructor with arguments
	RDCSingle(
		utility::vector1< std::pair< core::io::nmr::AtomSelection, core::io::nmr::AtomSelection > > const & spinsAB,
		pose::Pose const & pose,
		Real const rdc_exp,
		Real const rdc_err
	);

	/// @brief copy constructor
	RDCSingle(RDCSingle const & other);

	/// @brief assignment operator
	RDCSingle&
	operator=(RDCSingle const & rhs);

	/// @brief destructor
	~RDCSingle();

	/// @brief serialize a RDCSingle object to a json_spirit object
	utility::json_spirit::Value serialize() const;

	/// @brief deserialize a json_spirit object to a RDCSingle object
	void deserialize(utility::json_spirit::mObject data);

	// Getters
	utility::vector1< std::pair< id::AtomID, id::AtomID > > const & get_spinsAB() const { return spinsAB_; }
	Real get_rdc_exp() const { return rdc_exp_; }
	Real get_rdc_err() const { return rdc_err_; }
	Real get_rdc_calc() const { return rdc_calc_; }
	Real get_weight() const { return weight_; }
	utility::vector1< std::pair< Vector, Vector > > const & get_atom_derivatives() const { return atom_derivatives_; }
	RDC_TYPE get_rdc_type() const { return rdc_type_; }

	// Setters
	void set_rdc_exp(Real rdc) { rdc_exp_ = rdc; }
	void set_rdc_err(Real err) { rdc_err_ = err; }
	void set_rdc_calc(Real calc) { rdc_calc_ = calc; }
	void set_weight(Real weight) { weight_ = weight; }
	void set_atom_derivatives(Size idx, Real fdx, Real fdy, Real fdz)
	{
		// We don't set the derivatives to 0.0 and/or use operator += here
		// The accumulation happens in the energy method.
		// If we have for the same spin a second (or third ...) RDCSingle object
		// (e.g. if we have duplicated the original RDC data for symmetric subunits) the
		// derivative from all RDCSingle objects (i.e. from the other symmetric subunits)
		// gets summed up in the energy method.
		atom_derivatives_[idx].first.x() = fdx;
		atom_derivatives_[idx].first.y() = fdy;
		atom_derivatives_[idx].first.z() = fdz;
		atom_derivatives_[idx].second.x() = -fdx;
		atom_derivatives_[idx].second.y() = -fdy;
		atom_derivatives_[idx].second.z() = -fdz;
	}

	void show(std::ostream & TR) const;

	friend bool operator==(RDCSingle const & lhs, RDCSingle const & rhs);
	friend bool operator!=(RDCSingle const & lhs, RDCSingle const & rhs);

private: // Data

	// Vector of ambiguous AB spin pairs for which an average RDC is observed.
	// This can include spins from equivalent subunits or degenerate spin pairs (e.g. Ala CH3).
	// In case of Ala CH3 the vector would contain three pairs: [ (CB, 1HB), (CB, 2HB), (CB, 3HB) ]
	// Or in case of CA-HA RDCs for Gly. There are two bond vectors that contribute to the RDC CA-1HA and CA-2HA.
	utility::vector1< std::pair< id::AtomID, id::AtomID > > spinsAB_;
	Real rdc_exp_;
	Real rdc_err_;
	Real rdc_calc_;
	Real weight_;
	utility::vector1< std::pair< Vector, Vector > > atom_derivatives_;
	RDC_TYPE rdc_type_;

};

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_rdc_RDCSingle_HH
