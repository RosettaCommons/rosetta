// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PRESingle.hh
/// @brief   class that stores and handles data related to one single PRE value
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pre_PRESingle_HH
#define INCLUDED_core_scoring_nmr_pre_PRESingle_HH

// Unit headers
#include <core/scoring/nmr/pre/PRESingle.fwd.hh>

// Package headers
#include <core/io/nmr/AtomSelection.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/json_spirit/json_spirit_value.h>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

class PRESingle {

public: // Methods

	/// @brief default constructor
	PRESingle();

	/// @brief constructor with arguments
	PRESingle(
		utility::vector1< core::io::nmr::AtomSelection > const & spins,
		pose::Pose const & pose,
		Real const pre_exp,
		Real const pre_err
	);

	/// @brief copy constructor
	PRESingle(PRESingle const & other);

	/// @brief assignment operator
	PRESingle &
	operator=(PRESingle const & rhs);

	/// @brief destructor
	~PRESingle();

	/// @brief serialize a PRESingle object to a json_spirit object
	utility::json_spirit::Value serialize() const;

	/// @brief deserialize a json_spirit object to a PRESingle object
	void deserialize(utility::json_spirit::mObject data);

	// Getters
	utility::vector1< id::AtomID > const & get_protein_spins() const { return protein_spins_; }
	Real get_pre_exp() const { return pre_exp_; }
	Real get_pre_err() const { return pre_err_; }
	Real get_pre_calc() const { return pre_calc_; }
	Real get_weight() const { return weight_; }
	utility::vector1< Vector > const & get_atom_derivatives() const { return atom_derivatives_; }
	PRE_SPIN_TYPE get_pre_spin_type() const { return pre_type_; }

	// Setters
	void set_pre_exp(Real pre) { pre_exp_ = pre; }
	void set_pre_err(Real err) { pre_err_ = err; }
	void set_pre_calc(Real calc) { pre_calc_ = calc; }
	void set_weight(Real weight) { weight_ = weight; }
	void set_atom_derivatives(Size index, Real fdx, Real fdy, Real fdz)
	{
		// We don't set the derivatives to 0.0 and/or use operator += here
		// The accumulation happens in the energy method.
		// If we have for the same spin a second (or third ...) PRESingle object
		// (e.g. if we have duplicated the original PRE data for symmetric subunits) the
		// derivative from all PRESingle objects (i.e. from the other symmetric subunits)
		// gets summed up in the energy method.
		atom_derivatives_[index].x() = fdx;
		atom_derivatives_[index].y() = fdy;
		atom_derivatives_[index].z() = fdz;
	}

	void show(std::ostream & TR) const;

	friend bool operator==(PRESingle const & lhs, PRESingle const & rhs);
	friend bool operator!=(PRESingle const & lhs, PRESingle const & rhs);

private: // Data

	utility::vector1<id::AtomID> protein_spins_;
	Real pre_exp_;
	Real pre_err_;
	Real pre_calc_;
	Real weight_;
	utility::vector1< Vector > atom_derivatives_;
	PRE_SPIN_TYPE pre_type_;

};

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pre_PRESpinlabel_HH
