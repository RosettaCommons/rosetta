// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingle.hh
/// @brief   class that stores and handles data related to one single PCS
/// @details last Modified: 06/21/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pcs_PCSSingle_HH
#define INCLUDED_core_scoring_nmr_pcs_PCSSingle_HH

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingle.fwd.hh>

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

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

class PCSSingle {

public: // Methods

	/// @brief default constructor
	PCSSingle();

	/// @brief constructor with arguments
	PCSSingle(
		utility::vector1< core::io::nmr::AtomSelection > const & spin,
		pose::Pose const & pose,
		Real const pcs_exp,
		Real const pcs_err
	);

	/// @brief copy constructor
	PCSSingle(PCSSingle const & other);

	/// @brief assignment operator
	PCSSingle&
	operator=(PCSSingle const & rhs);

	/// @brief destructor
	~PCSSingle();

	/// @brief serialize a PCSSingle object to a json_spirit object
	utility::json_spirit::Value serialize() const;

	/// @brief deserialize a json_spirit object to a PCSSingle object
	void deserialize(utility::json_spirit::mObject data);

	// Getters
	utility::vector1< id::AtomID > const & get_protein_spins() const { return protein_spins_; }
	Real get_pcs_exp() const { return pcs_exp_; }
	Real get_pcs_err() const { return pcs_err_; }
	Real get_pcs_calc() const { return pcs_calc_; }
	Real get_weight() const { return weight_; }
	utility::vector1< Vector > const & get_atom_derivatives() const { return atom_derivatives_; }

	// Setters
	void set_pcs_exp(Real pcs) { pcs_exp_ = pcs; }
	void set_pcs_err(Real err) { pcs_err_ = err; }
	void set_pcs_calc(Real calc) { pcs_calc_ = calc; }
	void set_weight(Real weight) { weight_ = weight; }
	void set_atom_derivatives(Size idx, Real fdx, Real fdy, Real fdz )
	{
		// We don't set the derivatives to 0.0 and/or use operator += here
		// The accumulation happens in the energy method.
		// If we have for the same spin a second (or third ...) PCSSingle object
		// (e.g. if we have duplicated the original PCS data for symmetric subunits) the
		// derivative from all PCSSingle objects (i.e. from the other symmetric subunits)
		// gets summed up in the energy method.
		atom_derivatives_[idx].x() = fdx;
		atom_derivatives_[idx].y() = fdy;
		atom_derivatives_[idx].z() = fdz;
	}

	void show(std::ostream & TR) const;

	friend bool operator==(PCSSingle const & lhs, PCSSingle const & rhs);
	friend bool operator!=(PCSSingle const & lhs, PCSSingle const & rhs);

private: // Data

	utility::vector1< id::AtomID > protein_spins_;
	Real pcs_exp_;
	Real pcs_err_;
	Real pcs_calc_;
	Real weight_;
	utility::vector1< Vector > atom_derivatives_;

};

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pcs_PCSSingle_HH
