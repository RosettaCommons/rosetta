// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PRESingle.cc
/// @brief   Implementation of class PRESingle
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pre/PRESingle.hh>

// Package headers
#include <core/io/nmr/AtomSelection.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.json.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector.hh>

// C++ headers
#include <iostream>
#include <iomanip>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

static basic::Tracer TR( "core.scoring.nmr.pre.PRESingle" );

// @brief default constructor
PRESingle::PRESingle() :
	protein_spins_(),
	pre_exp_(0),
	pre_err_(0),
	pre_calc_(0),
	weight_(1.0),
	atom_derivatives_()
{}

/// @brief constructor with arguments
PRESingle::PRESingle(
	utility::vector1< core::io::nmr::AtomSelection > const & spins,
	pose::Pose const & pose,
	Real const pre_exp,
	Real const pre_err
) :
	protein_spins_(),
	pre_exp_(pre_exp),
	pre_err_(pre_err),
	pre_calc_(0),
	weight_(1.0),
	atom_derivatives_()
{
	runtime_assert_msg(spins.size() > 0, "ERROR during creation of PRESingle object. Vector of spin atom selections must not be empty.");

	using utility::to_string;

	pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	for ( Size i = 1, i_end = spins.size(); i <= i_end; ++i ) {
		Size rsd = spins[i].get_rsd();
		std::string atom = spins[i].get_atom();
		char const & chain = spins[i].get_chain();

		if ( pdbinfo ) {
			TR.Trace << "Converting AtomSelection \"" + to_string(rsd) + " " + atom + " " + to_string(chain) + "\" in PRE datafile from pdb to pose numbering." << std::endl;
			rsd = pdbinfo->pdb2pose(chain, rsd);
			if ( rsd == 0 ) { // Residue not found
				utility_exit_with_message( "ERROR: Cannot convert AtomSelection " + to_string(spins[i].get_rsd()) + " " + atom + " " + to_string(chain) + " from PRE datafile into pose numbering. Residue number was not found." );
			}
		} else { // Assume pose numbering instead
			TR.Warning << "Cannot convert residue numbers from pdb to pose numbering. No PDBInfo object. Make sure to provide residues in pose numbering in PRE datafile instead." << std::endl;
		}

		if ( i == 1 ) { pre_type_ = pre_spin_type_from_atom_name(spins[i]); }
		if ( pre_type_ != pre_spin_type_from_atom_name(spins[i]) ) {
			utility_exit_with_message("ERROR in creating PRESingle object. PRE atoms are not of the same type.");
		}
		// Look up pseudoprotons and convert to AtomIDs
		utility::vector1< id::AtomID >protein_spins = lookup_pseudoprotons(rsd, atom, pose);
		if ( protein_spins.empty() ) {
			utility_exit_with_message("Could not convert NMR AtomSelection into Rosetta AtomIDs.");
		}

		for ( Size j = 1; j <= protein_spins.size(); ++j ) {
			protein_spins_.push_back(protein_spins[j]);
		}
	}
	atom_derivatives_.resize(protein_spins_.size());
}

/// @brief copy constructor
PRESingle::PRESingle(PRESingle const & other) :
	protein_spins_(other.protein_spins_),
	pre_exp_(other.pre_exp_),
	pre_err_(other.pre_err_),
	pre_calc_(other.pre_calc_),
	weight_(other.weight_),
	atom_derivatives_(other.atom_derivatives_),
	pre_type_(other.pre_type_)
{}

/// @brief assignment operator
PRESingle&
PRESingle::operator=(PRESingle const & rhs) {
	if ( this != & rhs ) {
		protein_spins_ = rhs.protein_spins_;
		pre_exp_ = rhs.pre_exp_;
		pre_err_ = rhs.pre_err_;
		pre_calc_ = rhs.pre_calc_;
		weight_ = rhs.weight_;
		atom_derivatives_ = rhs.atom_derivatives_;
		pre_type_ = rhs.pre_type_;
	}
	return *this;
}

/// @brief destructor
PRESingle::~PRESingle() {}

void
PRESingle::show(std::ostream & tracer) const {
	// Store old iostream manipulator flags
	std::ios oldState(nullptr);
	oldState.copyfmt(tracer);
	tracer << "( ";

	for ( Size i = 1; i <= protein_spins_.size()-1; ++i ) {
		tracer << std::setw(4) << protein_spins_[i].rsd() << " " << std::setw(4) << protein_spins_[i].atomno() << " or ";
	}
	tracer << std::setw(4) << protein_spins_.back().rsd() << " " << std::setw(4) << protein_spins_.back().atomno() << " ) ";
	tracer << std::setw(6) << std::fixed << std::setprecision(3) << pre_exp_ << " +/- " << pre_err_ << std::endl;
	tracer.copyfmt(oldState);
}

bool operator==(PRESingle const & lhs, PRESingle const & rhs) {
	if ( lhs.protein_spins_.size() != rhs.protein_spins_.size() ) {
		return false;
	}
	for ( Size i = 1, i_end = lhs.protein_spins_.size(); i <= i_end; ++i ) {
		if ( lhs.protein_spins_[i] != rhs.protein_spins_[i] ) {
			return false;
		}
	}
	return true;
}

bool operator!=(PRESingle const & lhs, PRESingle const & rhs) {
	return !( lhs == rhs );
}

/// @brief serialize a PRESingle object to a json_spirit object
utility::json_spirit::Value
PRESingle::serialize() const {
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	std::vector<Value> pre_spins_vector;
	for ( Size i = 1; i <= protein_spins_.size(); ++i ) {
		pre_spins_vector.push_back(protein_spins_[i].serialize());
	}
	Pair pre_spins_record("pre_spins", pre_spins_vector);
	Pair pre_exp("pre_exp", Value(pre_exp_));
	Pair pre_err("pre_err", Value(pre_err_));
	Pair pre_calc("pre_calc", Value(pre_calc_));
	Pair weight("pre_weight", Value(weight_));
	std::vector<Value> atom_deriv_vector;
	for ( Size i = 1; i <= atom_derivatives_.size(); ++i ) {
		atom_deriv_vector.push_back(numeric::serialize(atom_derivatives_[i]));
	}
	Pair atom_deriv_record("derivatives", atom_deriv_vector);

	return Value(utility::tools::make_vector(pre_spins_record,pre_exp, pre_err, pre_calc, weight,atom_deriv_record));
}

/// @brief deserialize a json_spirit object to a PRESingle object
void
PRESingle::deserialize(utility::json_spirit::mObject data) {
	utility::json_spirit::mArray pre_spins_data(data["pre_spins"].get_array());
	protein_spins_.resize(pre_spins_data.size());
	Size i(1);
	for ( utility::json_spirit::mArray::iterator it = pre_spins_data.begin(); it != pre_spins_data.end(); ++it ) {
		core::id::AtomID id;
		id.deserialize(it->get_obj());
		protein_spins_[i++] = id;
	}

	pre_exp_  = data["pre_exp"].get_real();
	pre_err_  = data["pre_err"].get_real();
	pre_calc_ = data["pre_calc"].get_real();
	weight_   = data["pre_weight"].get_real();

	utility::json_spirit::mArray atom_deriv_data(data["derivatives"].get_array());
	atom_derivatives_.resize(atom_deriv_data.size());
	Size j(1);
	for ( utility::json_spirit::mArray::iterator it = atom_deriv_data.begin(); it != atom_deriv_data.end(); ++it ) {
		atom_derivatives_[j++] = numeric::deserialize<Real>(it->get_array());
	}
}

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core


