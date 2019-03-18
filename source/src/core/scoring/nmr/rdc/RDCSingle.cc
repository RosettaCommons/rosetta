// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCSingle.cc
/// @brief   Implementation of class RDCSingle
/// @details last Modified: 06/21/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/rdc/RDCSingle.hh>

// Package headers
#include <core/io/nmr/AtomSelection.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>

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
namespace rdc {

static basic::Tracer TR( "core.scoring.nmr.rdc.RDCSingle" );

/// @brief default constructor
RDCSingle::RDCSingle() :
	spinsAB_(),
	rdc_exp_(0),
	rdc_err_(0),
	rdc_calc_(0),
	weight_(1.0),
	atom_derivatives_()
{}

/// @brief constructor with arguments
RDCSingle::RDCSingle(
	utility::vector1< std::pair< core::io::nmr::AtomSelection, core::io::nmr::AtomSelection > > const & spinsAB,
	pose::Pose const & pose,
	Real const rdc_exp,
	Real const rdc_err
) :
	spinsAB_(),
	rdc_exp_(rdc_exp),
	rdc_err_(rdc_err),
	rdc_calc_(0),
	weight_(1.0),
	atom_derivatives_()
{
	runtime_assert_msg(spinsAB.size() > 0, "ERROR during creation of RDCSingle object. Vector of spins A and B must not be empty.");

	using utility::to_string;

	pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	for ( Size i = 1, i_end = spinsAB.size(); i <= i_end; ++i ) {
		Size rsdA = spinsAB[i].first.get_rsd();
		Size rsdB = spinsAB[i].second.get_rsd();
		std::string atomA = spinsAB[i].first.get_atom();
		std::string atomB = spinsAB[i].second.get_atom();
		char const & chain_spinA = spinsAB[i].first.get_chain();
		char const & chain_spinB = spinsAB[i].second.get_chain();

		if ( pdbinfo ) {
			TR.Trace << "Converting AtomSelection \"(" + to_string(rsdA) + " " + atomA + " " + to_string(chain_spinA) + ") ("
				+ to_string(rsdB) + " " + atomB + " " + to_string(chain_spinB) + ")\" in RDC datafile from pdb to pose numbering." << std::endl;
			rsdA = pdbinfo->pdb2pose(chain_spinA, rsdA);
			rsdB = pdbinfo->pdb2pose(chain_spinB, rsdB);
			if ( rsdA == 0 || rsdB == 0 ) { // Residue not found
				utility_exit_with_message( "ERROR: Cannot convert AtomSelection \"(" + to_string(spinsAB[i].first.get_rsd()) + " " + atomA + " " + to_string(chain_spinA) + ") ("
					+ to_string(spinsAB[i].second.get_rsd()) + " " + atomB + " " + to_string(chain_spinB) + ")\" from RDC datafile into pose numbering. Residue number was not found." );
			}
		} else { // Assume pose numbering instead
			TR.Warning << "Cannot convert residue numbers from pdb to pose numbering. No PDBInfo object. Make sure to provide residues in pose numbering in RDC datafile instead." << std::endl;
		}

		// Look up pseudoprotons and convert to AtomIDs
		utility::vector1< id::AtomID > spinsA = lookup_pseudoprotons(rsdA, atomA, pose);
		utility::vector1< id::AtomID > spinsB = lookup_pseudoprotons(rsdB, atomB, pose);
		if ( spinsA.empty() || spinsB.empty() ) {
			utility_exit_with_message("Could not convert NMR AtomSelection into Rosetta AtomIDs.");
		}

		if ( spinsA.size() < spinsB.size() ) { // e.g. spinsA include CA and spinsB include 1HA and 2HA, then resize spinsA and copy AtomID for CA
			spinsA.resize(spinsB.size(), spinsA.front());
		} else if ( spinsA.size() > spinsB.size() ) { // the order of spins A and B in the input file might also be reversed (spinsA -> 1HA, 2HA, spinsB -> CA)
			spinsB.resize(spinsA.size(), spinsB.front());
		}

		for ( Size j = 1; j <= spinsA.size(); ++j ) {
			spinsAB_.push_back(std::make_pair(spinsA[j], spinsB[j]));
		}
	}

	atom_derivatives_.resize(spinsAB_.size());
	rdc_type_ = rdc_type_from_atom_names(spinsAB.front());
}

/// @brief copy constructor
RDCSingle::RDCSingle(RDCSingle const & other) :
	spinsAB_(other.spinsAB_),
	rdc_exp_(other.rdc_exp_),
	rdc_err_(other.rdc_err_),
	rdc_calc_(other.rdc_calc_),
	weight_(other.weight_),
	atom_derivatives_(other.atom_derivatives_),
	rdc_type_(other.rdc_type_)
{}

/// @brief assignment operator
RDCSingle &
RDCSingle::operator=(RDCSingle const & rhs) {
	if ( this != & rhs ) {
		spinsAB_ = rhs.spinsAB_;
		rdc_exp_ = rhs.rdc_exp_;
		rdc_err_ = rhs.rdc_err_;
		rdc_calc_ = rhs.rdc_calc_;
		weight_ = rhs.weight_;
		atom_derivatives_ = rhs.atom_derivatives_;
		rdc_type_ = rhs.rdc_type_;
	}
	return *this;
}

/// @brief destructor
RDCSingle::~RDCSingle() {}

void
RDCSingle::show(std::ostream & tracer) const {
	// Store old iostream manipulator flags
	std::ios oldState(nullptr);
	oldState.copyfmt(tracer);
	tracer << "( ";
	for ( Size i = 1; i <= spinsAB_.size() - 1; ++i ) {
		tracer << std::setw(4) << spinsAB_[i].first.rsd() << " " << std::setw(4) << spinsAB_[i].first.atomno() << " or ";
	}
	tracer << std::setw(4) << spinsAB_.back().first.rsd() << " " << std::setw(4) << spinsAB_.back().first.atomno() << " ) ( ";
	for ( Size i = 1; i <= spinsAB_.size() - 1; ++i ) {
		tracer << std::setw(4) << spinsAB_[i].second.rsd() << " " << std::setw(4) << spinsAB_[i].second.atomno() << " or ";
	}
	tracer << std::setw(4) << spinsAB_.back().second.rsd() << " " << std::setw(4) << spinsAB_.back().second.atomno() << " ) ";
	tracer << std::setw(6) << std::fixed << std::setprecision(3) << rdc_exp_ << " +/- " << rdc_err_ << std::endl;
	tracer.copyfmt(oldState);
}

bool operator==(RDCSingle const & lhs, RDCSingle const & rhs) {
	if ( lhs.spinsAB_.size() != rhs.spinsAB_.size() ) {
		return false;
	}
	for ( Size i = 1, i_end = lhs.spinsAB_.size(); i <= i_end; ++i ) {
		if ( (lhs.spinsAB_[i].first != rhs.spinsAB_[i].first) || (lhs.spinsAB_[i].second != rhs.spinsAB_[i].second) ) {
			return false;
		}
	}
	return true;
}

bool operator!=(RDCSingle const & lhs, RDCSingle const & rhs) {
	return !( lhs == rhs );
}

/// @brief serialize a RDCSingle object to a json_spirit object
utility::json_spirit::Value
RDCSingle::serialize() const {
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	std::vector<Value> rdc_spins_vector;
	for ( Size i = 1; i <= spinsAB_.size(); ++i ) {
		std::vector<Value> spin_pair(utility::tools::make_vector(Value(spinsAB_[i].first.serialize()), Value(spinsAB_[i].second.serialize())));
		rdc_spins_vector.push_back(Value(spin_pair));
	}
	Pair rdc_spins_record("rdc_spins", rdc_spins_vector);
	Pair rdc_exp("rdc_exp", Value(rdc_exp_));
	Pair rdc_err("rdc_err", Value(rdc_err_));
	Pair rdc_calc("rdc_calc", Value(rdc_calc_));
	Pair weight("rdc_weight", Value(weight_));
	std::vector<Value> atom_deriv_vector;
	for ( Size i = 1; i <= atom_derivatives_.size(); ++i ) {
		std::vector<Value> deriv_pair(utility::tools::make_vector(numeric::serialize(atom_derivatives_[i].first),
			numeric::serialize(atom_derivatives_[i].second)));
		atom_deriv_vector.push_back(Value(deriv_pair));
	}
	Pair atom_deriv_record("derivatives", atom_deriv_vector);

	return Value(utility::tools::make_vector(rdc_spins_record,rdc_exp,rdc_err,rdc_calc,weight,atom_deriv_record));
}

/// @brief deserialize a json_spirit object to a RDCSingle object
void
RDCSingle::deserialize(utility::json_spirit::mObject data) {
	utility::json_spirit::mArray rdc_spins_data(data["pcs_spins"].get_array());
	spinsAB_.resize(rdc_spins_data.size());
	Size i(1);
	for ( utility::json_spirit::mArray::iterator it = rdc_spins_data.begin(); it != rdc_spins_data.end(); ++it ) {
		utility::json_spirit::mArray spin_pair(it->get_array());
		core::id::AtomID spinA, spinB;
		spinA.deserialize(spin_pair[0].get_obj());
		spinB.deserialize(spin_pair[1].get_obj());
		spinsAB_[i++] = std::make_pair(spinA, spinB);
	}

	rdc_exp_  = data["rdc_exp"].get_real();
	rdc_err_  = data["rdc_err"].get_real();
	rdc_calc_ = data["rdc_calc"].get_real();
	weight_   = data["rdc_weight"].get_real();

	utility::json_spirit::mArray atom_deriv_data(data["derivatives"].get_array());
	atom_derivatives_.resize(atom_deriv_data.size());
	Size j(1);
	for ( utility::json_spirit::mArray::iterator it = atom_deriv_data.begin(); it != atom_deriv_data.end(); ++it ) {
		utility::json_spirit::mArray deriv_pair(it->get_array());
		Vector A, B;
		A = numeric::deserialize<Real>(deriv_pair[0].get_array());
		B = numeric::deserialize<Real>(deriv_pair[1].get_array());
		atom_derivatives_[j++] = std::make_pair(A,B);
	}
}

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core
