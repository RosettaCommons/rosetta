// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingle.cc
/// @brief   Implementation of class PCSSingle
/// @details last Modified: 06/21/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingle.hh>

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
#include <vector>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "core.scoring.nmr.pcs.PCSSingle" );

/// @brief default constructor
PCSSingle::PCSSingle() :
	protein_spins_(),
	pcs_exp_(0),
	pcs_err_(0),
	pcs_calc_(0),
	weight_(1.0),
	atom_derivatives_()
{}

/// @brief constructor with arguments
PCSSingle::PCSSingle(
	utility::vector1< core::io::nmr::AtomSelection > const & spins,
	pose::Pose const & pose,
	Real const pcs_exp,
	Real const pcs_err
) :
	protein_spins_(),
	pcs_exp_(pcs_exp),
	pcs_err_(pcs_err),
	pcs_calc_(0),
	weight_(1.0),
	atom_derivatives_()
{
	runtime_assert_msg(spins.size() > 0, "ERROR during creation of PCSSingle object. Vector of spin atom selections must not be empty.");

	using utility::to_string;

	core::pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	for ( Size i = 1, i_end = spins.size(); i <= i_end; ++i ) {
		Size rsd = spins[i].get_rsd();
		std::string atom = spins[i].get_atom();
		char const & chain = spins[i].get_chain();

		if ( pdbinfo ) {
			TR.Trace << "Converting AtomSelection \"" + to_string(rsd) + " " + atom + " " + to_string(chain) + "\" in PCS datafile from pdb to pose numbering." << std::endl;
			rsd = pdbinfo->pdb2pose(chain, rsd);
			if ( rsd == 0 ) { // Residue not found
				utility_exit_with_message( "ERROR: Cannot convert AtomSelection " + to_string(spins[i].get_rsd()) + " " + atom + " " + to_string(chain) + " from PCS datafile into pose numbering. Residue number was not found." );
			}
		} else { // Assume pose numbering instead
			TR.Warning << "Cannot convert residue numbers from pdb to pose numbering. No PDBInfo object. Make sure to provide residues in pose numbering in PCS datafile instead." << std::endl;
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
PCSSingle::PCSSingle(PCSSingle const & other) :
	protein_spins_(other.protein_spins_),
	pcs_exp_(other.pcs_exp_),
	pcs_err_(other.pcs_err_),
	pcs_calc_(other.pcs_calc_),
	weight_(other.weight_),
	atom_derivatives_(other.atom_derivatives_)
{}

/// @brief assignment operator
PCSSingle&
PCSSingle::operator=(PCSSingle const & rhs) {
	if ( this != & rhs ) {
		protein_spins_ = rhs.protein_spins_;
		pcs_exp_ = rhs.pcs_exp_;
		pcs_err_ = rhs.pcs_err_;
		pcs_calc_ = rhs.pcs_calc_;
		weight_ = rhs.weight_;
		atom_derivatives_ = rhs.atom_derivatives_;
	}
	return *this;
}

/// @brief destructor
PCSSingle::~PCSSingle() {}

void
PCSSingle::show(std::ostream & tracer) const {
	// Store old iostream manipulator flags
	std::ios oldState(nullptr);
	oldState.copyfmt(tracer);
	tracer << "( ";
	for ( Size i = 1; i <= protein_spins_.size()-1; ++i ) {
		tracer << std::setw(4) << protein_spins_[i].rsd() << " " << std::setw(4) << protein_spins_[i].atomno() << " or ";
	}
	tracer << std::setw(4) << protein_spins_.back().rsd() << " " << std::setw(4) << protein_spins_.back().atomno() << " ) ";
	tracer << std::setw(6) << std::fixed << std::setprecision(3) << pcs_exp_ << " +/- " << pcs_err_ << std::endl;
	tracer.copyfmt(oldState);
}

bool
operator==(PCSSingle const & lhs, PCSSingle const & rhs) {
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

bool
operator!=(PCSSingle const & lhs, PCSSingle const & rhs) {
	return !( lhs == rhs );
}

/// @brief serialize a PCSSingle object to a json_spirit object
utility::json_spirit::Value
PCSSingle::serialize() const {
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	std::vector<Value> pcs_spins_vector;
	for ( Size i = 1; i <= protein_spins_.size(); ++i ) {
		pcs_spins_vector.push_back(protein_spins_[i].serialize());
	}
	Pair pcs_spins_record("pcs_spins", pcs_spins_vector);
	Pair pcs_exp("pcs_exp", Value(pcs_exp_));
	Pair pcs_err("pcs_err", Value(pcs_err_));
	Pair pcs_calc("pcs_calc", Value(pcs_calc_));
	Pair weight("pcs_weight", Value(weight_));
	std::vector<Value> atom_deriv_vector;
	for ( Size i = 1; i <= atom_derivatives_.size(); ++i ) {
		atom_deriv_vector.push_back(numeric::serialize(atom_derivatives_[i]));
	}
	Pair atom_deriv_record("derivatives", atom_deriv_vector);

	return Value(utility::tools::make_vector(pcs_spins_record,pcs_exp, pcs_err, pcs_calc, weight,atom_deriv_record));
}

/// @brief deserialize a json_spirit object to a PCSSingle object
void
PCSSingle::deserialize(utility::json_spirit::mObject data) {
	utility::json_spirit::mArray pcs_spins_data(data["pcs_spins"].get_array());
	protein_spins_.resize(pcs_spins_data.size());
	Size i(1);
	for ( utility::json_spirit::mArray::iterator it = pcs_spins_data.begin(); it != pcs_spins_data.end(); ++it ) {
		core::id::AtomID id;
		id.deserialize(it->get_obj());
		protein_spins_[i++] = id;
	}

	pcs_exp_  = data["pcs_exp"].get_real();
	pcs_err_  = data["pcs_err"].get_real();
	pcs_calc_ = data["pcs_calc"].get_real();
	weight_   = data["pcs_weight"].get_real();

	utility::json_spirit::mArray atom_deriv_data(data["derivatives"].get_array());
	atom_derivatives_.resize(atom_deriv_data.size());
	Size j(1);
	for ( utility::json_spirit::mArray::iterator it = atom_deriv_data.begin(); it != atom_deriv_data.end(); ++it ) {
		atom_derivatives_[j++] = numeric::deserialize<Real>(it->get_array());
	}
}

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core
