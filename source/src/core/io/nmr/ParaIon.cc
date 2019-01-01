// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParaIon.cc
/// @brief   Function definitions of ParaIon class
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)


// Unit headers
#include <core/io/nmr/ParaIon.hh>

// Package headers
#include <core/types.hh>

// C++ headers
#include <iostream>
#include <iomanip>
#include <string>

namespace core {
namespace io {
namespace nmr {

/// @brief default constructor
ParaIon::ParaIon() :
	ion_label_("Nitroxide"),
	charge_(0.0),
	S_(0.5),
	L_(0.0),
	J_(0.5),
	gJ_(2.0),
	tau_e_(100000.0)
{}

/// @brief constructor with arguments
ParaIon::ParaIon(
	std::string const & label,
	Real const charge,
	Real const s,
	Real const l,
	Real const j,
	Real const gj,
	Real const te
) :
	ion_label_(label),
	charge_(charge),
	S_(s),
	L_(l),
	J_(j),
	gJ_(gj),
	tau_e_(te)
{}

/// @brief copy constructor
/// @details Make a deep copy of this mover object
ParaIon::ParaIon( ParaIon const & src ) :
	ion_label_(src.ion_label_),
	charge_(src.charge_),
	S_(src.S_),
	L_(src.L_),
	J_(src.J_),
	gJ_(src.gJ_),
	tau_e_(src.tau_e_)
{}

/// @brief assignment operator
/// @details Make a deep copy of this mover object
ParaIon &
ParaIon::operator=( ParaIon const & src ) {
	if ( this != & src ) {
		ion_label_ = src.ion_label_;
		charge_ = src.charge_;
		S_ = src.S_;
		L_ = src.L_;
		J_ = src.J_;
		gJ_ = src.gJ_;
		tau_e_ = src.tau_e_;
	}
	return *this;
}

/// @brief destructor
ParaIon::~ParaIon() {}

/// @brief calculate Lande g-factor
void
ParaIon::calc_gJ() {
	gJ_ = (Real)(J_*(J_+1) - L_*(L_+1) + S_*(S_+1))/(2.0*J_*(J_+1));
}

/// @brief output ion data
void
ParaIon::show(std::ostream & out) const {
	std::ios::fmtflags f(out.flags());
	out << "ion type                    : " << ion_label_ << ", ";
	out << std::setprecision(3) << std::fixed;
	out << "spin quantum number S       : " << S_ << ", ";
	out << "orbital quantum number L    : " << L_ << ", ";
	out << "total quantum number J      : " << J_ << ", ";
	out << "Lande g-factor gJ           : " << gJ_ << ", ";
	out << "electron relaxation time te : " << std::scientific << tau_e_ << std::endl;
	out.flags(f);
}

} // namespace nmr
} // namespace io
} // namespace core
