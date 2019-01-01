// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/AtomSelection.cc
/// @brief   Implementation of class AtomSelection
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/io/nmr/AtomSelection.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace io {
namespace nmr {

AtomSelection::AtomSelection() :
	rsd_(1),
	atom_(""),
	chain_('^')
{ }

AtomSelection::AtomSelection(
	Size const residue,
	std::string const & atom,
	char const & chain
) :
	rsd_(residue),
	atom_(atom),
	chain_(toupper(chain))
{
	runtime_assert_msg( residue > 0 , "ERROR: Residue number must be positive");
}

AtomSelection::AtomSelection(AtomSelection const & other) :
	rsd_(other.rsd_),
	atom_(other.atom_),
	chain_(other.chain_)
{}

AtomSelection & AtomSelection::operator=(AtomSelection const & rhs) {
	if ( this != &rhs ) {
		rsd_ = rhs.rsd_;
		atom_ = rhs.atom_;
		chain_ = rhs.chain_;
	}
	return *this;
}

AtomSelection::~AtomSelection() {}

bool operator<(AtomSelection const & lhs, AtomSelection const & rhs) {
	if ( lhs.rsd_   < rhs.rsd_   ) return true;
	if ( lhs.rsd_   > rhs.rsd_   ) return false;

	if ( lhs.atom_  < rhs.atom_  ) return true;
	if ( lhs.atom_  > rhs.atom_  ) return false;

	if ( lhs.chain_ < rhs.chain_ ) return true;

	return false;
}

bool operator>(AtomSelection const & lhs, AtomSelection const & rhs) {
	return rhs < lhs;
}

bool operator<=(AtomSelection const & lhs, AtomSelection const & rhs) {
	return !(lhs > rhs);
}

bool operator>=(AtomSelection const & lhs, AtomSelection const & rhs) {
	return !(lhs < rhs);
}

bool operator==(AtomSelection const & lhs, AtomSelection const & rhs) {
	return ( (lhs.rsd_ == rhs.rsd_) && (lhs.atom_ == rhs.atom_) && (lhs.chain_ == rhs.chain_) );
}

bool operator!=(AtomSelection const & lhs, AtomSelection const & rhs) {
	return !( lhs == rhs );
}

} // namespace nmr
} // namespace io
} // namespace core
