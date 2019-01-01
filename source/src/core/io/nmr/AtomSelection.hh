// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/AtomSelection.hh
/// @brief   Class that holds selection identifiers for a single atom.
///          An atom list is read in from a general NMR input/assignment file.
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_AtomSelection_HH
#define INCLUDED_core_io_nmr_AtomSelection_HH

// Unit headers
#include <core/io/nmr/AtomSelection.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ headers
#include <string>
#include <iostream>


namespace core {
namespace io {
namespace nmr {

class AtomSelection {

public: // Methods

	/// @brief default constructor
	AtomSelection();

	/// @brief create from residue number, atom name and chain ID
	AtomSelection(
		Size const residue,
		std::string const & atom,
		char const & chain
	);

	/// @brief copy constructor
	AtomSelection(AtomSelection const & other);

	/// @brief copy assignment
	AtomSelection &
	operator=(AtomSelection const & rhs);

	/// @brief destructor
	~AtomSelection();

	// Setters and Getters
	Size get_rsd() const { return rsd_; }
	std::string get_atom() const { return atom_; }
	const char & get_chain() const { return chain_; }
	void set_rsd(Size const residue) { rsd_ = residue; }
	void set_atom(std::string const & atom) { atom_ = atom; }
	void set_chain(char const & chain) { chain_ = toupper(chain); }
	void show(std::ostream & TR) const { TR << rsd_ << " " << atom_ << " " << chain_ << std::endl; }

	friend bool operator<(AtomSelection const & lhs, AtomSelection const & rhs);
	friend bool operator>(AtomSelection const & lhs, AtomSelection const & rhs);
	friend bool operator<=(AtomSelection const & lhs, AtomSelection const & rhs);
	friend bool operator>=(AtomSelection const & lhs, AtomSelection const & rhs);
	friend bool operator==(AtomSelection const & lhs, AtomSelection const & rhs);
	friend bool operator!=(AtomSelection const & lhs, AtomSelection const & rhs);

private: // Data

	Size rsd_;
	std::string atom_;
	char chain_;

};



} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_AtomSelection_HH
