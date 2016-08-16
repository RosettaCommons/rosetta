// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/DunbrackRotamerLibrarySpecification.hh
/// @brief  The DunbrackRotamerLibrarySpecification class tells how to build DunbrackRotamers
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_DunbrackRotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_DunbrackRotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/DunbrackRotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AA.hh>

// C++ headers
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class DunbrackRotamerLibrarySpecification : public RotamerLibrarySpecification {
public:
	DunbrackRotamerLibrarySpecification();
	DunbrackRotamerLibrarySpecification( AA aa );
	DunbrackRotamerLibrarySpecification( std::istream & input );
	virtual ~DunbrackRotamerLibrarySpecification();

	/// @brief The AA for which we're building the Rotamer library.
	AA get_aa() const {
		assert( aa_ != aa_unk ); // uninitialized
		return aa_;
	}

	void aa( AA aa );

	/// @brief Which type of SingleResidueRotamerLibrary does this specification sub-type correspond to?
	virtual
	std::string
	keyname() const;

	/// @brief Return empty string, as Dunbrack rotamer caching happens external to the Factory caching.
	virtual
	std::string
	cache_tag(ResidueType const &) const { return ""; }

	/// @brief Static function for access to type_name, to have a single string which is used for both
	/// this class and for the SingleResidueRotamerLibraryCreator.
	static
	std::string
	library_name();

private:
	AA aa_;

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
