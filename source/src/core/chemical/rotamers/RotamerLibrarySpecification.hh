// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/RotamerLibrarySpecification.hh
/// @brief  The RotamerLibrarySpecification class tells how to build a Rotamer library for a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_RotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_RotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/RotamerLibrarySpecification.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>

// Basic headers

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class RotamerLibrarySpecification : public utility::pointer::ReferenceCount {
public:
	RotamerLibrarySpecification();

	virtual ~RotamerLibrarySpecification();

	/// @brief Which type of SingleResidueRotamerLibrary does this specification sub-type correspond to?
	virtual
	std::string
	keyname() const = 0;

	/// @brief How, if at all, should the corresponding SingleResidueRotamerLibrary be cached?
	///
	/// The default is to return an empty string, which turns off caching.
	///
	/// The SingleResidueRotamerLibraries are cached in the SingleResidueRotamerLibraryFactory
	/// based on keyname() and cache_tag() (as keys in a map< string, map< string, SRRL > > ).
	/// Two RotamerLibrarySpecifications with identical return values for keyname() and cache_tag()
	/// should correspond to (functionally) identical SingleResidueRotamerLibraries.
	///
	/// @details This has to be in the RotamerLibrarySpecification,
	/// as when reading we need to know the cache string before creating the library.
	///
	/// A note on writing RotamerLibrarySpecifications and SingleResidueRotamerLibrarys:
	/// The functions of a SingleResidueRotamerLibrary will normally have
	/// access to the actual RotamerLibrarySpecification from the passed Residue/ResidueType.
	/// Therefore, you don't need to store all the information from a RLS in the SRRL.
	/// Not doing so allows you to have more general cache_tag(), as the cache_tag() function
	/// only needs to disambiguate RotamerLibrarySpecifications which result in
	/// different SingleResidueRotamerLibrarys. (That is, cache_tag() only needs to encapsulate
	/// data used by SingleResidueRotamerLibraryCreator to *create* the SingleResidueRotamerLibrary.)
	///
	/// The ResidueType is passed to cache_tag() so that if the SingleResidueRotamerLibraryCreator
	/// needs details from the ResidueType in order to correctly create the SingleResidueRotamerLibrary,
	/// that information can be extracted.
	/// In general, though, you want to avoid keying off of information in ResidueType as much as possible.
	virtual
	std::string
	cache_tag(core::chemical::ResidueType const &) const { return ""; }

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
