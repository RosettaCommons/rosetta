// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/BasicRotamerLibrarySpecification.hh
/// @brief  The BasicRotamerLibrarySpecification class says to do simple rotamer libraries
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_BasicRotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_BasicRotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/BasicRotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>

// C++ headers
#include <iosfwd>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace rotamers {

class BasicRotamerLibrarySpecification : public RotamerLibrarySpecification {
public:
	BasicRotamerLibrarySpecification();
	virtual ~BasicRotamerLibrarySpecification();

	/// @brief Which type of SingleResidueRotamerLibrary does this specification sub-type correspond to?
	virtual
	std::string
	keyname() const;

	virtual
	std::string
	cache_tag(ResidueType const &) const;

	/// @brief Static function for access to type_name, to have a single string which is used for both
	/// this class and for the SingleResidueRotamerLibraryCreator.
	static
	std::string
	library_name();

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_rotamers_BasicRotamerLibrarySpecification )
#endif // SERIALIZATION


#endif
