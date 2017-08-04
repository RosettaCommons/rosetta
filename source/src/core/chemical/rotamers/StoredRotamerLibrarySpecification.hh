// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/StoredRotamerLibrarySpecification.hh
/// @brief  The StoredRotamerLibrarySpecification class says to build a rotamer library from a set of stored coordinates
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_StoredRotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_StoredRotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/types.hh>

// Basic headers

// Utility Headers
#include <utility/vector1.hh>

// C++ headers
#include <istream>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace rotamers {

/// @brief A class which stores atom coordinates for a rotamer library.
/// Internally, this is stored as a list of name:coordinate maps.
/// @details
/// This is intended as an analog to PDBRotamerLibrarySpecifications
/// for those instances where the coordinates are generated internally or
/// programmatically.

class StoredRotamerLibrarySpecification : public RotamerLibrarySpecification {
public:
	StoredRotamerLibrarySpecification();
	StoredRotamerLibrarySpecification(std::istream & input);
	virtual ~StoredRotamerLibrarySpecification();

	utility::vector1< std::map< std::string, core::Vector > > const &
	coordinates() const { return coordinates_; }

	/// @brief Add a particular rotamer to the list
	void
	add_rotamer( std::map< std::string, core::Vector > const & rotamer );

	virtual
	std::string
	keyname() const;

	static std::string library_name();

private:

	utility::vector1< std::map< std::string, core::Vector > > coordinates_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_rotamers_StoredRotamerLibrarySpecification )
#endif // SERIALIZATION


#endif
