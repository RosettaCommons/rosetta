// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/residue_selector/RotamerLibrarySpecificationFactory.hh
/// @brief  class for instantiating RotamerLibrarySpecifications from strings
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_RotamerLibrarySpecificationFactory_HH
#define INCLUDED_core_chemical_rotamers_RotamerLibrarySpecificationFactory_HH

// Package headers
#include <core/chemical/rotamers/RotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecificationCreator.hh>

// Program headers

// Basic headers

// Utility headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <map>
#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class RotamerLibrarySpecificationFactory : public utility::SingletonBase< RotamerLibrarySpecificationFactory > {
public:

	void factory_register( RotamerLibrarySpecificationCreatorOP creator );

	bool has_type( std::string const & ) const;

	core::chemical::rotamers::RotamerLibrarySpecificationOP
	get( std::string const & ) const;

	core::chemical::rotamers::RotamerLibrarySpecificationOP
	get( std::string const &, std::istream & ) const;

public:
	static RotamerLibrarySpecificationFactory * create_singleton_instance();

private:

	RotamerLibrarySpecificationCreatorCOP
	get_creator( std::string const & tag ) const;

private:
	RotamerLibrarySpecificationFactory();
	RotamerLibrarySpecificationFactory( RotamerLibrarySpecificationFactory const & ); // unimplemented
	RotamerLibrarySpecificationFactory const & operator = ( RotamerLibrarySpecificationFactory const & ); // unimplemented

private:
	typedef std::map< std::string, RotamerLibrarySpecificationCreatorOP > CreatorMap;
	CreatorMap creator_map_;

};

} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
