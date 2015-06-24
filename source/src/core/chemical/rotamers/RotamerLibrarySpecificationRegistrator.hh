// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/RotamerLibrarySpecificationRegistrator.hh
/// @brief  Declaration of the template class for registrating RotamerLibrarySpecificationCreators with
///         the RotamerLibrarySpecificationFactory
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_rotamers_RotamerLibrarySpecificationRegistrator_hh
#define INCLUDED_core_chemical_rotamers_RotamerLibrarySpecificationRegistrator_hh

// Package headers
#include <core/chemical/rotamers/RotamerLibrarySpecificationFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace chemical {
namespace rotamers {

/// @brief This templated class will register an instance of an
/// RotamerLibrarySpecificationCreator (class T) with the RotamerLibrarySpecificationFactory.  It will ensure
/// that no RotamerLibrarySpecification creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class RotamerLibrarySpecificationRegistrator : public utility::factory::WidgetRegistrator< RotamerLibrarySpecificationFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< RotamerLibrarySpecificationFactory, T > parent;
public:
	RotamerLibrarySpecificationRegistrator() : parent() {}
};

} // rotamers
} // chemical
} // core

#endif
