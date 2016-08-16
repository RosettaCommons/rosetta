// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibraryRegistrator.hh
/// @brief  Declaration of the template class for registrating SingleResidueRotamerLibraryCreators with
///         the SingleResidueRotamerLibraryFactory
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_pack_rotamers_SingleResidueRotamerLibraryRegistrator_hh
#define INCLUDED_core_pack_rotamers_SingleResidueRotamerLibraryRegistrator_hh

// Package headers
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace pack {
namespace rotamers {

/// @brief This templated class will register an instance of an
/// SingleResidueRotamerLibraryCreator (class T) with the SingleResidueRotamerLibraryFactory.  It will ensure
/// that no SingleResidueRotamerLibrary creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class SingleResidueRotamerLibraryRegistrator : public utility::factory::WidgetRegistrator< SingleResidueRotamerLibraryFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< SingleResidueRotamerLibraryFactory, T > parent;
public:
	SingleResidueRotamerLibraryRegistrator() : parent() {}
};

} // rotamers
} // pack
} // core

#endif
