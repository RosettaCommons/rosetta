// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryRegistrator.hh
/// @brief  Declaration of the template class for registrating Chemistry objects with
///         the ChemistryFactory
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_ChemistryRegistrator_hh
#define INCLUDED_protocols_chemistries_ChemistryRegistrator_hh

// Package headers
#include <protocols/chemistries/ChemistryFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace protocols {
namespace chemistries {

/// @brief This templated class will register an instance of an
/// ChemistryCreator (class T) with the ChemistryFactory.  It will ensure
/// that no ChemistryCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class ChemistryRegistrator : public utility::factory::WidgetRegistrator< ChemistryFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ChemistryFactory, T > parent;
public:
	ChemistryRegistrator() : parent() {}
};

}
}

#endif
