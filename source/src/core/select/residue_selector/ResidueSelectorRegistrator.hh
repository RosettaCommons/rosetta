// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSelectorRegistrator.hh
/// @brief  Declaration of the template class for registrating ResidueSelectorCreators with
///         the ResidueSelectorFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_select_residue_selector_ResidueSelectorRegistrator_hh
#define INCLUDED_core_select_residue_selector_ResidueSelectorRegistrator_hh

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace select {
namespace residue_selector {

/// @brief This templated class will register an instance of an
/// ResidueSelectorCreator (class T) with the ResidueSelectorFactory.  It will ensure
/// that no ResidueSelector creator is registered twice, and centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class ResidueSelectorRegistrator : public utility::factory::WidgetRegistrator< ResidueSelectorFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResidueSelectorFactory, T > parent;
public:
	ResidueSelectorRegistrator() : parent() {}
};

}
}
}

#endif
