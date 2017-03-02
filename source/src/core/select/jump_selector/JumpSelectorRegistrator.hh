// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpSelectorRegistrator.hh
/// @brief  Declaration of the template class for registrating JumpSelectorCreators with
///         the JumpSelectorFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_select_jump_selector_JumpSelectorRegistrator_hh
#define INCLUDED_core_select_jump_selector_JumpSelectorRegistrator_hh

// Package headers
#include <core/select/jump_selector/JumpSelectorFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace select {
namespace jump_selector {

/// @brief This templated class will register an instance of an
/// JumpSelectorCreator (class T) with the JumpSelectorFactory.  It will ensure
/// that no JumpSelector creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class JumpSelectorRegistrator : public utility::factory::WidgetRegistrator< JumpSelectorFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< JumpSelectorFactory, T > parent;
public:
	JumpSelectorRegistrator() : parent() {}
};

}
}
}

#endif
