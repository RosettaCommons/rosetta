// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/factory/WidgetFactory.fwd.hh
/// @brief  WidgetFactory forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_factory_WidgetFactory_fwd_hh
#define INCLUDED_utility_factory_WidgetFactory_fwd_hh


namespace utility {
namespace factory {


/// @brief Factory base class holds a map between strings and owning pointers of
/// the creator classes.
template< class Creator >
class WidgetFactory;

} // namespace factory
} // namespace utility


#endif // INCLUDED_utility_factory_WidgetFactory_FWD_HH
