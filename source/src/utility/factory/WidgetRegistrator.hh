// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/factory/WidgetRegistrator.hh
/// @brief  WidgetRegistrator class which registers WidgetCreator classes at load time
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_factory_WidgetRegistrator_hh
#define INCLUDED_utility_factory_WidgetRegistrator_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace factory {


/// @brief This class will register a Creator with a Factory at load time.
///
/// No forward header for this class since it should never appear listed
/// in a parameter for a function or contained in an owning pointer in some
/// other class.
template< class FACTORY, class CREATOR >
class WidgetRegistrator
{
	typedef utility::pointer::owning_ptr< CREATOR > CREATOROP;

public:

	WidgetRegistrator() {
		FACTORY::get_instance()->factory_register( CREATOROP( new CREATOR ) );
	}

};

} // namespace factory
} // namespace utility


#endif // INCLUDED_utility_factory_WidgetRegistrator_HH
