// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EnergyMethodRegistrator.hh
/// @brief  Declaration of the base class for EnergyMethod factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_methods_EnergyMethodRegistrator_hh
#define INCLUDED_core_scoring_methods_EnergyMethodRegistrator_hh

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace scoring {
namespace methods {

/// @brief This templated class will register an instance of an
/// EnergyMethodCreator (class T) with the ScoringManager.  It will ensure
/// that no energy method creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class EnergyMethodRegistrator : public utility::factory::WidgetRegistrator< ScoringManager, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ScoringManager, T > parent;
public:
	EnergyMethodRegistrator() : parent() {}
	/*
	{
	register_with_scoring_manager();
	}

	private:

	static bool registered_;
	static void register_with_scoring_manager() {
	if ( registered_ ) return;
	/// NOT THREADSAFE!
	registered_ = true;
	ScoringManager::get_instance()->factory_register( new T );
	}
	*/
};

/*template < class T >
bool EnergyMethodRegistrator< T >::registered_( false );*/

}
}
}

#endif
