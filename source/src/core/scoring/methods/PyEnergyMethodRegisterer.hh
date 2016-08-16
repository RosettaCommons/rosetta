// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/PyEnergyMethodRegisterer.hh
/// @brief  EnergyMethod factory registration and creation in PyRosetta
/// @author Sergey Lyskov


#ifndef INCLUDED_core_scoring_methods_PyEnergyMethodRegisterer_hh
#define INCLUDED_core_scoring_methods_PyEnergyMethodRegisterer_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief This class will register an instance of an EnergyMethodCreator (class T) with the ScoringManager.
/// It will ensure that no energy method creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in one place

class PyEnergyMethodRegistrator
{
public:
	PyEnergyMethodRegistrator(EnergyMethodCreatorOP CreatorOP)
	{
		ScoringManager::get_instance()->factory_register( CreatorOP );
	}

};

}
}
}

#endif
