// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/electron_density/ElectronDensityOptionsCreator.hh
/// @brief  Header for ElectronDensityOptions Creator for the load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_core_scoring_electron_density_ElectronDensityOptionsCreator_hh
#define INCLUDED_core_scoring_electron_density_ElectronDensityOptionsCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <core/scoring/electron_density/ElectronDensityOptions.fwd.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>

namespace core {
namespace scoring {
namespace electron_density {

/// @brief creator for the ElectronDensityOptions class
class ElectronDensityOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:
	ElectronDensityOptionsCreator();
	virtual ~ElectronDensityOptionsCreator();

	virtual basic::resource_manager::ResourceOptionsOP create_options() const;
	virtual std::string options_type() const;
};

} //namespace
} //namespace
} //namespace

#endif
