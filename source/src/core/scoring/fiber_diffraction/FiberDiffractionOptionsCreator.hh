// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionOptionsCreator.hh
/// @brief  Options for fiber diffraction data
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionOptionsCreator_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionOptionsCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <core/scoring/fiber_diffraction/FiberDiffractionOptions.fwd.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>

namespace core {
namespace scoring {
namespace fiber_diffraction {

/// @brief creator for the FiberDiffractionOptions class
class FiberDiffractionOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:
	FiberDiffractionOptionsCreator();
	virtual ~FiberDiffractionOptionsCreator();

	virtual basic::resource_manager::ResourceOptionsOP create_options() const;
	virtual std::string options_type() const;
};

} //namespace
} //namespace
} //namespace

#endif
