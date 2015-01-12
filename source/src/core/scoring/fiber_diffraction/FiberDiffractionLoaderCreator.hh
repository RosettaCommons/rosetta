// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionLoaderCreator.hh
/// @brief  Fiber diffraction data creator
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionLoaderCreator_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionLoaderCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>
#include <core/scoring/fiber_diffraction/FiberDiffractionLoader.fwd.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>

namespace core {
namespace scoring {
namespace fiber_diffraction {

class FiberDiffractionLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:
	FiberDiffractionLoaderCreator();
	virtual ~FiberDiffractionLoaderCreator();

	virtual basic::resource_manager::ResourceLoaderOP create_resource_loader() const;
	virtual std::string loader_type() const;
};

} //namespace
} //namespace
} //namespace

#endif
