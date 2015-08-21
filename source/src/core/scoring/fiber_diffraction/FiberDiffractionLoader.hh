// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionLoader.hh
/// @brief  Options for constructing fiber diffraction layer lines
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionLoader_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionLoader_hh


// Unit Headers
#include <core/scoring/fiber_diffraction/FiberDiffractionLoader.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

// Platform Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <string>
#include <istream>

namespace core {
namespace scoring {
namespace fiber_diffraction {

class FiberDiffractionLoader : public basic::resource_manager::ResourceLoader
{
public:
	FiberDiffractionLoader();

	~FiberDiffractionLoader();

	FiberDiffractionLoader(
		FiberDiffractionLoader const & src);

	basic::resource_manager::ResourceOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
	) const;

	basic::resource_manager::ResourceOptionsOP
	default_options() const;

};


} // namespace
} // namespace
} // namespace


#endif
