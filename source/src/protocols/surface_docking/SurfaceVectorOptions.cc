// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/surface_docking/SurfaceVectorOptions.cc
/// @brief
/// @author Michael Pacella (mpacella88@gmail.com)

//unit headers
#include <protocols/surface_docking/SurfaceVectorOptions.hh>
#include <protocols/surface_docking/SurfaceVectorOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>

//C++ headers

namespace protocols {
namespace surface_docking {

std::string
SurfaceVectorOptionsCreator::options_type() const { return "SurfaceVectorOptions"; }

basic::resource_manager::ResourceOptionsOP
SurfaceVectorOptionsCreator::create_options() const { return new SurfaceVectorOptions; }


SurfaceVectorOptions::SurfaceVectorOptions() : basic::resource_manager::ResourceOptions() {}
SurfaceVectorOptions::~SurfaceVectorOptions() {}

void
SurfaceVectorOptions::parse_my_tag(
	utility::tag::TagCOP
)
{
	return;
}

std::string
SurfaceVectorOptions::type() const
{
	return "SurfaceVectorOptions";
}
} // namespace surface_docking
} // namespace protocols
