// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/AntibodyInfoRMOptions.cc
/// @brief
/// @author 

//unit headers
#include <protocols/antibody/AntibodyInfoRMOptions.hh>
#include <protocols/antibody/AntibodyInfoRMOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>

//C++ headers

namespace protocols {
namespace antibody {

std::string
AntibodyInfoRMOptionsCreator::options_type() const { return "AntibodyInfoRMOptions"; }

basic::resource_manager::ResourceOptionsOP
AntibodyInfoRMOptionsCreator::create_options() const { return basic::resource_manager::ResourceOptionsOP( new AntibodyInfoRMOptions ); }


AntibodyInfoRMOptions::AntibodyInfoRMOptions() : basic::resource_manager::ResourceOptions() {}
AntibodyInfoRMOptions::~AntibodyInfoRMOptions() {}

void
AntibodyInfoRMOptions::parse_my_tag(
	utility::tag::TagCOP
)
{
	return;
}

std::string
AntibodyInfoRMOptions::type() const
{
	return "AntibodyInfoRMOptions";
}
} // namespace antibody
} // namespace protocols
