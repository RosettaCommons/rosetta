// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/AntibodyInfoRMOptionsCreator.hh
/// @brief
/// @author 

#ifndef INCLUDED_protocols_antibody_AntibodyInfoRMOptionsCreator_hh
#define INCLUDED_protocols_antibody_AntibodyInfoRMOptionsCreator_hh

//unit headers

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers

namespace protocols {
namespace antibody {

class AntibodyInfoRMOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:

	virtual std::string options_type() const;
	virtual basic::resource_manager::ResourceOptionsOP create_options() const;

};

} // namespace antibody
} // namespace protocols

#endif //INCLUDED_protocols_antibody_AntibodyInfoRMOptionsCreator_hh
