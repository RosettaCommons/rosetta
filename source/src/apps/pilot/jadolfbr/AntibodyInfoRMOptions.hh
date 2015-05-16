// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/AntibodyInfoRMOptions.hh
/// @brief
/// @author 

#ifndef INCLUDED_protocols_antibody_AntibodyInfoRMOptions_hh
#define INCLUDED_protocols_antibody_AntibodyInfoRMOptions_hh

//unit headers
#include <protocols/antibody/AntibodyInfoRMOptions.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceOptions.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers

namespace protocols {
namespace antibody {

class AntibodyInfoRMOptions : public basic::resource_manager::ResourceOptions
{
public:
	AntibodyInfoRMOptions();
	virtual ~AntibodyInfoRMOptions();

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	);

	virtual
	std::string
	type() const;

};

} // namespace antibody
} // namespace protocols



#endif //INCLUDED_protocols_antibody_AntibodyInfoRMOptions_hh
