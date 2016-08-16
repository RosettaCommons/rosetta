// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/AntibodyInfoRMLoader.hh
/// @brief A hacky way to create AntibodyInfo objects through the RM
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_AntibodyInfoRMLoader_HH
#define INCLUDED_protocols_surface_docking_AntibodyInfoRMLoader_HH

//unit headers
#include <protocols/antibody/AntibodyInfoRMLoader.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace surface_docking {

/// @brief AntibodyInfoRMLoader constructs an AntibodyInfo instance from data provided by the resource manager
/// @details The AntibodyInfoRMLoader ignores the given istream
class AntibodyInfoRMLoader : public basic::resource_manager::ResourceLoader
{
public:
	AntibodyInfoRMLoader();
	virtual ~AntibodyInfoRMLoader();

	/// @brief Returns a SurfaceParametersOP which is constructed from the given input
	/// stream (istream).
	virtual
	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
	) const;
	
	/// @brief Returns the default options for AntibodyInfo
	virtual
	basic::resource_manager::ResourceOptionsOP
	default_options() const;

};

} // namespace antibody
} // namespace protocols

#endif //INCLUDED_protocols_antibody_AntibodyInfoRMLoader_HH
