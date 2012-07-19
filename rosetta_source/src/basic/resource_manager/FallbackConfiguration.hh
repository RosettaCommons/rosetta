// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfiguration.hh
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_HH

// Unit Headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace resource_manager {

class FallbackConfiguration : public utility::pointer::ReferenceCount {

public:
	FallbackConfiguration() {}
	
	virtual ~FallbackConfiguration();
	
	virtual
	ResourceTag const &
	get_resource_tag_from_description( ResourceDescription const & desc ) const = 0;

	virtual
	LocatorTag const &
	get_locator_tag_from_description( ResourceDescription const & desc ) const = 0;

	virtual
	LocatorID const &
	get_locator_id_from_description( ResourceDescription const & desc ) const = 0;

	virtual
	LoaderType const &
	get_loader_type_from_description( ResourceDescription const & desc ) const = 0;

	virtual
	ResourceOptionsTag const &
	get_resource_options_tag_from_description( ResourceDescription const & desc ) const = 0;

};

} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_fallback_configuration_HH
