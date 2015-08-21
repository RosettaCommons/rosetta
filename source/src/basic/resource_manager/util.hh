// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/util.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_util_hh
#define INCLUDED_basic_resource_manager_util_hh

//unit headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/types.hh>

//project headers
#include <utility/excn/Exceptions.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

// numeric headers

//C++ headers
#include <sstream>

namespace basic {
namespace resource_manager {

using std::stringstream;

template< class ResourceType >
utility::pointer::shared_ptr< ResourceType >
get_resource(
	ResourceDescription const & resource_description
) {
	utility::pointer::shared_ptr< ResourceType > resource(
		utility::pointer::dynamic_pointer_cast< ResourceType > ( ResourceManager::get_instance()->get_resource(resource_description) ));

	if ( !resource ) {
		stringstream err_msg;
		err_msg
			<< "The '" << resource_description << "' "
			<< "cannot be cast to the given type.";
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}
	return resource;
}

} // namespace resource_manager
} // namespace basic

#endif
