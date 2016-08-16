// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/methods/ResFilterCreator.hh
/// @brief  Declaration of the base class for ResFilter factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_pack_task_operation_ResFilterCreator_hh
#define INCLUDED_core_pack_task_operation_ResFilterCreator_hh

// Unit headers
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>

// Package headers
#include <core/pack/task/operation/ResFilter.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief The ResFilterCreator class's responsibilities are to create
/// on demand a new ResFilter class.
/// The ResFilterCreator must register itself with the ResFilterFactory
/// at load time (before main() begins) so that the ResFilterFactory is ready
/// to start creating ResFilters by the time any protocol
/// requests one.
class ResFilterCreator : public utility::pointer::ReferenceCount
{
public:
	/// @brief Instantiate a new ResFilter
	virtual ResFilterOP create_res_filter() const = 0;
	virtual std::string keyname() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;

};

}
}
}
}

#endif
