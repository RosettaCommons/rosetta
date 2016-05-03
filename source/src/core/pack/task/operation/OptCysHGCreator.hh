// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_OptCysHGCreator_hh
#define INCLUDED_core_pack_task_operation_OptCysHGCreator_hh

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <core/pack/task/operation/TaskOperation.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class OptCysHGCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
