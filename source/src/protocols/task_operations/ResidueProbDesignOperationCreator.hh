// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Sample a set of mutations each time packer is generated.
/// @author Samuel Schmitz (Samuel.Schmitz@vanderbilt.edu)

#ifndef INCLUDED_protocols_task_operations_ResidueProbDesignOperationCreator_hh
#define INCLUDED_protocols_task_operations_ResidueProbDesignOperationCreator_hh


#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>

namespace protocols {
namespace task_operations {

class ResidueProbDesignOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	core::pack::task::operation::TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace task_operations
} //namespace protocols
#endif

