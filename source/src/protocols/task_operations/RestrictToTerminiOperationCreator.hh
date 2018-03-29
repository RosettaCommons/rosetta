// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictToTerminiOperation.fwd.hh
/// @brief  RestrictToTerminiOperation class creator
/// @author Arpit Tandon
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_task_operations_RestrictToTerminiOperationCreator_hh
#define INCLUDED_protocols_task_operations_RestrictToTerminiOperationCreator_hh

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>


namespace protocols {
namespace task_operations {

class RestrictToTerminiOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace
} //namespace

#endif // include guard

