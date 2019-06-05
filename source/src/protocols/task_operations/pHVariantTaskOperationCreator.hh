// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/pHVariantTaskOperation.cc
/// @brief  During repacking, allow interchangeability of protonated and deprotonated variants
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_task_operations_pHVariantTaskOperationCreator_HH
#define INCLUDED_protocols_task_operations_pHVariantTaskOperationCreator_HH

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>

namespace protocols {
namespace task_operations {

class pHVariantTaskOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	core::pack::task::operation::TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


} //protocols
} //task_operations


#endif //INCLUDED_protocols/task_operations_pHVariantTaskOperationCreator_HH
