// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_toolbox_task_operations_ProteinCoreResFilter_Creator_hh
#define INCLUDED_protocols_toolbox_task_operations_ProteinCoreResFilter_Creator_hh

#include <core/pack/task/operation/ResFilterCreator.hh>

#include <string>


namespace protocols {
namespace toolbox {
namespace task_operations {

class ProteinCoreFilterCreator : public core::pack::task::operation::ResFilterCreator {
public:
	virtual core::pack::task::operation::ResFilterOP create_res_filter() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif

