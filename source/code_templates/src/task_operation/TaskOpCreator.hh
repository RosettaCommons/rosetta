// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   --path--/--class--.cc
/// @brief  Creator for --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--Creator_HH
#define INCLUDED_--path_underscore--_--class--Creator_HH

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>


--namespace--

class --class--Creator : public core::pack::task::operation::TaskOperationCreator {
public:
	core::pack::task::operation::TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};



--end_namespace--


#endif //INCLUDED_--path--_--class--Creator_HH
