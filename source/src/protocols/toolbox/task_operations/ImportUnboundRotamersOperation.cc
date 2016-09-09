// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/ImportUnboundRotamersOperation.cc
/// @brief  Import unbound (or native) rotamers into Rosetta!
/// @author Dave La ( davela@uw.edu )


// Unit Headers
#include <protocols/toolbox/task_operations/ImportUnboundRotamersOperation.hh>
#include <protocols/toolbox/task_operations/ImportUnboundRotamersOperationCreator.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <utility/tag/Tag.hh>


#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#ifdef WIN32
#include <utility/graph/Graph.hh>
#endif

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

///////////////////////////////////////////////////////////////////////////////////////////
TaskOperationOP
ImportUnboundRotamersOperationCreator::create_task_operation() const
{
	return TaskOperationOP( new ImportUnboundRotamersOperation );
}

void ImportUnboundRotamersOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ImportUnboundRotamersOperation::provide_xml_schema( xsd );
}

std::string ImportUnboundRotamersOperationCreator::keyname() const
{
	return ImportUnboundRotamersOperation::keyname();
}

/// @brief default constructor
ImportUnboundRotamersOperation::ImportUnboundRotamersOperation():
	TaskOperation()
{}

/// @brief destructor
ImportUnboundRotamersOperation::~ImportUnboundRotamersOperation(){}

/// @brief clone
TaskOperationOP
ImportUnboundRotamersOperation::clone() const {
	return TaskOperationOP( new ImportUnboundRotamersOperation( *this ) );
}

//mjo commenting out 'pose' because it is unused and causes a warning
/// @brief
void
ImportUnboundRotamersOperation::apply( Pose const & /*pose*/, PackerTask & task ) const
{
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot_( new core::pack::rotamer_set::UnboundRotamersOperation );
	unboundrot_->initialize_from_command_line();
	task.append_rotamerset_operation( unboundrot_ );

}

void
ImportUnboundRotamersOperation::parse_tag( TagCOP /*tag*/ , DataMap & )
{
}

void ImportUnboundRotamersOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

} // TaskOperations
} // toolbox
} // protocols

