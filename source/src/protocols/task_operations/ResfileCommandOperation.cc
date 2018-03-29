// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/ResfileCommandOperation.cc
/// @brief Applies the equivalent of a resfile line (without the resnums) to residues specified in a residue selector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/task_operations/ResfileCommandOperation.hh>
#include <protocols/task_operations/ResfileCommandOperationCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pack/task/operation/task_op_schemas.hh>



static basic::Tracer TR( "protocols.task_operations.ResfileCommandOperation" );

namespace protocols {
namespace task_operations {
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace core::select::residue_selector;

ResfileCommandOperation::ResfileCommandOperation():
	TaskOperation()

{
	command_map_ = create_command_map();;
}

ResfileCommandOperation::ResfileCommandOperation( ResidueSelectorCOP selector, std::string const & command):
	TaskOperation()

{
	command_map_ = create_command_map();
	set_residue_selector(selector);
	set_command(command);
}

ResfileCommandOperation::~ResfileCommandOperation() {}

TaskOperationOP
ResfileCommandOperation::clone() const {
	return TaskOperationOP( new ResfileCommandOperation( *this ) );
}

ResfileCommandOperation::ResfileCommandOperation( ResfileCommandOperation const & src ):
	TaskOperation(src)

{
	command_map_ = create_command_map();

	if ( src.selector_ ) {
		set_residue_selector(src.selector_);
	}

	set_command(src.command_);
}

void
ResfileCommandOperation::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap& data){

	set_command(tag->getOption< std::string >("command"));
	set_residue_selector(protocols::rosetta_scripts::parse_residue_selector(tag, data));
}


void
ResfileCommandOperation::set_command( std::string const & command ){
	command_ = command;
	commands_.clear();

	commands_ = parse_res_agnostic_commands( command_, command_map_ );
}

void
ResfileCommandOperation::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector;
}

void
ResfileCommandOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task) const {

	if ( (selector_ == nullptr) | ( command_ == "" ) ) {
		utility_exit_with_message("residue_selector and command must be set!");
	}

	utility::vector1< bool > subset = selector_->apply( pose );
	for ( core::Size const & position : selection_positions( subset ) ) {
		for ( ResfileCommandOP cmd : commands_ ) {
			cmd->residue_action(task, position );
		}
	}

}

std::string
ResfileCommandOperation::keyname() {
	return "ResfileCommandOperation";
}


void
ResfileCommandOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	std::string command_doc =  "@brief A resfile command string without any numbers in the front.\n"
		"Example:\n"
		" POLAR\n"
		" EMPTY NC R2 NC T6 NC OP5 \n";

	attlist + XMLSchemaAttribute("command", xs_string, command_doc);

	attributes_for_parse_residue_selector( attlist );
	task_op_schema_w_attributes( xsd, keyname(), attlist, "Applies the equivalent of a resfile line (without the resnums) to residues specified in a residue selector." );
}

TaskOperationOP
ResfileCommandOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ResfileCommandOperation );
}

std::string
ResfileCommandOperationCreator::keyname() const
{
	return ResfileCommandOperation::keyname();
}

void
ResfileCommandOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResfileCommandOperation::provide_xml_schema( xsd );
}


} //protocols
} //task_operations
