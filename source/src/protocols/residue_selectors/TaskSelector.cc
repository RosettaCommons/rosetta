// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/TaskSelector.hh
/// @brief  The TaskSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

// Unit headers
#include <protocols/residue_selectors/TaskSelector.hh>
#include <protocols/residue_selectors/TaskSelectorCreator.hh>

// Core headers
#include <core/select/residue_selector/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.residue_selectors.TaskSelector" );

namespace protocols {
namespace residue_selectors {

core::select::residue_selector::ResidueSelectorOP
TaskSelectorCreator::create_residue_selector() const
{
	return core::select::residue_selector::ResidueSelectorOP( new TaskSelector );
}

std::string
TaskSelectorCreator::keyname() const
{
	return TaskSelector::class_name();
}

void
TaskSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TaskSelector::provide_xml_schema( xsd );
}

TaskSelector::TaskSelector() :
	ResidueSelector(),
	tf_(),
	select_designable_( true ),
	select_packable_( true ),
	select_fixed_( false )
{}

TaskSelector::TaskSelector(
	core::pack::task::TaskFactoryOP tf,
	bool const select_designable,
	bool const select_packable,
	bool const select_fixed ):
	tf_(std::move( tf )),
	select_designable_( select_designable ),
	select_packable_( select_packable ),
	select_fixed_( select_fixed )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP
TaskSelector::clone() const
{
	return core::select::residue_selector::ResidueSelectorOP( new TaskSelector(*this) );
}

TaskSelector::~TaskSelector() = default;

void
quit_no_tf()
{
	std::stringstream msg;
	msg << "TaskSelector: no task factory was specified!" << std::endl;
	throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
}

core::select::residue_selector::ResidueSubset
TaskSelector::apply( core::pose::Pose const & pose ) const
{
	if ( !tf_ ) quit_no_tf();

	core::pack::task::PackerTaskOP task = tf_->create_task_and_apply_taskoperations( pose );
	debug_assert( task );
	TR.Debug << "PackerTask:" << std::endl;
	task->show( TR.Debug );
	TR.Debug.flush();
	core::select::residue_selector::ResidueSubset subset;
	for ( core::Size resid=1; resid<=task->total_residue(); ++resid ) {
		bool selected = false;
		if ( select_designable_ && task->being_designed( resid ) ) {
			selected = true;
		}
		if ( select_packable_ && task->being_packed( resid ) ) {
			selected = true;
		}
		if ( select_fixed_ && ( !task->being_packed( resid ) && !task->being_designed( resid ) ) ) {
			selected = true;
		}
		subset.push_back( selected );
	}
	return subset;
}
void
TaskSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	std::string const task_operations_str = tag->getOption< std::string >( "task_operations", "" );
	if ( task_operations_str.empty() ) quit_no_tf();

	utility::vector1< std::string > const task_op_names = utility::string_split( task_operations_str, ',' );
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory );
	for (const auto & task_op_name : task_op_names) {
		if ( data.has( "task_operations", task_op_name ) ) {
			tf->push_back( data.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", task_op_name ) );
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + task_op_name + " not found in basic::datacache::DataMap.");
		}
		TR << "Adding task operation " << task_op_name << std::endl;
	}
	set_task_factory( tf );

	set_select_designable( tag->getOption< bool >( "designable", select_designable_ ) );
	set_select_packable( tag->getOption< bool >( "packable", select_packable_ ) );
	if ( tag->hasOption( "repackable" ) ) {
		set_select_packable( tag->getOption< bool >( "repackable" ) );
	}
	set_select_fixed( tag->getOption< bool >( "fixed", select_fixed_ ) );
}

// APL TO DO!
void
TaskSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "task_operations", xs_string )
		+ XMLSchemaAttribute( "designable", xs_boolean )
		+ XMLSchemaAttribute( "packable", xs_boolean )
		+ XMLSchemaAttribute( "repackable", xs_boolean )
		+ XMLSchemaAttribute( "fixed", xs_boolean );

	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), attributes );

}

void
TaskSelector::set_task_factory( core::pack::task::TaskFactoryOP tf )
{
	tf_ = tf;
}

void
TaskSelector::set_select_designable( bool const sel_designable )
{
	select_designable_ = sel_designable;
}

void
TaskSelector::set_select_packable( bool const sel_packable )
{
	select_packable_ = sel_packable;
}

void
TaskSelector::set_select_fixed( bool const sel_fixed )
{
	select_fixed_ = sel_fixed;
}

std::string
TaskSelector::get_name() const
{
	return TaskSelector::class_name();
}

std::string
TaskSelector::class_name()
{
	return "Task";
}

} //namespace residue_selectors
} //namespace protocols

