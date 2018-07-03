// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ClashBasedRepackShell.cc
/// @brief  ClashBasedRepackShell header file.
/// @author Kale Kundert (kale@thekunderts.net)

// Unit headers
#include <core/pack/task/operation/ClashBasedRepackShell.hh>
#include <core/pack/task/operation/ClashBasedRepackShellCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelector.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>


using namespace std;
using core::Size;
using core::Real;
using core::select::residue_selector::ResidueSubset;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using core::pack::task::residue_selector::ClashBasedShellSelector;
using core::pack::task::residue_selector::ClashBasedShellSelectorOP;

namespace core {
namespace pack {
namespace task {
namespace operation {

static basic::Tracer TR( "core.pack.task.operation.ClashBasedRepackShell" );

ClashBasedRepackShell::ClashBasedRepackShell()
: shell_selector_( ClashBasedShellSelectorOP( new ClashBasedShellSelector ) )
{}

TaskOperationOP ClashBasedRepackShell::clone() const {
	return TaskOperationOP( new ClashBasedRepackShell( *this ) );
}

void ClashBasedRepackShell::apply( Pose const & pose, PackerTask & task ) const {
	ResidueSubset shell = shell_selector_->apply(pose);
	// Freeze any position that's not part of the shell.
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( shell[i] == false ) {
			task.nonconst_residue_task( i ).prevent_repacking();
		}
	}
}

void ClashBasedRepackShell::parse_tag(
	utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ) {

	shell_selector_->parse_my_tag(tag, datamap);
}

void
ClashBasedRepackShell::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd ) {

	utility::tag::AttributeList attributes;
	ClashBasedShellSelector::provide_xml_schema_attributes(attributes);

	task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"The ClashBasedShellSelector identifies all residues that clash with at least one rotamer of a given residue selection.");
}

string ClashBasedRepackShell::keyname() {
	return "ClashBasedRepackShell";
}

ClashBasedShellSelectorOP ClashBasedRepackShell::selector() const {
	return shell_selector_;
}

void ClashBasedRepackShell::selector(ClashBasedShellSelectorOP selector) {
	shell_selector_ = selector;
}

TaskOperationOP ClashBasedRepackShellCreator::create_task_operation() const {
	return TaskOperationOP( new ClashBasedRepackShell );
}

std::string ClashBasedRepackShellCreator::keyname() const {
	return ClashBasedRepackShell::keyname();
}

void ClashBasedRepackShellCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ClashBasedRepackShell::provide_xml_schema( xsd );
}


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core

