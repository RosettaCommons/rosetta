// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/SequenceMotifTaskOperation.cc
/// @brief A TaskOp that takes a regex-like pattern and turns it into a set of design residues.  The string should identify what to do for each position.  An X indicates any residue, and is the same as [A_Z].  Anything other than a charactor should be in [].  The ^ denotes a not.  An example Regex for glycosylation is NX[ST].  This would design an ASN into the first position, skip the second, and allow S and T mutations only at the third position.  N[^P][ST] would denote that at the second position, we do not allow proline.  Sets of charactors using _ can be denoted, even though this doesnt really help us in design.  In the future one can imagine having sets of polars, etc. etc.  This TaskOp is like a simple RESFILE in a string.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/task_operations/SequenceMotifTaskOperation.hh>
#include <protocols/task_operations/SequenceMotifTaskOperationCreator.hh>

#include <core/pack/task/util.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/ResfileReader.hh>

#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/sequence/sequence_motif.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/datacache/DataMap.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

static basic::Tracer TR( "protocols.task_operations.SequenceMotifTaskOperation" );

namespace protocols {
namespace task_operations {
using namespace core::select::residue_selector;
using namespace core::pack::task;
using namespace core::pack::task::operation;

SequenceMotifTaskOperation::SequenceMotifTaskOperation():
	TaskOperation()
{

}

SequenceMotifTaskOperation::SequenceMotifTaskOperation(
	ResidueSelectorCOP selector,
	std::string const & motif ):

	TaskOperation()
{
	set_residue_selector(selector);
	set_motif(motif);
}

SequenceMotifTaskOperation::~SequenceMotifTaskOperation() {}

core::pack::task::operation::TaskOperationOP
SequenceMotifTaskOperation::clone() const {
	return TaskOperationOP( new SequenceMotifTaskOperation( *this ) );
}

SequenceMotifTaskOperation::SequenceMotifTaskOperation( SequenceMotifTaskOperation const & src ):
	TaskOperation(src)

{
	motif_ = src.motif_;
	if ( motif_ != "" ) {
		setup_commands();
	}
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
}

void
SequenceMotifTaskOperation::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap& data){

	if ( tag->hasOption("motif") ) {
		set_motif( tag->getOption< std::string >("motif") );
	} else {
		utility_exit_with_message("SequenceMotifTaskOperation: Need required option 'motif'");
	}

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(parse_residue_selector(tag, data));
	} else {
		utility_exit_with_message("SequenceMotifTaskOperation: Need required option 'residue_selector'");
	}
}


void
SequenceMotifTaskOperation::set_motif(std::string const & motif){
	motif_ = motif;
	setup_commands();
}

void
SequenceMotifTaskOperation::setup_commands(){
	debug_assert( motif_ != "");

	commands_.clear();
	commands_ = get_resfile_commands( motif_ );

}


void
SequenceMotifTaskOperation::set_residue_selector( ResidueSelectorCOP selector ){
	selector_ = selector;
}

void
SequenceMotifTaskOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const {

	//First, we need a list of AA for each position.
	if ( selector_ == nullptr ) {
		utility_exit_with_message("SequenceMotifTaskOperation requires a residue selector!");
	}

	if ( motif_ == "" ) {
		utility_exit_with_message("SequenceMotifTaskOperation requires a motif to be set!");
	}


	utility::vector1< core::Size > subset = selector_->apply( pose );

	TR.Debug << "Creating motifs at " << utility::to_string( core::select::residue_selector::selection_positions(subset ));

	for ( core::Size i : selection_positions( subset ) ) {

		for ( core::Size motif_pos = 1; motif_pos <= commands_.size(); ++motif_pos ) {

			for ( ResfileCommandOP cmd : commands_[motif_pos] ) {
				cmd->residue_action(task, i + motif_pos -1 );
			}
		}
	}
}

std::string
SequenceMotifTaskOperation::keyname() {
	return "SequenceMotifTaskOperation";
}

void
SequenceMotifTaskOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	std::string motif_str = core::sequence::get_design_sequence_motif_syntax();
	attlist + XMLSchemaAttribute("motif", xs_string, motif_str);

	attributes_for_parse_residue_selector( attlist );

	core::pack::task::operation::task_op_schema_w_attributes( xsd, keyname(), attlist, "A TaskOp that takes a regex-like pattern and turns it into a set of design residues.  The string should identify what to do for each position.  Does not TURN OFF design or packing for ANY residue other than those specified in the motif as '-' or with specific resfile command!" );


}

core::pack::task::operation::TaskOperationOP
SequenceMotifTaskOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SequenceMotifTaskOperation );
}

std::string
SequenceMotifTaskOperationCreator::keyname() const
{
	return SequenceMotifTaskOperation::keyname();
}

void
SequenceMotifTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SequenceMotifTaskOperation::provide_xml_schema( xsd );


}


} //protocols
} //task_operations
