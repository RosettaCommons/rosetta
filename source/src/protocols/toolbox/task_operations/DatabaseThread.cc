// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/DatabaseThread.cc
/// @brief  picks a sequence from database by start and end position on the pose
/// @author Assaf Alon assafalon@gmail.com

// Unit Headers
#include <protocols/toolbox/task_operations/DatabaseThread.hh>
#include <protocols/toolbox/task_operations/DatabaseThreadCreator.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>

//Package headers
#include <utility/io/izstream.hh>
#include <core/pose/Pose.hh>
#include <utility/string_util.hh>

#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/random/random.hh>
#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.DatabaseThread" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

core::pack::task::operation::TaskOperationOP
DatabaseThreadCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DatabaseThread );
}

void DatabaseThreadCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DatabaseThread::provide_xml_schema( xsd );
}

std::string DatabaseThreadCreator::keyname() const
{
	return DatabaseThread::keyname();
}

DatabaseThread::DatabaseThread() : parent(),
	target_sequence_(""),
	template_file_(""),
	database_fname_( "" ),
	template_pose_(/* NULL */),
	start_res_(0),
	end_res_(0 ),
	allow_design_around_(true )
{
	design_.clear(); revert_to_template_.clear(); full_database_.clear(); designable_.clear(); leave_as_is_.clear();
}

DatabaseThread::~DatabaseThread() {}

core::pack::task::operation::TaskOperationOP DatabaseThread::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new DatabaseThread( *this ) );
}


void
DatabaseThread::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	protocols::toolbox::task_operations::ThreadSequenceOperation thread;
	std::string sequence;
	if ( target_sequence()=="" ) {
		sequence=pick_sequence_from_database(pose);
	} else {
		sequence=target_sequence();
	}
	core::Size const nearest_to_start_on_pose( rosetta_scripts::find_nearest_res( pose, *template_pose_, start_res(), 1/*chain*/));
	if ( !designable_.empty() ) { //if user entered positions to mark for design, label those with 'X'
		mark_designable(sequence,pose);
	}
	if ( !leave_as_is_.empty() ) { //if user entered positions to avoid threading mark those with '_'
		mark_leave_as_is(sequence, pose);
	}
	thread.start_res(nearest_to_start_on_pose);
	thread.allow_design_around(allow_design_around_);
	TR<<"Threading the following sequence:"<<std::endl;
	TR<<sequence<<std::endl;
	thread.target_sequence(sequence);
	thread.apply(pose,task);
}

core::Size
DatabaseThread::find_length(const core::pose::Pose &pose) const // find the lenght of the pose according to template positions - if there was an insertion or deletion along the way.
{
	core::Size const user_start(start_res()); //get the user defined start residue on the template
	core::Size const nearest_to_start_on_pose( rosetta_scripts::find_nearest_res( pose, *template_pose_, user_start, 1/*chain*/ ) );
	core::Size const user_end(end_res()); //get the user defined end residue on the template
	core::Size const nearest_to_end_on_pose( rosetta_scripts::find_nearest_res( pose, *template_pose_, user_end, 1/*chain*/ ) );
	int const start_end_delta( nearest_to_end_on_pose - nearest_to_start_on_pose + 1);
	if ( start_end_delta<0 ) {
		utility_exit_with_message("Difference between start and end residues is negative - Aborting!!!");
	} else {
		TR<<"The position on the pose to start is "<<nearest_to_start_on_pose<<std::endl;
		TR<<"The position on the pose to end is "<<nearest_to_end_on_pose<<std::endl;
		return start_end_delta;
	}
}


std::string
DatabaseThread::pick_sequence_from_database( core::pose::Pose const & pose ) const{
	core::Size const segment_length(find_length( pose ));
	//std::string line;
	utility::vector1< std::string > sized_database;
	for ( std::string const & line : full_database_ ) {
		if ( line.length()==segment_length ) { // if length of line is the same as segment length, incorporate into vector of strings.
			sized_database.push_back( line );
		}
	}
	if ( sized_database.size()==0 ) {
		utility_exit_with_message("no entries with correct length were found in the database: " + database_fname_ );
	}
	TR<<"Finished reading database "<<database_fname_<<" with "<<sized_database.size()<<" entries of length "<<segment_length<<std::endl;
	core::Size const entry = numeric::random::rg().uniform() * sized_database.size() + 1;
	TR<<"Picked the sequence:"<<std::endl;
	TR<<sized_database[entry]<<std::endl;
	return sized_database[entry];
}

void //mark residue positions on the pose (not positions in the sequence!) as designlable
DatabaseThread::mark_designable(std::string & sequence, core::pose::Pose const & pose ) const{
	for ( core::Size const i : designable_ ) {
		core::Size const near_to_i(rosetta_scripts::find_nearest_res(pose,*template_pose_,i,1/*chain*/));
		core::Size const position(near_to_i - start_res_ +1);
		sequence[position]='X';
	}
}


void //mark residue positions on the pose (not positions in the sequence!) to leave as in template
DatabaseThread::mark_leave_as_is(std::string &sequence, core::pose::Pose const & pose ) const{
	for ( core::Size const i : leave_as_is_ ) {
		core::Size const near_to_i(rosetta_scripts::find_nearest_res(pose,*template_pose_,i,1/*chain*/));
		core::Size const position(near_to_i - start_res_ +1);
		sequence[position]='_';
	}
}

void
DatabaseThread::parse_tag(TagCOP tag, DataMap &)
{
	target_sequence( tag->getOption< std::string >( "target_sequence","" ) );
	template_file( tag->getOption< std::string >( "template_file") );
	template_pose_ = core::pose::PoseOP( new core::pose::Pose );
	core::import_pose::pose_from_file( *template_pose_, template_file_ , core::import_pose::PDB_file);
	database_fname( tag->getOption< std::string >( "database","" ) );
	if ( target_sequence()=="" ) {
		if ( database_fname()=="" ) {
			utility_exit_with_message("Please provide either a database file or target sequence! Aborting!");
		} else {
			utility::io::izstream database( database_fname_ );
			if ( !database ) {
				utility_exit_with_message("cannot open database " + database_fname_ +"\n");
			}
			std::string line;
			while ( getline( database, line ) ) { // if length of line is the same as segment length, incorporate into vector of strings.
				full_database_.push_back( line );
			}
		}
	}
	start_res( tag->getOption< core::Size >( "start_res" ) );
	end_res( tag->getOption< core::Size >( "end_res" ) );
	allow_design_around( tag->getOption< bool >( "allow_design_around", true ) );
	designable(utility::string_split(tag->getOption<std::string>("design_residues",""),',',core::Size()));
	leave_as_is(utility::string_split(tag->getOption<std::string>("keep_original_identity",""),',',core::Size()));
}

// AMW: Comma separated string list... candidate for common_simple_types?
void DatabaseThread::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"target_sequence", xs_string,
		"The desired sequence if there is only one desired sequence "
		"(this can happen if the pose is changed during design such that the start "
		"and end positions are not constant. In such cases ThreadSequence is not useful). "
		"The task operation expects either a database or a target sequence and will "
		"fail if neither are provided. If both are provided the database will be ignored.",
		"")
		+ XMLSchemaAttribute::required_attribute(
		"template_file", xs_string,
		"a pdb that serves as a constant template to map the start and end residues onto "
		"the pose in case that the length of the pose is altered during design" )
		+ XMLSchemaAttribute::attribute_w_default(
		"database", xs_string,
		"The database should be a text file with a list of single letter amino acids (not fasta)",
		""  )

		+ XMLSchemaAttribute::required_attribute(
		"start_res", xsct_non_negative_integer,
		"the residue to start threading from. This is a residue in the template pdb. "
		"It is used to find the closest residue on the source pdb." )
		+ XMLSchemaAttribute::required_attribute(
		"end_res", xsct_non_negative_integer,
		"the residue to end the threading. This is a residue in the template pdb. "
		"It is used to find the closest residue on the source pdb. "
		"The delta between the end and start residue is used to find the desired "
		"sequence length in the database." )
		+ XMLSchemaAttribute::attribute_w_default(
		"allow_design_around", xsct_rosetta_bool,
		"if set to false, only design the region that is threaded. The rest is set to repack.",
		"true"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"design_residues", xs_string,
		"the same as placing 'X' in the target sequence. This trumps the sequence in the "
		"database so if a residue has a different identity in the database it is changed to 'X'.",
		""  )
		+ XMLSchemaAttribute::attribute_w_default(
		"keep_original_identity", xs_string,
		"the same as placing a ' ' or a '_' in the sequence. The pose residue is marked "
		"for repacking only. This trumps both the database sequence and the list from design_residues.",
		""  );

	task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"This task operation is designed to deal with situations in which a pose is "
		"changed in a way that adds or removes residues. This creates a problem for "
		"normal threading that requires a constant start and stop positions. This "
		"task operation can use a database of sequences or a single target sequence. "
		"It also need a template pdb to find on the pose a user defined start and end residues. "
		"A sequence length and threading start position are calculated and then a "
		"correct length sequence is randomly chosen from the "
		"database and threaded onto the pose.");
}


} //namespace protocols
} //namespace toolbox
} //namespace task_operations
