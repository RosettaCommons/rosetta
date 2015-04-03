// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/AlignedThreadOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/AlignedThreadOperation.hh>
#include <protocols/toolbox/task_operations/AlignedThreadOperationCreator.hh>

#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/io/izstream.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.AlignedThreadOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace std;

AlignedThreadOperation::AlignedThreadOperation() :
	parent(),
	alignment_file_( "" ),
	query_name_( "" ),
	template_name_( "" ),
	start_res_( 1 )
{}

AlignedThreadOperation::~AlignedThreadOperation() {}

core::pack::task::operation::TaskOperationOP
AlignedThreadOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new AlignedThreadOperation );
}

core::pack::task::operation::TaskOperationOP AlignedThreadOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new AlignedThreadOperation( *this ) );
}

void
AlignedThreadOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	std::string query_seq, template_seq; // query is the sequence we want to model; template is the sequence of the PDB structure

  utility::io::izstream data(alignment_file());
  runtime_assert( data );
  string line;
  getline(data, line);
  while( data ) {
    if (line.length() == 0) {
			continue;
    }
		if( line.substr(1, query_name().length() ) == query_name() ){
			while( data ){
				getline( data, line );
				if( line[0] == '>' )
					break;
				query_seq += line;
			}
		}
		else if( line.substr(1, template_name().length() ) == template_name() ){
			while( data ){
				getline( data, line );
				if( line[0] == '>' )
					break;
				template_seq += line;
			}
		}
		else {
			while( data ){
				getline( data, line );
				if( line[0] == '>' )
					break;
			}
		}
	}
	data.close();
	TR<<"template seq:\n"<<template_seq<<"\nquery seq:\n"<<query_seq<<std::endl;
	ThreadSequenceOperation tso;
	std::string target_sequence("");
	tso.start_res( start_res() );

	for( core::Size i = 0; i < template_seq.length(); ++i ){
		if( template_seq[i] == '-' )
			continue;
		if( query_seq[i] == '-' ){
			target_sequence += template_seq[i];
			continue;
		}
		target_sequence += query_seq[i];
	}
	TR<<"sequence for threading: \n"<<target_sequence<<std::endl;
	tso.target_sequence( target_sequence );
	tso.apply( pose, task );
}

void
AlignedThreadOperation::parse_tag( TagCOP tag , DataMap & )
{
	alignment_file( tag->getOption< std::string >("alignment_file" ) );
	query_name( tag->getOption< std::string >("query_name" ));
	template_name( tag->getOption< std::string >("template_name" ));
	start_res( tag->getOption< core::Size >( "start_res", 1 ) );

	TR<<"Aligned thread with options: alignment_file: "<<alignment_file()<<" query_name: "<<query_name()<<" start_res: "<<start_res()<<std::endl;
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
