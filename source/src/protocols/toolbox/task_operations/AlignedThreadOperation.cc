// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.AlignedThreadOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
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

void AlignedThreadOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignedThreadOperation::provide_xml_schema( xsd );
}

std::string AlignedThreadOperationCreator::keyname() const
{
	return AlignedThreadOperation::keyname();
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
	while ( data ) {
		if ( line.length() == 0 ) {
			continue;
		}
		if ( line.substr(1, query_name().length() ) == query_name() ) {
			while ( data ) {
				getline( data, line );
				if ( line[0] == '>' ) {
					break;
				}
				query_seq += line;
			}
		} else if ( line.substr(1, template_name().length() ) == template_name() ) {
			while ( data ) {
				getline( data, line );
				if ( line[0] == '>' ) {
					break;
				}
				template_seq += line;
			}
		} else {
			while ( data ) {
				getline( data, line );
				if ( line[0] == '>' ) {
					break;
				}
			}
		}
	}
	data.close();
	TR<<"template seq:\n"<<template_seq<<"\nquery seq:\n"<<query_seq<<std::endl;
	ThreadSequenceOperation tso;
	std::string target_sequence("");
	tso.start_res( start_res() );

	for ( core::Size i = 0; i < template_seq.length(); ++i ) {
		if ( template_seq[i] == '-' ) {
			continue;
		}
		if ( query_seq[i] == '-' ) {
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

void AlignedThreadOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::required_attribute( "alignment_file", xs_string , "The name of the alignment file in FASTA format. Should be in the usual -'[right-pointing-angle-bracket]name_of_sequence' followed by the amino acid single letter sequence on the next line or lines- for this to work." )
		+ XMLSchemaAttribute::required_attribute( "query_name", xs_string , "The name of the query sequence, as written in the alignment file." )
		+ XMLSchemaAttribute::required_attribute( "template_name", xs_string , "The name of the template sequence, as written in the alignment file. the same sequence as that of the structure passed with -s." )
		+ XMLSchemaAttribute::attribute_w_default(  "start_res", xsct_positive_integer, "The residue at which to start threading. Useful for threading the non-first chain.",  "1"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "A task operation that enables threading of aligned residues between a query and a template. receives a FASTA format sequence alignment (file may hold multiple sequences), and allows the threading only of residues that are aligned between query and structure. positions where either the template structure or the query sequence have a gap '-' are skipped. suitable for when you wish to model a sequence over a structure, and they are of different lengths" );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
