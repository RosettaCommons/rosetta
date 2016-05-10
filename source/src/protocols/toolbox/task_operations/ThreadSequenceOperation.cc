// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/ThreadSequenceOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/operation/TaskOperations.hh>

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
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.ThreadSequenceOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

core::pack::task::operation::TaskOperationOP
ThreadSequenceOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ThreadSequenceOperation );
}

void ThreadSequenceOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ThreadSequenceOperation::provide_xml_schema( xsd );
}

std::string ThreadSequenceOperationCreator::keyname() const
{
	return ThreadSequenceOperation::keyname();
}

ThreadSequenceOperation::ThreadSequenceOperation() : parent()
{
	target_sequence( "" ); start_res( 1 ); allow_design_around( true );
}

ThreadSequenceOperation::ThreadSequenceOperation( std::string const seq ) : parent()
{
	target_sequence( seq );
}

ThreadSequenceOperation::~ThreadSequenceOperation() {}

core::pack::task::operation::TaskOperationOP ThreadSequenceOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ThreadSequenceOperation( *this ) );
}

void
ThreadSequenceOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using namespace core::chemical;

	runtime_assert( target_sequence().length() + start_res() -1 <= pose.total_residue() );
	core::pack::task::operation::RestrictResidueToRepacking rtr;
	bool activated_rtr( false );
	for ( core::Size resi( 1 ); resi <= pose.total_residue(); ++resi ) {
		if ( resi >= start_res() && resi <= start_res() + target_sequence_.length() - 1 ) {
			if ( target_sequence_[ resi - start_res() ] == 'x' ) continue; // allows for 'wildcard' residues that can be allowed to design to anything within the threaded sequence
			if ( target_sequence_[ resi - start_res() ] == ' ' || target_sequence_[ resi - start_res() ] == '_' ) { /// only repack this residue
				rtr.include_residue( resi );
				activated_rtr = true;
				continue;
			}
			RestrictAbsentCanonicalAAS racaas;
			utility::vector1< bool > keep_aas( num_canonical_aas, false );
			keep_aas[ aa_from_oneletter_code( target_sequence_[ resi-start_res() ] ) ] = true;
			racaas.keep_aas( keep_aas );
			racaas.include_residue( resi );
			racaas.apply( pose, task );
		} else if ( !allow_design_around() ) {
			rtr.include_residue( resi );
		}
	}
	if ( !allow_design_around() || activated_rtr ) {
		rtr.apply( pose, task );
	}
}

void
ThreadSequenceOperation::target_sequence( std::string const seq )
{
	target_sequence_ = seq;
}

std::string
ThreadSequenceOperation::target_sequence() const
{
	return( target_sequence_ );
}

core::Size
ThreadSequenceOperation::start_res() const{
	return start_res_;
}

void
ThreadSequenceOperation::start_res( core::Size const s ){
	start_res_ = s;
}

void
ThreadSequenceOperation::parse_tag( TagCOP tag , DataMap & )
{
	target_sequence( tag->getOption< std::string >( "target_sequence" ) );
	start_res( tag->getOption< core::Size >( "start_res", 1 ) );
	allow_design_around( tag->getOption< bool >( "allow_design_around", true ) );
	TR<<"Threading with sequence: "<<target_sequence()<<" starting at residue #"<<start_res()<<" allow design around "<<allow_design_around()<<std::endl;
}

void ThreadSequenceOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::required_attribute( "target_sequence", xs_string )
		+ XMLSchemaAttribute::attribute_w_default(  "start_res", xsct_non_negative_integer, "1" )
		+ XMLSchemaAttribute::attribute_w_default(  "allow_design_around", xs_boolean, "true" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
