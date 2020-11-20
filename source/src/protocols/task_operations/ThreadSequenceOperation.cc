// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/ThreadSequenceOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/task_operations/ThreadSequenceOperation.hh>
#include <protocols/task_operations/ThreadSequenceOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <utility/vector0.hh>

// C++ headers
#include <vector>


static basic::Tracer TR( "protocols.TaskOperations.ThreadSequenceOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

core::pack::task::operation::TaskOperationOP
ThreadSequenceOperationCreator::create_task_operation() const
{
	return utility::pointer::make_shared< ThreadSequenceOperation >();
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
	target_sequence( "" ); start_res( 1 ); allow_design_around( true ); chain_num( 0 ); filter_non_aas( false );
}

ThreadSequenceOperation::ThreadSequenceOperation( std::string const & seq ) : parent()
{
	target_sequence( seq );
}

ThreadSequenceOperation::~ThreadSequenceOperation() = default;

core::pack::task::operation::TaskOperationOP ThreadSequenceOperation::clone() const
{
	return utility::pointer::make_shared< ThreadSequenceOperation >( *this );
}

void
ThreadSequenceOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using namespace core::chemical;

	// Since we can not assume the pose hasn't changed (although it is weird)
	// we must reselect residues everytime we apply
	utility::vector1<std::pair<core::Size, char>> thread_residues;
	core::pack::task::operation::RestrictResidueToRepacking rtr;

	core::Size const thread_begin = ( chain_num() == 0 ? start_res()
		: pose.conformation().chain_begin( chain_num() ) + start_res() - 1 );
	core::Size const chain_end = ( chain_num() == 0 ? pose.size() : pose.conformation().chain_end( chain_num() ));

	// Sanity checks - make sure we won't break anything down the road
	runtime_assert( target_sequence().length() -1 + thread_begin <= pose.size() );

	if ( filter_non_aas() ) {
		core::Size thread_offset = 0;
		core::Size sequence_idx = 0;
		while ( thread_begin + thread_offset < chain_end and sequence_idx < target_sequence_.length() ) {
			if ( !pose.residue_type(thread_begin + thread_offset ).is_protein() ) {
				rtr.include_residue(thread_begin + thread_offset );
				thread_offset++;
				continue;
			}
			thread_residues.emplace_back(thread_begin + thread_offset, target_sequence_[sequence_idx]);
			++thread_offset;
			++sequence_idx;
		}
	} else {
		for ( core::Size i = 0; i < target_sequence_.length(); ++i ) {
			thread_residues.emplace_back(thread_begin + i, target_sequence_[i]);
		}
	}

	core::Size const thread_end = thread_residues.back().first;
	// Sanity checks - make sure we won't break anything down the road
	runtime_assert( thread_end <= pose.size() );

	if ( !allow_design_around() ) {
		for ( core::Size i( 1 ); i < thread_begin; ++i ) {
			rtr.include_residue(i);
		}
		for ( core::Size i(thread_end + 1); i <= pose.size(); ++i ) {
			rtr.include_residue(i);
		}
	}

	core::Size resi;
	char aa;
	for ( auto & thread_residue : thread_residues ) {
		resi = thread_residue.first;
		aa = thread_residue.second;
		if ( aa == 'x' ) continue; // allows for 'wildcard' residues that can be allowed to design to anything within the threaded sequence
		if ( aa == ' ' or aa == '_' ) {
			/// Only repack this residue
			rtr.include_residue(resi);
			continue;
		}
		restrict_aas(resi, aa, pose, task );
	}
	rtr.apply( pose, task );
}

void
ThreadSequenceOperation::restrict_aas(core::Size const resi, char const aa, core::pose::Pose const & pose, core::pack::task::PackerTask & task) {
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using namespace core::chemical;

	RestrictAbsentCanonicalAAS racaas;
	utility::vector1< bool > keep_aas( num_canonical_aas, false );
	keep_aas[ aa_from_oneletter_code( aa ) ] = true;
	racaas.keep_aas( keep_aas );
	racaas.include_residue( resi );
	racaas.apply(pose, task);
}

void
ThreadSequenceOperation::target_sequence( std::string const & seq )
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
	chain_num( tag->getOption< core::Size >( "chain_num", 0 ) );
	filter_non_aas( tag->getOption< bool >( "filter_non_aas", false ) );
	TR<<"Threading chain: "<<chain_num()<<" with sequence: "<<target_sequence()<<" starting at residue #"<<start_res()<<" allow design around "<<allow_design_around()<<std::endl;
}

void ThreadSequenceOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::required_attribute( "target_sequence", xs_string , "The target sequence can contain two types of 'wildcards'. Placing 'x' in the sequence results in design at this position: target_sequence='TFYxxxHFS' will thread the two specified tripeptides and allow design in the intervening tripeptide. Placing ' ' (space) or '_' (underscore), however, restricts this position to repacking: the string 'TFY HFS' (three spaces between the two triplets) will thread the two tripeptides and will repack the pose's original intervening tripeptide. The string 'TFY___HFS' -three underscores between the two triplets- will also only repack the original intervening tripeptide." )
		+ XMLSchemaAttribute::attribute_w_default(  "start_res", xsct_positive_integer, "Residue at which to start." , "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "allow_design_around", xsct_rosetta_bool, "If set to false, only design the region that is threaded. The rest is set to repack.",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "chain_num", xsct_non_negative_integer, "Which chain to apply the sequence to. If 0 - the entire structure",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "filter_non_aas", xsct_rosetta_bool, "If set to true non amino acid residues will not be threaded",  "false"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "Threads a single letter sequence onto the source pdb." );
}

} //namespace protocols
} //namespace task_operations
