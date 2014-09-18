// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/OperateOnCertainResidues.cc
/// @brief  class to support the general case of configuring PackerTask at the ResidueLevelTask level in some way, with an optional filter to limit the effects to certain residues
/// @author ashworth

// Unit Headers
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OperateOnCertainResiduesCreator.hh>


#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/ResFilter.hh>
#include <core/pack/task/operation/ResFilterFactory.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

using basic::t_info;
using basic::t_debug;
static thread_local basic::Tracer TR( "core.pack.task.operation.OperateOnCertainResidues", t_info );

OperateOnCertainResidues::OperateOnCertainResidues()
	: parent(),
		op_(0),
		filter_(0)
{}

OperateOnCertainResidues::OperateOnCertainResidues(
	ResLvlTaskOperationOP rlto,
	ResFilterOP filter
)
	: parent(),
		op_( rlto ),
		filter_( filter )
{}

OperateOnCertainResidues::OperateOnCertainResidues( OperateOnCertainResidues const & src )
	: TaskOperation( src )
{
	*this = src;
}

OperateOnCertainResidues &
OperateOnCertainResidues::operator = ( OperateOnCertainResidues const & src )
{
	residue_indices_ = src.residue_indices();
	if ( src.op_ ) op_ = src.op_->clone();
	else op_ = 0;
	if ( src.filter_ ) filter_ = src.filter_->clone();
	else filter_ = 0;
	return *this;
}

OperateOnCertainResidues::~OperateOnCertainResidues() {}

TaskOperationOP OperateOnCertainResiduesCreator::create_task_operation() const
{
	return new OperateOnCertainResidues;
}

TaskOperationOP OperateOnCertainResidues::clone() const
{
	return new OperateOnCertainResidues( *this );
}

void
OperateOnCertainResidues::apply( Pose const & pose, PackerTask & ptask ) const
{
	Size const nres( pose.total_residue() );
	runtime_assert( nres == ptask.total_residue() );
	// local nonconst version of residue indices container allows default applicability over
	// all residues in input pose (while avoiding code duplication or extra function call)
	ResidueIndices residue_indices_local;
	if ( residue_indices_.empty() ) {
		for ( Size i(1); i <= nres; ++i ) residue_indices_local.push_back(i);
	} else {
		residue_indices_local = residue_indices_;
	}
	for ( ResidueIndices::const_iterator index( residue_indices_local.begin() ),
		end( residue_indices_local.end() ); index != end; ++index ) {
		runtime_assert( *index > 0 && *index <= nres );
		// skip this residue if there is a filter and it returns false
		if ( filter_ && ! (*filter_)( pose, *index ) ) continue;
		op_->apply( ptask.nonconst_residue_task( *index ) );
	}
}

void OperateOnCertainResidues::residue_indices( utility::vector1< Size > const & indices )
{
	residue_indices_ = indices;
}

void OperateOnCertainResidues::op( ResLvlTaskOperationCOP op_in )
{
	runtime_assert( op_in );
	op_ = op_in->clone();
}

void OperateOnCertainResidues::filter( ResFilterCOP filter_in )
{
	runtime_assert( filter_in );
	filter_ = filter_in->clone();
}


///@brief tag parsing for factory construction of this class and its children
/*!
Example Tag syntax for parser as of Summer 2009

<OperateOnCertainResidues name=PROTEINnopack>
	<PreventRepackingRLT/>
	<ResidueHasProperty property=PROTEIN/>
</OperateOnCertainResidues>

*/
void OperateOnCertainResidues::parse_tag( TagCOP tag , DataMap & )
{
  utility::vector0< TagCOP > const & subtags( tag->getTags() );
  for ( utility::vector0< TagCOP >::const_iterator subtag( subtags.begin() ), end( subtags.end() );
		subtag != end; ++subtag ) {
		std::string const type( (*subtag)->getName() );
		ResLvlTaskOperationFactory * rltof = ResLvlTaskOperationFactory::get_instance();
		if ( rltof && rltof->has_type( type ) ) {
			op_ = rltof->newRLTO( type );
			op_->parse_tag( *subtag );
			TR(t_debug) << "using ResLvlTaskOperation of type " << type << std::endl;
			continue;
		}
		ResFilterFactory * res_filter_factory = ResFilterFactory::get_instance();
		if ( res_filter_factory && res_filter_factory->has_type( type ) ) {
			filter_ = res_filter_factory->newResFilter( type );
			filter_->parse_tag( *subtag );
			TR(t_debug) << "using ResFilter of type " << type << std::endl;
			continue;
		}
		utility_exit_with_message( type + " is not known to the factories passed to OperateOnCertainResidues." );
	}
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
