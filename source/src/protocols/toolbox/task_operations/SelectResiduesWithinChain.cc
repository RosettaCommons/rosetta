// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/SelectResiduesWithinChainOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/SelectResiduesWithinChain.hh>
#include <protocols/toolbox/task_operations/SelectResiduesWithinChainCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <utility/string_util.hh>

// C++ Headers
#include <utility/vector1.hh>
#include <algorithm>

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.SelectResiduesWithinChainOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;

SelectResiduesWithinChainOperation::SelectResiduesWithinChainOperation() :
	chain_( 1 ),
	allow_design_( true )
{
	resid_.clear();
}

SelectResiduesWithinChainOperation::~SelectResiduesWithinChainOperation() {}

core::pack::task::operation::TaskOperationOP
SelectResiduesWithinChainOperationCreator::create_task_operation() const
{
	return new SelectResiduesWithinChainOperation;
}

core::pack::task::operation::TaskOperationOP SelectResiduesWithinChainOperation::clone() const
{
	return new SelectResiduesWithinChainOperation( *this );
}

///@brief restricts to repacking all residues outside of design_shell_ around each residue
void
SelectResiduesWithinChainOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	utility::vector1< core::Size > packing_residues, prevent_repacking_residues;
	packing_residues.clear(); prevent_repacking_residues.clear();

	if( allow_design() )
		TR<<"Residues set to design (all others are prevented from repacking): ";
	else
		TR<<"Residues set to repacking (all others are prevented from repacking): ";

	for( core::Size i = 1; i <= pose.total_residue(); ++i ){
		if( pose.residue( i ).chain() != chain() ){
			prevent_repacking_residues.push_back( i );
			continue;
		}
		core::Size const within_chain_idx( std::max( 0, (int)( i - pose.conformation().chain_begin( chain() ) + 1 )) );
		if( std::find( resid_.begin(), resid_.end(), within_chain_idx ) != resid_.end() ){
			TR<<i<<',';
			if( !allow_design() )
				packing_residues.push_back( i );
		}//fi std::find
		else
			prevent_repacking_residues.push_back( i );
	}//for i
	TR<<std::endl;
	OperateOnCertainResidues oocr_repacking, oocr_prevent_repacking;
	if( packing_residues.size() ){
		oocr_repacking.op( new RestrictToRepackingRLT );
		oocr_repacking.residue_indices( packing_residues );
		oocr_repacking.apply( pose, task );
	}
	if( prevent_repacking_residues.size() ){
		oocr_prevent_repacking.op( new PreventRepackingRLT );
		oocr_prevent_repacking.residue_indices( prevent_repacking_residues );
		oocr_prevent_repacking.apply( pose, task );
	}
}

void
SelectResiduesWithinChainOperation::parse_tag( TagPtr tag )
{
	chain( tag->getOption< core::Size >( "chain", 1 ) );
	std::string const res( tag->getOption< std::string >( "resid" ) );
	allow_design( tag->getOption< bool >( "allow_design", 1 ) );

	resid_ = utility::string_split< core::Size >( res, ',', core::Size() );
	TR<<"chain: "<<chain()<<" allow_design: "<<allow_design()<<" over residues: ";
	foreach( core::Size const r, resid() ){
		TR<<r<<", ";
	}
	TR<<std::endl;
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
