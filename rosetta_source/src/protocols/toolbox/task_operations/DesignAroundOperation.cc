// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/DesignAroundOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperationCreator.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project Headers
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
// AUTO-REMOVED #include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
// AUTO-REMOVED #include <core/pack/task/operation/NoRepackDisulfides.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictToInterface.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>
// Auto-header: duplicate removed #include <core/pack/task/operation/TaskOperations.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <set>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.DesignAroundOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace std;

DesignAroundOperation::DesignAroundOperation() :
	design_shell_( 8.0 ),
	string_resnums_( "" ),
	repack_on_( true )
{
	resid_.clear();
}

DesignAroundOperation::~DesignAroundOperation() {}

core::pack::task::operation::TaskOperationOP
DesignAroundOperationCreator::create_task_operation() const
{
	return new DesignAroundOperation;
}

core::pack::task::operation::TaskOperationOP DesignAroundOperation::clone() const
{
	return new DesignAroundOperation( *this );
}

///@brief restricts to repacking all residues outside of design_shell_ around each residue
void
DesignAroundOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task::operation;

	set< core::Size > focus_residues;// all of the residues that were input (notice that the method is const, so I can't change resid_)
	focus_residues.clear();
	focus_residues.insert( resid_.begin(), resid_.end() );
	set< core::Size > const res_vec( protocols::rosetta_scripts::get_resnum_list( string_resnums_, pose ) );
	focus_residues.insert( res_vec.begin(), res_vec.end() );

	TR.Debug<<"Design will be allowed around the following residues (others will be allowed to repack only): ";
	foreach( core::Size const res, focus_residues )
		TR.Debug<<res<<", ";
	TR.Debug<<std::endl;

	utility::vector1< core::Size > residues;
	residues.clear();
	for( core::Size i=1; i<=pose.total_residue(); ++i ){
		bool allow_design( false );
		foreach( core::Size const res, focus_residues ){
			core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )) );
			if( distance <= design_shell_ ){
				allow_design = true;
				break;
			}// fi distance
		}//foreach res
		if( !allow_design )
			residues.push_back( i );
	}//for i

	if( residues.size() ){
		TR.Debug<<"The following residues will be repacked only: ";
		foreach( core::Size const res, residues )
			TR.Debug<<res<<", ";
		TR.Debug<<std::endl;
		OperateOnCertainResidues oocr;
		if( repack_on() )
			oocr.op( new RestrictToRepackingRLT );
		else
			oocr.op( new PreventRepackingRLT );
		oocr.residue_indices( residues );
		oocr.apply( pose, task );
	}
}

void
DesignAroundOperation::design_shell( core::Real const radius )
{
	design_shell_ = radius;
}

void
DesignAroundOperation::include_residue( core::Size const resid )
{
	resid_.insert( resid );
}

void
DesignAroundOperation::parse_tag( TagPtr tag )
{
	string_resnums_ = tag->getOption< std::string >( "resnums" );// these are kept in memory until the pose is available (at apply time)
  design_shell( tag->getOption< core::Real >( "design_shell", 8.0 ) );
	repack_on( tag->getOption< bool >( "repack_on", 1 ) );
}

void
DesignAroundOperation::repack_on( bool const repack_on )
{
	repack_on_ = repack_on;
}

bool
DesignAroundOperation::repack_on() const
{
	return repack_on_;
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
