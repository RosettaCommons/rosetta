// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.cc
/// @brief
/// @author Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperationCreator.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <boost/foreach.hpp>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
// AUTO-REMOVED #include <set>

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>



using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.PreventResiduesFromRepackingOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;

PreventResiduesFromRepackingOperation::PreventResiduesFromRepackingOperation() {}

PreventResiduesFromRepackingOperation::PreventResiduesFromRepackingOperation( utility::vector1 < core::Size > residues )
      : parent(), residues_( residues )
{
}
	
PreventResiduesFromRepackingOperation::~PreventResiduesFromRepackingOperation() {}

core::pack::task::operation::TaskOperationOP
PreventResiduesFromRepackingOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new PreventResiduesFromRepackingOperation );
}

core::pack::task::operation::TaskOperationOP PreventResiduesFromRepackingOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new PreventResiduesFromRepackingOperation( *this ) );
}

void
PreventResiduesFromRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	runtime_assert( residues_.size() != 0 );
	
	core::pack::task::operation::PreventRepacking pp;
	for( core::Size i = 1; i<=residues_.size(); ++i )
		pp.include_residue( residues_[i] );
	
	pp.apply( pose, task );
}

	
utility::vector1< core::Size >
PreventResiduesFromRepackingOperation::get_residues() const
{
	return residues_;
}

void
PreventResiduesFromRepackingOperation::set_residues( utility::vector1  < core::Size > residues_vec )
{
        runtime_assert( residues_vec.size() != 0 );
        residues_.clear();
        BOOST_FOREACH( core::Size const item, residues_vec )
          residues_.push_back( item );
}

	
void
PreventResiduesFromRepackingOperation::parse_tag( TagCOP tag , DataMap & )
{
        unparsed_residues_ = tag->getOption< std::string >( "residues" ) ;
        if( unparsed_residues_ != "" ){

          utility::vector1< std::string > const res_keys( utility::string_split( unparsed_residues_ , ',' ) );
          utility::vector1< core::Size > residues_vector;
          residues_.clear();

          BOOST_FOREACH( std::string const key, res_keys ){
            Size const res( utility::string2int( key ) );
            TR.Debug<<"not designing residue: "<< key  <<std::endl;
            residues_.push_back( res );
          }
        }
      }
		} //namespace protocols
	} //namespace toolbox
} //namespace task_operations
