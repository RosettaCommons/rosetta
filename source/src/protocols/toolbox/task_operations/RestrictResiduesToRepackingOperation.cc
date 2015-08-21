// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictResiduesToRepackingOperation.cc
/// @brief
/// @author Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictResiduesToRepackingOperation.hh>
#include <protocols/toolbox/task_operations/RestrictResiduesToRepackingOperationCreator.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>

// Project Headers

#include <core/pack/task/operation/TaskOperations.hh>


// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <boost/foreach.hpp>

// C++ Headers

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <boost/functional/hash.hpp>


using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictResiduesToRepackingOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;


RestrictResiduesToRepackingOperation::RestrictResiduesToRepackingOperation(): parent() {
	residues_.clear();
	reference_pdb_id_ = "";
}

RestrictResiduesToRepackingOperation::RestrictResiduesToRepackingOperation( utility::vector1 < core::Size > residues )
: parent(), residues_( residues ), reference_pdb_id_( "" )
{
}

RestrictResiduesToRepackingOperation::~RestrictResiduesToRepackingOperation() {}

core::pack::task::operation::TaskOperationOP
RestrictResiduesToRepackingOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictResiduesToRepackingOperation );
}

core::pack::task::operation::TaskOperationOP RestrictResiduesToRepackingOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictResiduesToRepackingOperation( *this ) );
}


void
RestrictResiduesToRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	if ( residues_.size() == 0 ) { // do nothing
		return;
	}

	core::pack::task::operation::RestrictResidueToRepacking rrtr;
	for ( core::Size i = 1; i<=residues_.size(); ++i ) {
		TR.Debug << "TASKOPERATION: restrict to repacking: " << residues_[i] << std::endl;
		rrtr.include_residue( residues_[i] );
	}
	rrtr.apply( pose, task );
}


utility::vector1< core::Size >
RestrictResiduesToRepackingOperation::get_residues() const
{
	return residues_;
}

void
RestrictResiduesToRepackingOperation::set_residues( utility::vector1  < core::Size > residues_vec )
{
	runtime_assert( residues_vec.size() != 0 );
	residues_.clear();
	BOOST_FOREACH ( core::Size const item, residues_vec ) {
		residues_.push_back( item );
	}
}

void
RestrictResiduesToRepackingOperation::parse_tag( TagCOP tag , DataMap & )
{
	reference_pdb_id_ = tag->getOption< std::string >( "reference_pdb_id", "" );
	if ( tag->getOption< std::string >( "residues" ) == "-1" ) {
		residues_.clear();
		TR<<"No residues specified to restrict repacking. I'm doing nothing"<<std::endl;
		return;
	}

	unparsed_residues_ = tag->getOption< std::string >( "residues" ) ;
	if ( unparsed_residues_ != "" ) {
		core::pose::Pose reference_pose;
		if ( reference_pdb_id_ != "" ) {
			core::import_pose::pose_from_pdb( reference_pose, reference_pdb_id_ );
		}

		utility::vector1< std::string > const res_keys( utility::string_split( unparsed_residues_ , ',' ) );
		utility::vector1< core::Size > residues_vector;
		residues_.clear();

		BOOST_FOREACH ( std::string const key, res_keys ) {
			Size res( utility::string2int( key ) );
			if ( reference_pdb_id_ != "" ) {
				res = core::pose::parse_resnum( key, reference_pose );
			}
			TR.Debug<<"not designing residue: "<< key<<" which is translated to residue "<<res  <<std::endl;
			residues_.push_back( res );
		}
	}
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
