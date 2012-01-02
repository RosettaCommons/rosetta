// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToAlignedSegmentsOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToAlignedSegments.hh>
#include <protocols/toolbox/task_operations/RestrictToAlignedSegmentsCreator.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/import_pose/import_pose.hh>

// Project Headers
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <set>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictToAlignedSegmentsOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace std;

RestrictToAlignedSegmentsOperation::RestrictToAlignedSegmentsOperation()
{
	source_pdb_.clear();
	start_res_.clear();
	stop_res_.clear();
}

RestrictToAlignedSegmentsOperation::~RestrictToAlignedSegmentsOperation() {}

core::pack::task::operation::TaskOperationOP
RestrictToAlignedSegmentsOperationCreator::create_task_operation() const
{
	return new RestrictToAlignedSegmentsOperation;
}

core::pack::task::operation::TaskOperationOP RestrictToAlignedSegmentsOperation::clone() const
{
	return new RestrictToAlignedSegmentsOperation( *this );
}

void
RestrictToAlignedSegmentsOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace protocols::rosetta_scripts;

	std::set< core::Size > designable;
	designable.clear();
	core::pose::Pose source_pose;
	for( core::Size count = 1; count <= source_pdb_.size(); ++count ){
		if( count == 1 || source_pdb_[ count ] != source_pdb_[ count - 1 ]){ // scrimp on reading from disk
			core::import_pose::pose_from_pdb( source_pose, source_pdb_[ count ] );
		}
		core::Size const parsed_start( parse_resnum( start_res_[ count ], source_pose ) );
		core::Size const parsed_stop ( parse_resnum( stop_res_[ count  ], source_pose ) );
		core::Size const nearest_to_from = find_nearest_res( pose, source_pose, parsed_start );
		core::Size const nearest_to_to = find_nearest_res( pose, source_pose, parsed_stop );

		if( nearest_to_from == 0 || nearest_to_to == 0 ){
			TR<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
			continue;
		}
		for( core::Size position = nearest_to_from; position <= nearest_to_to; ++position )
			designable.insert( position );
  }

	utility::vector1< core::Size > restrict_to_repacking;
	restrict_to_repacking.clear();
	TR<<"Repackable residues: ";
	for( core::Size i = 1; i<=pose.total_residue(); ++i ){
		if( std::find( designable.begin(), designable.end(), i ) == designable.end() ){
			TR<<i<<",";
			restrict_to_repacking.push_back( i );
		}
	}
///for some unfathomable reason OperateOnCertainResidues defaults to applying to all residues if none are defined, so you have to be careful here...
	OperateOnCertainResidues oocr_repacking;
	if( restrict_to_repacking.size() ){
		oocr_repacking.op( new RestrictToRepackingRLT );
		oocr_repacking.residue_indices( restrict_to_repacking );
		oocr_repacking.apply( pose, task );
	}
}

void
RestrictToAlignedSegmentsOperation::parse_tag( TagPtr tag )
{
	using namespace protocols::rosetta_scripts;
	if( tag->hasOption( "source_pdb" ) )
		source_pdb_.push_back( tag->getOption< std::string >( "source_pdb" ) );
	if( tag->hasOption( "start_res" ) )
		start_res_.push_back( tag->getOption< std::string >( "start_res" ) );
	if( tag->hasOption( "stop_res" ) )
		stop_res_.push_back( tag->getOption< std::string >( "stop_res" ) );

	if( tag->hasOption( "source_pdb" ) || tag->hasOption( "start_res" ) || tag->hasOption( "stop_res" ) ){
		runtime_assert( tag->hasOption( "source_pdb" ) && tag->hasOption( "start_res" ) && tag->hasOption( "stop_res" ) );
		runtime_assert( start_res_[ 1 ] < stop_res_[ 1 ] );
	}

	utility::vector0< TagPtr > const btags( tag->getTags() );
	foreach( TagPtr const btag, btags ){
		source_pdb_.push_back( btag->getOption< std::string >( "source_pdb" ) );
		start_res_.push_back( btag->getOption< std::string >( "start_res" ) );
		stop_res_.push_back( btag->getOption< std::string >( "stop_res" ) );
	}
}
} //namespace protocols
} //namespace toolbox
} //namespace task_operations
