// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/splice/RestrictToAlignedSegmentsOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToAlignedSegments.hh>
#include <protocols/toolbox/task_operations/RestrictToAlignedSegmentsCreator.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pack/task/TaskFactory.hh>
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
static basic::Tracer TR( "protocols.toolbox.task_operations.RestrictToAlignedSegmentsOperation" );

namespace protocols {
namespace toolbox{
namespace task_operations {

using namespace core::pack::task::operation;
using namespace std;

RestrictToAlignedSegmentsOperation::RestrictToAlignedSegmentsOperation()
{
	source_pose_.clear();
	start_res_.clear();
	stop_res_.clear();
	chain_ = 1;
	repack_shell_ = 6.0;
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
	//using namespace protocols::rosetta_scripts; //cannot be used because of TR in this namespace

	std::set< core::Size > designable;
	designable.clear();
	for( core::Size count = 1; count <= source_pose_.size(); ++count ){
		core::Size const nearest_to_from = protocols::rosetta_scripts::find_nearest_res( pose, *source_pose_[ count ], start_res_[ count ], chain() );
		core::Size const nearest_to_to = protocols::rosetta_scripts::find_nearest_res( pose, *source_pose_[ count ], std::min( stop_res_[ count ], source_pose_[ count ]->total_residue() ), chain() );

		TR<<"Finding nearest residue to residue pair "<<start_res_[ count ]<<','<<stop_res_[ count ]<<" in source pose"<<std::endl;
		if( nearest_to_from == 0 || nearest_to_to == 0 ){
			TR<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
			continue;
		}
		for( core::Size position = nearest_to_from; position <= nearest_to_to; ++position )
			designable.insert( position );
  }
/// in the following we use dao to compute the residues that surround the aligned region. We then go over these residues to make sure they're within the target chain
	protocols::toolbox::task_operations::DesignAroundOperationOP dao = new protocols::toolbox::task_operations::DesignAroundOperation;
	dao->design_shell( 0.1 );
	dao->repack_shell( repack_shell() );
	foreach( core::Size const d, designable ){
		dao->include_residue( d );
	}
	core::pack::task::TaskFactoryOP dao_tf = new core::pack::task::TaskFactory;
	dao_tf->push_back( dao );

	utility::vector1< core::Size > const surrounding_shell( protocols::rosetta_scripts::residue_packer_states( pose, dao_tf, false/*designable*/, true/*packable*/ ) );
	utility::vector1< core::Size > const designed_residues( protocols::rosetta_scripts::residue_packer_states( pose, dao_tf, true/*designable*/, false/*packable*/ ) );

	utility::vector1< core::Size > repackable, immutable;
	repackable.clear(); immutable.clear();
	for( core::Size i = pose.conformation().chain_begin( chain() ); i<=pose.conformation().chain_end( chain() ); ++i ){
		if( std::find( designed_residues.begin(), designed_residues.end(), i ) != designed_residues.end() ) // don't change designed residues
			continue;
		if( std::find( surrounding_shell.begin(), surrounding_shell.end(), i ) == surrounding_shell.end() )
			immutable.push_back( i );
		else
			repackable.push_back( i );
	}
///for some unfathomable reason OperateOnCertainResidues defaults to applying to all residues if none are defined, so you have to be careful here...
	OperateOnCertainResidues oocr_repackable, oocr_immutable;
	oocr_immutable.op( new PreventRepackingRLT );
	oocr_repackable.op( new RestrictToRepackingRLT );
	if( repackable.size() ){
		oocr_repackable.residue_indices( repackable );
		oocr_repackable.apply( pose, task );
		TR<<"allowing repacking in: ";
		foreach( core::Size const r, repackable )
			TR<<r<<' ';
		TR<<std::endl;
	}
	if( immutable.size() ){
		oocr_immutable.residue_indices( immutable );
		oocr_immutable.apply( pose, task );
		TR<<"no repack in: ";
		foreach( core::Size const i, immutable )
			TR<<i<<' ';
		TR<<std::endl;
	}
}

void
RestrictToAlignedSegmentsOperation::parse_tag( TagPtr tag )
{
	using namespace protocols::rosetta_scripts;
	utility::vector1< std::string > pdb_names, start_res, stop_res;
	pdb_names.clear(); start_res.clear(); stop_res.clear();
	if( tag->hasOption( "source_pdb" ) )
		pdb_names.push_back( tag->getOption< std::string >( "source_pdb" ) );
	if( tag->hasOption( "start_res" ) )
		start_res.push_back( tag->getOption< std::string >( "start_res" ) );
	if( tag->hasOption( "stop_res" ) )
		stop_res.push_back( tag->getOption< std::string >( "stop_res" ) );

	chain( tag->getOption< core::Size >( "chain", 1 ) );
	if( tag->hasOption( "source_pdb" ) || tag->hasOption( "start_res" ) || tag->hasOption( "stop_res" ) ){
		runtime_assert( tag->hasOption( "source_pdb" ) && tag->hasOption( "start_res" ) && tag->hasOption( "stop_res" ) );
	}

	utility::vector0< TagPtr > const btags( tag->getTags() );
	foreach( TagPtr const btag, btags ){
		pdb_names.push_back( btag->getOption< std::string >( "source_pdb" ) );
		start_res.push_back( btag->getOption< std::string >( "start_res" ) );
		stop_res.push_back( btag->getOption< std::string >( "stop_res" ) );
	}

	for( core::Size i = 1; i <= pdb_names.size(); ++i ){
		if( i == 1 || pdb_names[ i ] != pdb_names[ i - 1 ]){ // scrimp on reading from disk
			source_pose_.push_back( new core::pose::Pose );
			core::import_pose::pose_from_pdb( *source_pose_[ i ], pdb_names[ i ] );
		}
		else
			source_pose_.push_back( source_pose_[ i - 1 ] );
		core::Size const parsed_start( core::pose::parse_resnum( start_res[ i ], *source_pose_[ i ] ) );
		core::Size const parsed_stop ( core::pose::parse_resnum( stop_res[ i ], *source_pose_[ i ] ) );
		start_res_.push_back( parsed_start );
		stop_res_. push_back( parsed_stop );
	}
	repack_shell( tag->getOption< core::Real >( "repack_shell", 6.0 ));
}

} //namespace task_operations
} // namespace toolbox
} //namespace protocols
