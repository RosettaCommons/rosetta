// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/DesignInterfacesOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <devel/splice/DesignInterfacesOperation.hh>
#include <devel/splice/DesignInterfacesOperationCreator.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

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
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/simple_moves/SwitchChainOrderMover.hh>
#include <protocols/simple_moves/CutChainMover.hh>
#include <core/pose/util.hh>
#include <protocols/protein_interface_design/movers/SetAtomTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <set>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "devel.splice.DesignInterfacesOperation" );

namespace devel {
namespace splice {

using namespace core::pack::task::operation;
using namespace std;

DesignInterfacesOperation::DesignInterfacesOperation() :
	design_shell_( 6.0 ),
	repack_shell_( 8.0 )
{
}

DesignInterfacesOperation::~DesignInterfacesOperation() {}

core::pack::task::operation::TaskOperationOP
DesignInterfacesOperationCreator::create_task_operation() const
{
	return new DesignInterfacesOperation;
}

core::pack::task::operation::TaskOperationOP DesignInterfacesOperation::clone() const
{
	return new DesignInterfacesOperation( *this );
}
///@brief restricts to repacking all residues outside of design_shell_ around each residue
void
DesignInterfacesOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task;
    using namespace protocols::toolbox::task_operations;
    using namespace protocols::rosetta_scripts;

    TaskFactoryOP tf = new TaskFactory;
    ProteinInterfaceDesignOperationOP pido = new ProteinInterfaceDesignOperation;
    pido->interface_distance_cutoff( design_shell() );
    pido->design_chain1(true);
    pido->design_chain2(true);
    tf->push_back( pido );

///create a subpose containing only the two segments as two separate chains
    protocols::simple_moves::SwitchChainOrderMover sco;
    sco.chain_order( "1" );
    core::pose::Pose chain1( pose );
    sco.apply( chain1 );
    protocols::protein_interface_design::movers::SetAtomTree sat;
    sat.two_parts_chain1( true );
    sat.apply( chain1 );
    core::pose::Pose segment1, segment2, combined;
    core::pose::partition_pose_by_jump( chain1, 1/*jump*/, segment1, segment2 );
    combined = segment1;
    core::pose::append_pose_to_pose(combined, segment2);
		combined.conformation().detect_disulfides();
//		core::pose::Pose part1, part2;
//		part1.clear(); part2.clear();
//		for (core::Size resj = 1; resj <= cut; ++resj) {
//			core::conformation::Residue const & rsd( copy_pose.residue( resj) );
//			part1.append_residue_by_bond( rsd );
//	}



//		core::pose::Pose part1( core::pose::create_subpose( chain1, 1,
//		core::pose::append_subpose_to_pose( chain

    utility::vector1< core::Size > const two_segment_interface( residue_packer_states( combined, tf, true/*designable*/, true/*packable*/ ) );
    utility::vector1< core::Size > const two_chain_interface( residue_packer_states( pose, tf, true, true ) );
    DesignAroundOperation dao;
    dao.repack_shell(repack_shell());
    dao.design_shell(design_shell());
    TR<<"residues found in two_segment interface: ";
    foreach( core::Size const res, two_segment_interface ){
        dao.include_residue( res );
        TR<<res<<'+';
    }
    TR<<std::endl<<"residues found in two_chain_interface: ";
    foreach( core::Size const res, two_chain_interface ){
        dao.include_residue(res);
        TR<<res<<'+';
    }
    TR<<std::endl<<"total residues: "<<two_segment_interface.size() + two_chain_interface.size()<<std::endl;
    dao.apply( pose, task );
}

void
DesignInterfacesOperation::design_shell( core::Real const radius )
{
	design_shell_ = radius;
	if( radius >= repack_shell() )
		repack_shell_ = radius;
}

void
DesignInterfacesOperation::parse_tag( TagPtr tag )
{
    tag->getOption< core::Real >( "design_shell", 6.0 );
    tag->getOption<core::Real >("repack_shell", 8.0);
		TR<<"design shell: "<<design_shell()<<" repack shell: "<<repack_shell()<<std::endl;
}
} //namespace splice
} //namespace devel
