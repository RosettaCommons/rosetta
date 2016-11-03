// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
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
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>

// C++ Headers
#include <set>

#include <utility/vector0.hh>

using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "devel.splice.DesignInterfacesOperation" );

namespace devel {
namespace splice {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

DesignInterfacesOperation::DesignInterfacesOperation() :
	design_shell_( 6.0 ),
	repack_shell_( 8.0 ),
	restrict_to_repacking_chain1_( false ),
	restrict_to_repacking_chain2_( true )
{
}

DesignInterfacesOperation::~DesignInterfacesOperation() = default;

core::pack::task::operation::TaskOperationOP
DesignInterfacesOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DesignInterfacesOperation );
}

void DesignInterfacesOperationCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	DesignInterfacesOperation::provide_xml_schema( xsd );
}

std::string DesignInterfacesOperationCreator::keyname() const {
	return DesignInterfacesOperation::keyname();
}

void DesignInterfacesOperation::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	task_op_schema_empty( xsd, keyname() );
}

core::pack::task::operation::TaskOperationOP DesignInterfacesOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new DesignInterfacesOperation( *this ) );
}
/// @brief restricts to repacking all residues outside of design_shell_ around each residue
void
DesignInterfacesOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task;
	using namespace protocols::toolbox::task_operations;
	using namespace protocols::rosetta_scripts;

	TaskFactoryOP tf( new TaskFactory );
	ProteinInterfaceDesignOperationOP pido( new ProteinInterfaceDesignOperation );
	pido->interface_distance_cutoff( design_shell() );
	pido->design_chain1(true);
	pido->design_chain2(true);
	tf->push_back( pido );

	///create a subpose containing only the two segments as two separate chains
	protocols::simple_moves::SwitchChainOrderMover sco;
	sco.chain_order( "1" );
	core::pose::Pose chain1( pose );
	sco.apply( chain1 ); // remove chain2 from pose
	protocols::protein_interface_design::movers::SetAtomTree sat;
	sat.two_parts_chain1( true );
	sat.apply( chain1 ); // create a fold tree around the cut in chain1
	core::pose::Pose segment1, segment2, combined;
	core::pose::partition_pose_by_jump( chain1, 1/*jump*/, segment1, segment2 ); // partition the pose into two separate chains around teh cut
	combined = segment1;
	core::pose::append_pose_to_pose(combined, segment2); // concatenate the segments as two separate chains
	combined.conformation().detect_disulfides(); // otherwise things go awry downstream

	// figure out which residues are at the two interfaces
	utility::vector1< core::Size > const two_segment_interface( residue_packer_states( combined, tf, true/*designable*/, true/*packable*/ ) );
	utility::vector1< core::Size > const two_chain_interface( residue_packer_states( pose, tf, true, true ) );

	// now design around each of those residues
	DesignAroundOperation dao;
	dao.repack_shell(repack_shell());
	dao.design_shell(design_shell());
	TR<<"residues found in two_segment interface: ";
	for ( core::Size const res : two_segment_interface ) {
		dao.include_residue( res );
		TR<<res<<'+';
	}
	TR<<std::endl<<"residues found in two_chain_interface: ";
	for ( core::Size const res : two_chain_interface ) {
		dao.include_residue(res);
		TR<<res<<'+';
	}
	TR<<std::endl<<"total residues: "<<two_segment_interface.size() + two_chain_interface.size()<<std::endl;
	dao.apply( pose, task );
	RestrictChainToRepackingOperation rctr;
	if ( restrict_to_repacking_chain1() ) {
		TR<<"Restricting chain1 to repacking"<<std::endl;
		rctr.chain( 1 );
		rctr.apply( pose, task );
	}
	if ( restrict_to_repacking_chain2() && pose.conformation().num_chains()>1 ) {
		TR<<"Restricting chain2 to repacking"<<std::endl;
		rctr.chain( 2 );
		rctr.apply( pose, task );
	}
}

void
DesignInterfacesOperation::design_shell( core::Real const radius )
{
	design_shell_ = radius;
	if ( radius >= repack_shell() ) {
		repack_shell_ = radius;
	}
}

void
DesignInterfacesOperation::parse_tag( TagCOP tag , DataMap & )
{
	design_shell( tag->getOption< core::Real >( "design_shell", 6.0 ) );
	repack_shell( tag->getOption<core::Real >("repack_shell", 8.0) );
	restrict_to_repacking_chain1( tag->getOption< bool >( "restrict_to_repacking_chain1", false ) );
	restrict_to_repacking_chain2( tag->getOption< bool >( "restrict_to_repacking_chain2", true  ) );
	TR<<"design shell: "<<design_shell()<<" repack shell: "<<repack_shell()<<" restrict_to_repacking_chain1: "<<restrict_to_repacking_chain1()<<" restrict_to_repacking_chain2: "<<restrict_to_repacking_chain2()<<std::endl;
}
} //namespace splice
} //namespace devel
