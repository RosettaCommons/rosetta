// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/SpinMover
/// @author Eva-Maria Strauch ( evas01@u.washinigton.edu )

//unit header
#include <protocols/protein_interface_design/movers/SpinMover.hh>
#include <protocols/protein_interface_design/movers/SpinMoverCreator.hh>
//project header
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <string>

//#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.SpinMover" );

std::string SpinMoverCreator::keyname() const
{
        return SpinMoverCreator::mover_name();
}

protocols::moves::MoverOP
SpinMoverCreator::create_mover() const {
        return protocols::moves::MoverOP( new SpinMover );
}

std::string
SpinMoverCreator::mover_name() {
        return "SpinMover";
}

SpinMover::SpinMover( ) :
	protocols::moves::Mover( SpinMoverCreator::mover_name()  )
	//protocols::moves::Mover ( "SpinMover" )
{ }

SpinMover::SpinMover( core::Size jump_num ) :
	protocols::moves::Mover( SpinMoverCreator::mover_name()  ),
	// protocols::moves::Mover ( "SpinMover" ),
   jump_num_(jump_num)
{ }

std::string SpinMover::get_name() const {
 	       return SpinMoverCreator::mover_name();
}

void
SpinMover::apply( core::pose::Pose & pose )
{
	core::kinematics::Jump const start_jump = pose.jump( jump_num_ );
	core::kinematics::Jump curr_jump( start_jump );
	TR<<"current jump: " <<curr_jump<<std::endl;

	TR<<"using fold-tree: "<< pose.fold_tree()<<std::endl;

	//first determine where the jump is
	core::Size const upstream_res(pose.fold_tree().jump_edge(jump_num_).start());
	core::Size const downstream_res(pose.fold_tree().jump_edge(jump_num_).stop());
	std::string upstream_atom(pose.fold_tree().jump_edge(jump_num_).upstream_atom());
	std::string downstream_atom(pose.fold_tree().jump_edge(jump_num_).downstream_atom());
	if( upstream_atom == "" ) upstream_atom = "C";
	if( downstream_atom == "" ) downstream_atom = "C";
	TR<<"upstream residue: "<<upstream_res<<" and atom "<<upstream_atom<<std::endl;

	core::kinematics::Stub downstream_stub = pose.conformation().upstream_jump_stub( jump_num_ );
	//calculate the rotation axis
	//looking down the axis from the upstream to downstream atom, positive rotations are counterclockwise
	core::Vector axis( pose.residue(upstream_res).atom(upstream_atom).xyz()//minus
										 - pose.residue(downstream_res).atom(downstream_atom).xyz() );


	numeric::xyzVector<double> reference_center = pose.residue(downstream_res).atom(downstream_atom).xyz();
	curr_jump.rotation_by_axis( downstream_stub, axis, reference_center, 360.0f*numeric::random::rg().uniform() /*degrees*/ );
	TR<<"new jump: " << curr_jump<<std::endl;
	 TR<<"new fold-tree: "<< pose.fold_tree()<<std::endl;
	pose.set_jump( jump_num_, curr_jump );
	  TR<<"new jump: " << curr_jump<<std::endl;
   TR<<"new fold-tree: "<< pose.fold_tree()<<std::endl;
}
//mjo commenting out 'data' and 'pose' because they are unused and cause warnings
void
SpinMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & /*data*/, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & /*pose*/ )
{
        jump_num_ = tag->getOption<core::Size>( "jump_num", 1);
        TR<<"SpinMover was instantiated "<<std::endl;
}


}//movers
}//protein_interface_design
}//protocols
