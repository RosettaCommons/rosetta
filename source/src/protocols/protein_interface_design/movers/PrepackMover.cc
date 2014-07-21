// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PrepackMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PrepackMover.hh>
#include <protocols/protein_interface_design/movers/PrepackMoverCreator.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.PrepackMover" );

std::string
PrepackMoverCreator::keyname() const
{
	return PrepackMoverCreator::mover_name();
}

protocols::moves::MoverOP
PrepackMoverCreator::create_mover() const {
	return new PrepackMover;
}

std::string
PrepackMoverCreator::mover_name()
{
	return "Prepack";
}

PrepackMover::PrepackMover() :
	protocols::simple_moves::PackRotamersMover( PrepackMoverCreator::mover_name() ),
	scorefxn_( NULL ),
	jump_num_( 0 ),
	min_bb_( false ),
	mm_( NULL )
{}

PrepackMover::PrepackMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size jump_num
) :
	protocols::simple_moves::PackRotamersMover( PrepackMoverCreator::mover_name() ),
	scorefxn_( scorefxn ),
	jump_num_( jump_num )
{}


protocols::moves::MoverOP
PrepackMover::clone() const {
	return( protocols::moves::MoverOP( new PrepackMover( *this ) ) );
}

protocols::moves::MoverOP
PrepackMover::fresh_instance() const {
	return protocols::moves::MoverOP( new PrepackMover );
}

PrepackMover::~PrepackMover() {}

/// @details separate bound partners (if any), minimize, do rotamer trials, and re-minimize.
void PrepackMover::apply( pose::Pose & pose )
{
	// make a local packertask, reading resfiles and including current rotamers, excluding disulfides
	TR << "Performing repack..." << std::endl;
	using namespace core::pack::task;
	TaskFactoryOP tf;
	if( task_factory() ) tf = new TaskFactory( *task_factory() );
	else tf = new TaskFactory;
	tf->push_back( new operation::InitializeFromCommandline );
	tf->push_back( new operation::IncludeCurrent );
	tf->push_back( new operation::RestrictToRepacking );
	tf->push_back( new operation::NoRepackDisulfides );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new core::pack::rotamer_set::UnboundRotamersOperation();
	unboundrot->initialize_from_command_line();
	core::pack::task::operation::AppendRotamerSetOP unboundrot_operation = new core::pack::task::operation::AppendRotamerSet( unboundrot );
	tf->push_back( unboundrot_operation );
	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	using namespace protocols::toolbox::task_operations;
	if (basic::options::option[ basic::options::OptionKeys::docking::norepack1 ]()) tf->push_back( new DockingNoRepack1( jump_num_) );
	if (basic::options::option[ basic::options::OptionKeys::docking::norepack2 ]()) tf->push_back( new DockingNoRepack2( jump_num_) );

	//in case there is a resfile, information in this resfile overrides the computed task
	if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		tf->push_back( new operation::ReadResfile );
	}
	PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );

	TR << "Pre-minimizing structure..." << std::endl;
	core::kinematics::MoveMapOP mm_general;
	if( min_bb() && mm() ){
		mm_general = mm()->clone();
		for ( core::Size i = 1; i <= pose.total_residue(); ++i) {
			if ( !pose.residue(i).is_protein() ) {
				mm_general->set_chi( i, false );
				continue;
			}
			// Check for disulfide bonded cysteines
			if( pose.residue(i).type().name() == "CYD" ) mm_general->set_chi( i, false );
		}
  }
	else{
		mm_general = new core::kinematics::MoveMap;
		mm_general->clear();
	}
	if( min_bb() ){ //premin bb+sc
		if( !mm() ) mm_general->set_bb( true );
		protocols::simple_moves::MinMover min_bb_mover( mm_general, scorefxn_, "dfpmin_armijo_nonmonotone", 1e-5, true/*nblist*/, false/*deriv_check*/  );
		min_bb_mover.apply( pose );
	}

	// separate any bound partners
	protocols::rigid::RigidBodyTransMoverOP translate;
	if( (jump_num_ > 0) && (pose.conformation().num_chains() > 1) ) {
		TR<<"Translating along jump #"<<jump_num_<<std::endl;
		translate = new protocols::rigid::RigidBodyTransMover( pose, jump_num_ ) ;
		translate->step_size( 1000.0 );
		translate->apply( pose );
	}
  mm_general->set_bb( false );
	mm_general->set_jump( false );
	protocols::simple_moves::MinMover min_mover( mm_general, scorefxn_, "dfpmin_armijo_nonmonotone", 1e-5, true/*nblist*/, false/*deriv_check*/  );
	// pre-minimize sidechains
	min_mover.apply( pose );

	if( basic::options::option[basic::options::OptionKeys::docking::dock_rtmin].user() ) {
		protocols::simple_moves::RotamerTrialsMinMover rtmin( scorefxn_, tf );
		rtmin.apply( pose );
	}
	else {
		protocols::simple_moves::PackRotamersMover pack( scorefxn_, task );
		pack.apply( pose );
	}

	// post-minimize
	// using packer include_current() will make sure that these minimized rotamers are used.
	TR << "Post-minimizing structure..." << std::endl;
	min_mover.apply( pose );
	TR << "Done!\n";

	// move back together
	if( (jump_num_ > 0) && (pose.conformation().num_chains() > 1) ) {
		translate->trans_axis().negate();
		translate->apply( pose );
	}

	//final rescore to get everyone on the same page
	(*scorefxn_)(pose);
	TR.flush();
}

std::string
PrepackMover::get_name() const {
	return PrepackMoverCreator::mover_name();
}

void
PrepackMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	jump_num_ = tag->getOption<core::Size>("jump_number", 1 );
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	min_bb( tag->getOption< bool >( "min_bb", 0 ));
	if( min_bb() ) {
		mm_ = new core::kinematics::MoveMap;
		mm()->clear();
		protocols::rosetta_scripts::parse_movemap( tag, pose, mm_, data );
	}
	TR << "Prepack mover with scorefxn " << rosetta_scripts::get_score_function_name(tag) << " over jump number " << jump_num_ << "with min_bb "<<min_bb()<<std::endl;
}

void
PrepackMover::min_bb( bool const m ){
	min_bb_ = m;
}

bool
PrepackMover::min_bb() const{
	return min_bb_;
}

core::kinematics::MoveMapOP
PrepackMover::mm() const{
	if( !min_bb() ) TR<<"Warning: movemap requested but min_bb is set to false. This is probably wrong!"<<std::endl;
	return mm_;
}

void
PrepackMover::mm( core::kinematics::MoveMapOP mm ){
	mm_ = mm;
}

} //movers
} //protein_interface_design
} //protocols
