// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingPrepackProtocol.cc
/// @brief Prepacking of the bound sturcture before docking.
/// @author Robin A Thottungal (raugust1@jhu.edu)
///   added to: JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <protocols/docking/DockingPrepackProtocol.hh>

// Package headers
#include <protocols/docking/util.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>


#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using namespace protocols::docking;
using namespace protocols::moves;
using namespace core;
using namespace pack::task;
using protocols::jd2::JobDistributor;
using namespace protocols::membrane;

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.DockingPrepackProtocol" );

namespace protocols {
namespace docking {

DockingPrepackProtocol::DockingPrepackProtocol(): DockingHighRes()
{
	Mover::type( "DockingPrepackProtocol" );
	setup_defaults();
	register_options();
	init_from_options();
}


void DockingPrepackProtocol::setup_defaults()
{
	trans_magnitude_ = 1000.0;
	pack_operations_ = SequenceMoverOP( new SequenceMover() );
	dock_ppk_ = false;
	membrane_ = false;
    movers_setup_ = false;
}

DockingPrepackProtocol::~DockingPrepackProtocol(){
}

void DockingPrepackProtocol::init_from_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::docking::dock_rtmin ].user() ) {
		set_rt_min(option[ OptionKeys::docking::dock_rtmin ]());
	}

	if ( option[ OptionKeys::docking::sc_min ].user() ) {
		set_sc_min(option[ OptionKeys::docking::sc_min ]());
	}

	if ( option[ OptionKeys::docking::partners ].user() ) {
		set_partners(option[ OptionKeys::docking::partners ]());
	}

	if ( option[ OptionKeys::docking::dock_ppk ].user() ) {
		set_dock_ppk(option[ OptionKeys::docking::dock_ppk ]());
	}

	if ( option[ OptionKeys::mp::setup::spanfiles ].user() ) {
		membrane_ = true;
	}
}

void DockingPrepackProtocol::set_dock_ppk(bool dock_ppk)
{
	dock_ppk_ = dock_ppk;
}

void DockingPrepackProtocol::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::dock_rtmin );
	option.add_relevant( OptionKeys::docking::sc_min );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::mp::setup::spanfiles );
}

void DockingPrepackProtocol::score_and_output(std::string outfilename,
	core::pose::Pose & pose )
{
	using namespace core::scoring;

	// Creating a job compatible with JD2
	static protocols::jd2::JobOP job = jd2::JobDistributor::get_instance()->current_job();
	//job_ = jd2::JobDistributor::get_instance()->current_job();

	core::Real score_pose  = ( *scorefxn() )( pose ); // scoring the pose

	// Getting the name of current job
	//std::string job_name (JobDistributor::get_instance()->job_outputter()->output_name( job ) );
	job->add_string_real_pair("E"+outfilename, score_pose);
	JobDistributor::get_instance()->job_outputter()->other_pose( job, pose, outfilename+"_");

}

void DockingPrepackProtocol::setup_pack_operation_movers()
{
	prepack_full_repack_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover() );
	prepack_full_repack_->score_function( scorefxn_pack() );
	prepack_full_repack_->task_factory( task_factory() );
	pack_operations_->add_mover(prepack_full_repack_);

	if ( rt_min() ) {
		rtmin_mover_ = protocols::simple_moves::RotamerTrialsMinMoverOP( new protocols::simple_moves::RotamerTrialsMinMover( ) );
		rtmin_mover_->score_function( scorefxn_pack() );
		rtmin_mover_->task_factory( task_factory() );
		pack_operations_->add_mover(rtmin_mover_);
	}
	if ( sc_min() ) {
		scmin_mover_ = SidechainMinMoverOP( new SidechainMinMover() );
		scmin_mover_->set_scorefxn( scorefxn_pack() );
		scmin_mover_->set_task_factory( task_factory() );
		pack_operations_->add_mover( scmin_mover_ );
	}
    movers_setup_ = true;
}

void DockingPrepackProtocol::finalize_setup( pose::Pose & pose ) {

	// create a membrane protein from the pose
	if ( membrane_ ) {
		membrane::AddMembraneMoverOP add_mem( new membrane::AddMembraneMover() );
		add_mem->apply( pose );

		// set the foldtree for membrane proteins
		core::Size dock_jump = create_membrane_docking_foldtree_from_partners( pose, partners() );

		// set DockJumps in DockingHighres protocol
		DockJumps dock_jumps;
		dock_jumps.push_back( static_cast< int >( dock_jump ) );
		set_movable_jumps( dock_jumps );

		// get correct scorefunction for membrane proteins and set it in the parent class
		core::scoring::ScoreFunctionOP mem_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
		set_scorefxn( mem_sfxn );
	} else {

		// Robin Notes!!!!!!!!!!!!!!!!!!!!!!!!!!
		// setup_foldtree asserts movable jumps to be atleast 1
		// I am not sure how the partner flag & movable_jumps communicate
		// What if the the partner flag is A_B_C and movable jumps is 1
		/// May need to write a method to change the number of movable
		// jumps based on the partner flag!!!!
		// Possible breaking of code, needs to get changed later
		docking::setup_foldtree( pose, partners(), movable_jumps() );
	}
	tf2()->set_prepack_only(true);
	tf2()->create_and_attach_task_factory( this, pose );
    // JRJ Says this should only be done once... if we pass a list this is done more than once... bad.
    if ( !movers_setup_ ) {
	    setup_pack_operation_movers();
    }
}

void DockingPrepackProtocol::apply( core::pose::Pose & pose )
{
	finalize_setup(pose);
	score_and_output("initial",pose);

	//Move each partners away from the others
	for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {

		// if membrane protein: translate in membrane plane
		if ( membrane_ ) {

			// get membrane axis
			core::Vector trans_axis( membrane_axis( pose, *jump ) );

			// create new translation mover
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(trans_axis, *jump) );

			// do actual translation
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		} else {
			// if not membrane protein
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		}
	}
	score_and_output("away",pose);

	// packing the unbound structures
	pack_operations_->apply( pose );
	score_and_output("away_packed",pose);

	//bringing the packed structures together
	for ( DockJumps::const_iterator jump= movable_jumps().begin(); jump != movable_jumps().end(); ++jump ) {

		// for membrane protein, translate in membrane plane
		if ( membrane_ ) {

			core::Vector trans_axis( protocols::membrane::membrane_axis( pose, *jump ) );
			rigid::RigidBodyTransMoverOP translate_back( new rigid::RigidBodyTransMover(trans_axis, *jump) );
			translate_back->step_size( trans_magnitude_ );

			// why this needs to be commented out to work is a mystery
			//   translate_back->trans_axis().negate();
			translate_back->apply(pose);
		} else {
			// if not membrane protein
			rigid::RigidBodyTransMoverOP translate_back ( new rigid::RigidBodyTransMover(pose, *jump) );

			translate_back->step_size( trans_magnitude_ );
			translate_back->trans_axis().negate();
			translate_back->apply(pose);

			//fa_dock_slide_into_contact_ = new FaDockingSlideIntoContact(*jump);
			//fa_dock_slide_into_contact_-> apply(pose);
		}
	}

	if ( dock_ppk_ ) {
		pack_operations_->apply( pose );
	}
	score_and_output("prepack",pose);
}

std::string DockingPrepackProtocol::get_name() const {
	return "DockingPrepackProtocol";
}

}
}
