// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/CDRsMinPackMin.cc
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#include <protocols/antibody/CDRsMinPackMin.hh>

#include <core/chemical/VariantType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.antibody.CDRsMinPackMin" );
using namespace core;

namespace protocols {
namespace antibody {


CDRsMinPackMin::CDRsMinPackMin(AntibodyInfoOP antibody_info) : Mover() {
	user_defined_ = false;

	ab_info_ = antibody_info;

	init();
}

CDRsMinPackMin::CDRsMinPackMin(
	AntibodyInfoOP            antibody_info,
	core::scoring::ScoreFunctionOP  scorefxn,
	pack::task::TaskFactoryOP  tf,
	kinematics::MoveMapOP      movemap
) : Mover() {
	user_defined_ = true;

	ab_info_               = antibody_info;
	loop_scorefxn_highres_ = scorefxn;
	tf_                    = tf;
	allcdr_map_            = movemap;

	init();
}


CDRsMinPackMin::~CDRsMinPackMin() = default;


void CDRsMinPackMin::init() {

	benchmark_ = false;
	sc_min_ = false;
	rt_min_ = false;
	turnoff_minimization_ = false;

	min_type_ = "lbfgs_armijo_nonmonotone";
	Temperature_ = 0.8;
	min_tolerance_ = 0.1;
	update_rounds_ = 0;

	if ( !user_defined_ ) {
		tf_ = nullptr;
		allcdr_map_ = nullptr;

		// setup all the scoring functions
		loop_scorefxn_highres_ = core::scoring::get_score_function();
		loop_scorefxn_highres_->set_weight( core::scoring::chainbreak, 10. / 3. );
		loop_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10. / 3. );
	}


}

void CDRsMinPackMin::finalize_setup( pose::Pose & pose ) {
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;
	using namespace protocols;
	using namespace protocols::toolbox::task_operations;
	using namespace protocols::moves;

	cdr_sequence_move_ = protocols::moves::SequenceMoverOP( new moves::SequenceMover() );

	// **************** FoldTree ****************
	pose.fold_tree( * ab_info_->get_FoldTree_AllCDRs(pose)  );
	TR<<pose.fold_tree()<<std::endl;

	// adding cutpoint variants for chainbreak score computation
	loops::remove_cutpoint_variants( pose, true ); //remove first
	loops::add_cutpoint_variants( pose );

	// must score first
	( *loop_scorefxn_highres_ )( pose );

	//**************** MoveMap ****************
	if ( !allcdr_map_ ) { // use this if, because sometimes a user may input a movemap at the beginning
		allcdr_map_ = core::kinematics::MoveMapOP( new kinematics::MoveMap() );
		*allcdr_map_=ab_info_->get_MoveMap_for_Loops(pose, *ab_info_->get_AllCDRs_in_loopsop(), false, true, 10.0);
	} else {
		if ( update_rounds_ > 0 ) {
			allcdr_map_->clear();
			*allcdr_map_=ab_info_->get_MoveMap_for_Loops(pose, *ab_info_->get_AllCDRs_in_loopsop(), false, true, 10.0);
		}
	}


	//**************** TaskFactory ****************
	if ( !tf_ ) { //use this if, because sometimes a user may input a taskfactory at the beginning
		tf_ = ab_info_->get_TaskFactory_AllCDRs(pose);

		//core::pack::task::PackerTaskOP my_task2(tf_->create_task_and_apply_taskoperations(pose));
		//TR<<*my_task2<<std::endl; //exit(-1);
	} else {
		if ( update_rounds_ > 0 ) {
			tf_->clear();
			tf_ = ab_info_->get_TaskFactory_AllCDRs(pose);
		}
	}

	// 1. rotamer_trial
	minimization_packing::RotamerTrialsMoverOP rotamer_trial_mover( new minimization_packing::RotamerTrialsMover( loop_scorefxn_highres_, tf_ ) );
	cdr_sequence_move_->add_mover(rotamer_trial_mover);

	// 2. all_cdr_min_moves
	minimization_packing::MinMoverOP  all_cdr_min_moves( new minimization_packing::MinMover( allcdr_map_,loop_scorefxn_highres_, min_type_, min_tolerance_, true ) );
	if ( !turnoff_minimization_ ) cdr_sequence_move_ -> add_mover(all_cdr_min_moves);


	moves::MonteCarloOP mc( new moves::MonteCarlo(pose, *loop_scorefxn_highres_, Temperature_ ) );


	// 3. PackRotamer and Trial
	minimization_packing::PackRotamersMoverOP repack( new minimization_packing::PackRotamersMover( loop_scorefxn_highres_ ) );
	repack->task_factory( tf_ );
	moves::TrialMoverOP repack_trial( new moves::TrialMover(repack, mc) );
	cdr_sequence_move_ -> add_mover(repack_trial);


	// 4. optional, rt_min_ or sc_min_
	if ( rt_min_ ) {
		minimization_packing::RotamerTrialsMinMoverOP rtmin( new minimization_packing::RotamerTrialsMinMover( loop_scorefxn_highres_, tf_ ) );
		moves::TrialMoverOP rtmin_trial( new moves::TrialMover( rtmin, mc ) );
		cdr_sequence_move_ -> add_mover(rtmin_trial);
	}
	if ( sc_min_ ) {
		core::pack::task::TaskFactoryCOP my_tf( tf_); // input must be COP, weird
		docking::SidechainMinMoverOP scmin_mover( new docking::SidechainMinMover( loop_scorefxn_highres_, my_tf ) );
		moves::TrialMoverOP scmin_trial( new moves::TrialMover( scmin_mover, mc ) );
		cdr_sequence_move_ -> add_mover(scmin_trial);
	}


}//finalize_setup


void CDRsMinPackMin::apply( pose::Pose & pose ) {
	finalize_setup(pose);
	cdr_sequence_move_ -> apply(pose);
	update_rounds_++;

}


/// @details  Show the complete setup of the antibody modeler protocol
void CDRsMinPackMin::show( std::ostream & out ) const {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const CDRsMinPackMin & ab_m_2 ) {
	using namespace ObjexxFCL::format;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, "Rosetta 3 Antibody Modeler" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;
	out << line_marker << "  sc_min                : " << ab_m_2.sc_min_     << std::endl;
	out << line_marker << "  rt_min                : " << ab_m_2.rt_min_    << std::endl;
	out << line_marker << std::endl;

	// Display the state of the antibody modeler protocol that will be used
	out << line_marker << std::endl;
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}


std::string CDRsMinPackMin::get_name() const {
	return "CDRsMinPackMin";
}

void CDRsMinPackMin::set_task_factory(core::pack::task::TaskFactoryCOP tf) {
	tf_ = core::pack::task::TaskFactoryOP( new pack::task::TaskFactory(*tf) );
}

void CDRsMinPackMin::set_move_map(kinematics::MoveMapCOP movemap) {
	allcdr_map_ = core::kinematics::MoveMapOP( new kinematics::MoveMap(*movemap) );
}


} // end antibody
} // end protocols

