// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/CDRsMinPackMin.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#include <protocols/antibody2/CDRsMinPackMin.hh>

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
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/antibody2/AntibodyInfo.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>




using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.CDRsMinPackMin");
using namespace core;

namespace protocols {
namespace antibody2 {

// default constructor
CDRsMinPackMin::CDRsMinPackMin() :
	Mover(),
	sc_min_( false ),
	rt_min_( false ),
	cen_cst_( 0.0 ),
	high_cst_( 0.0 ),
	benchmark_( false ),
	min_type_( "" ),
	Temperature_( 0.0 ),
	min_tolerance_( 0.0 )
{
	user_defined_ = false;
}


CDRsMinPackMin::CDRsMinPackMin(AntibodyInfoOP antibody_info) : Mover()
{
	user_defined_ = false;

	ab_info_ = antibody_info;

	init();
}

CDRsMinPackMin::CDRsMinPackMin(
	AntibodyInfoOP            antibody_info,
	core::scoring::ScoreFunctionOP  scorefxn,
	pack::task::TaskFactoryOP  tf,
	kinematics::MoveMapOP      movemap
) : Mover()
{
	user_defined_ = true;

	ab_info_               = antibody_info;
	loop_scorefxn_highres_ = scorefxn;
	tf_                    = tf;
	allcdr_map_            = movemap;

	init();
}



CDRsMinPackMin::~CDRsMinPackMin() {}



void CDRsMinPackMin::init(){

	benchmark_ = false;
	sc_min_ = false;
	rt_min_ = false;

	min_type_ = "dfpmin_armijo_nonmonotone";
	Temperature_ = 0.8;
	min_tolerance_ = 0.1;


	if (!user_defined_){
		tf_ = NULL;
		allcdr_map_ = NULL;

		// setup all the scoring functions
		loop_scorefxn_highres_ = core::scoring::ScoreFunctionFactory::create_score_function("standard", "score12" );
		loop_scorefxn_highres_->set_weight( core::scoring::chainbreak, 10. / 3. );
		loop_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10. / 3. );
	}

	cdr_sequence_move_ = new moves::SequenceMover();


}

void CDRsMinPackMin::finalize_setup( pose::Pose & pose )
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;
	using namespace protocols;
	using namespace protocols::toolbox::task_operations;
	using namespace protocols::moves;

	// **************** FoldTree ****************
	ab_info_->all_cdr_fold_tree( pose );
	TR<<pose.fold_tree()<<std::endl;

	// adding cutpoint variants for chainbreak score computation
	loops::remove_cutpoint_variants( pose, true ); //remove first
	loops::add_cutpoint_variants( pose );

	// must score first
	( *loop_scorefxn_highres_ )( pose );

	//**************** MoveMap ****************
	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
	utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );

	select_loop_residues( pose, *(ab_info_->get_all_cdr_loops()), false /*include_neighbors*/, bb_is_flexible );
	select_loop_residues( pose, *(ab_info_->get_all_cdr_loops()), true /*include_neighbors*/, sc_is_flexible );

	//for (Size kk=1;kk<=sc_is_flexible.size();kk++){
	//    TR<<kk<<"    "<<sc_is_flexible[kk]<<std::endl;
	//}

	if(!allcdr_map_){
		allcdr_map_ = new kinematics::MoveMap();
		allcdr_map_->clear();
		allcdr_map_->set_chi( false );
		allcdr_map_->set_bb( false );


		allcdr_map_->set_bb( bb_is_flexible );
		allcdr_map_->set_chi( sc_is_flexible );


		for( Size ii = 1; ii <= ab_info_->get_all_cdr_loops()->num_loop(); ii++ ){
			allcdr_map_->set_jump( ii, false );
		}//TODO: start from 1 or 2? should have a set function to handle this!
	}


	//**************** TaskFactory ****************
	if(!tf_){
	  tf_ = setup_packer_task(pose);
	  tf_->push_back( new RestrictToInterface( sc_is_flexible ) );//TODO: check this, no rb_jump here
	  //core::pack::task::PackerTaskOP my_task2(tf_->create_task_and_apply_taskoperations(pose));
	  //TR<<*my_task2<<std::endl; //exit(-1);
	}

	// 1. rotamer_trial
	simple_moves::RotamerTrialsMoverOP rotamer_trial_mover = new simple_moves::RotamerTrialsMover( loop_scorefxn_highres_, tf_ );
	cdr_sequence_move_->add_mover(rotamer_trial_mover);

	// 2. all_cdr_min_moves
	simple_moves::MinMoverOP  all_cdr_min_moves = new simple_moves::MinMover( allcdr_map_,loop_scorefxn_highres_, min_type_, min_tolerance_, true );
	cdr_sequence_move_ -> add_mover(all_cdr_min_moves);



	moves::MonteCarloOP mc = new moves::MonteCarlo(pose, *loop_scorefxn_highres_, Temperature_ );


	// 3. PackRotamer and Trial
	simple_moves::PackRotamersMoverOP repack = new simple_moves::PackRotamersMover( loop_scorefxn_highres_ );
	  repack->task_factory( tf_ );
	moves::TrialMoverOP repack_trial = new moves::TrialMover(repack, mc);
	cdr_sequence_move_ -> add_mover(repack_trial);


	// 4. optional, rt_min_ or sc_min_
	if ( rt_min_ ){
	  simple_moves::RotamerTrialsMinMoverOP rtmin = new simple_moves::RotamerTrialsMinMover( loop_scorefxn_highres_, tf_ );
	  moves::TrialMoverOP rtmin_trial = new moves::TrialMover( rtmin, mc );
	  cdr_sequence_move_ -> add_mover(rtmin_trial);
	}
	if ( sc_min_ ){
	  core::pack::task::TaskFactoryCOP my_tf( tf_); // input must be COP, weird
	  docking::SidechainMinMoverOP scmin_mover = new docking::SidechainMinMover( loop_scorefxn_highres_, my_tf );
	  moves::TrialMoverOP scmin_trial = new moves::TrialMover( scmin_mover, mc );
	  cdr_sequence_move_ -> add_mover(scmin_trial);
	}



}//finalize_setup




void CDRsMinPackMin::apply( pose::Pose & pose ) {
	finalize_setup(pose);
	cdr_sequence_move_ -> apply(pose);
}







/// @details  Show the complete setup of the antibody modeler protocol
void CDRsMinPackMin::show( std::ostream & out ) {
    out << *this;
}

std::ostream & operator<<(std::ostream& out, const CDRsMinPackMin & ab_m_2 ){
    using namespace ObjexxFCL::fmt;

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

void CDRsMinPackMin::set_task_factory(core::pack::task::TaskFactoryCOP tf){
    tf_ = new pack::task::TaskFactory(*tf);
}

void CDRsMinPackMin::set_move_map(kinematics::MoveMapCOP movemap){
    allcdr_map_ = new kinematics::MoveMap(*movemap);
}


} // end antibody2
} // end protocols

