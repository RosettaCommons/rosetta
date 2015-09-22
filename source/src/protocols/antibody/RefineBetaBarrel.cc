// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/RefineBetaBarrel.cc
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/antibody/RefineBetaBarrel.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/LHRepulsiveRamp.hh>
#include <protocols/antibody_legacy/LHSnugFitLegacy.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/util.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.RefineBetaBarrel" );

using namespace core;
namespace protocols {
namespace antibody {


// default constructor
RefineBetaBarrel::RefineBetaBarrel() : Mover() {}
RefineBetaBarrel::~RefineBetaBarrel() {}


RefineBetaBarrel::RefineBetaBarrel(AntibodyInfoOP antibody_info) : Mover() {
	user_defined_ = false;
	ab_info_ = antibody_info;

	init();
}

RefineBetaBarrel::RefineBetaBarrel(AntibodyInfoOP antibody_info,
	core::scoring::ScoreFunctionCOP dock_scorefxn,
	core::scoring::ScoreFunctionCOP pack_scorefxn) : Mover() {
	user_defined_ = true;
	ab_info_ = antibody_info;
	dock_scorefxn_ = dock_scorefxn->clone();
	pack_scorefxn_ = pack_scorefxn->clone();

	init();
}


void RefineBetaBarrel::init( ) {

	if ( !user_defined_ ) {
		dock_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		dock_scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
		dock_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
		pack_scorefxn_ = core::scoring::get_score_function_legacy(  core::scoring::PRE_TALARIS_2013_STANDARD_WTS  );
	}

	repulsive_ramp_ = true;
	sc_min_ = false;
	rt_min_ = false;
	LH_dock_jump_.push_back(1);

}


void RefineBetaBarrel::finalize_setup(pose::Pose & pose ) {
	TR<<"   start finalize_setup function ..."<<std::endl;

	// add scores to map
	( *dock_scorefxn_ )( pose );

	// ************ MoveMap *************
	cdr_dock_map_ = kinematics::MoveMapOP( new kinematics::MoveMap() );
	*cdr_dock_map_=ab_info_->get_MoveMap_for_LoopsandDock(pose, *ab_info_->get_AllCDRs_in_loopsop(), false, true, 10.0);


	// ************ TaskFactory ************
	//set up sidechain movers for rigid body jump and loop & neighbors
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	// selecting movable c-terminal residues
	utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );
	select_loop_residues( pose, *(ab_info_->get_AllCDRs_in_loopsop()), true/*include_neighbors*/, sc_is_flexible);

	ObjexxFCL::FArray1D_bool loop_residues( pose.total_residue(), false );
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		loop_residues(i) = sc_is_flexible[i];
	} // check mapping

	using namespace protocols::toolbox::task_operations;
	if ( !tf_ ) {
		tf_= setup_packer_task(pose);
		tf_->push_back( TaskOperationCOP( new RestrictToInterface( LH_dock_jump_, loop_residues ) ) );
	}

	core::pack::task::PackerTaskOP my_task2(tf_->create_task_and_apply_taskoperations(pose));
	//TR<<*my_task2<<std::endl; //exit(-1);


	//************  FoldTree ************
	pose.fold_tree( * ab_info_->get_FoldTree_AllCDRs_LHDock(pose)   );
	TR<<pose.fold_tree()<<std::endl;


	//************  Variants ************
	// JQX:
	// 1. setting up the fold_tree doesn't automatically
	//    update the variants in residue_type
	// 2. the variants you saw from the PackTask are not the
	//    same as the variants in residue_type
	// 3. access the variants by a). pose.residue_type(i).variant_types()[1]
	//                           b). pose.residue(i).type().variant_types()[1]

	loops::remove_cutpoint_variants( pose, true ); //remove first
	loops::add_cutpoint_variants( pose ); // add back, based on the cutpoints defined by fold_tree

	/*
	for (Size i=1; i<=pose.total_residue();i++) {
	if (pose.residue(i).type().variant_types().size()>0){
	TR<<"residue "<<i<<"    "<< pose.residue(i).type().variant_types()[1]<<std::endl;
	}
	else{ TR<<"residue "<<i<<std::endl; }
	}
	exit(-1);
	*/
	TR<<"   finish finalize_setup function !!!"<<std::endl;

}


void RefineBetaBarrel::apply( pose::Pose & pose ) {

	finalize_setup(pose);

	//JQX:
	// the repulsive_ramp_ docking mover is very general now based on Jeff's request!
	// it will be moved to DockingProtocol soon
	// one must specify fold_tree and variants before using this mover
	if ( repulsive_ramp_ ) {
		lh_repulsive_ramp_ = LHRepulsiveRampOP( new LHRepulsiveRamp(LH_dock_jump_, dock_scorefxn_, pack_scorefxn_) );
		lh_repulsive_ramp_ -> set_move_map(cdr_dock_map_);
		lh_repulsive_ramp_ -> set_task_factory(tf_);
		if ( sc_min_ ) lh_repulsive_ramp_ -> set_sc_min(true);
		if ( rt_min_ ) lh_repulsive_ramp_ -> set_rt_min(true);
		lh_repulsive_ramp_->apply(pose);
		TR<<"   finish repulsive ramping !"<<std::endl;
	}


	dock_mcm_protocol_ = docking::DockMCMProtocolOP( new docking::DockMCMProtocol( LH_dock_jump_, dock_scorefxn_, pack_scorefxn_ ) );
	dock_mcm_protocol_ -> set_task_factory(tf_);
	dock_mcm_protocol_ -> set_move_map(cdr_dock_map_);
	if ( sc_min_ ) dock_mcm_protocol_ -> set_sc_min(true);
	if ( rt_min_ ) dock_mcm_protocol_ -> set_rt_min(true);
	dock_mcm_protocol_ -> apply(pose);

	TR<<"   finish L_H Docking !"<<std::endl;
	TR<<"FINISH BETA BARREL REFINEMENT STEP !! "<<std::endl;
}


std::string RefineBetaBarrel::get_name() const {
	return "RefineBetaBarrel";
}

void RefineBetaBarrel::set_task_factory(core::pack::task::TaskFactoryCOP tf) {
	tf_ = pack::task::TaskFactoryOP( new pack::task::TaskFactory(*tf) );
}

void RefineBetaBarrel::set_dock_score_func(core::scoring::ScoreFunctionCOP dock_scorefxn ) {
	dock_scorefxn_ = dock_scorefxn->clone();
}

void RefineBetaBarrel::set_pack_score_func(core::scoring::ScoreFunctionCOP pack_scorefxn) {
	pack_scorefxn_ = pack_scorefxn->clone();
}


}// namespace antibody
}// namespace protocols

