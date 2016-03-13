// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/LHRepulsiveRampLegacy.cc
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody_legacy/LHRepulsiveRampLegacy.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::format;

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/antibody/util.hh>


#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.LHRepulsiveRampLegacy" );


using namespace core;
namespace protocols {
namespace antibody {


// default constructor
LHRepulsiveRampLegacy::LHRepulsiveRampLegacy() : Mover() {

}


LHRepulsiveRampLegacy::LHRepulsiveRampLegacy(loops::Loops loops_in ) : Mover() {
	user_defined_ = false;

	init(loops_in, false);
}


LHRepulsiveRampLegacy::LHRepulsiveRampLegacy(AntibodyInfoOP antibody_in) : Mover() {
	user_defined_ = false;

	init( *(antibody_in->get_AllCDRs_in_loopsop()),false);
}

LHRepulsiveRampLegacy::LHRepulsiveRampLegacy(AntibodyInfoOP antibody_in, bool camelid) : Mover() {
	user_defined_ = false;

	init(*(antibody_in->get_AllCDRs_in_loopsop()), camelid);
}


LHRepulsiveRampLegacy::LHRepulsiveRampLegacy( loops::Loops loops_in,
	core::scoring::ScoreFunctionCOP dock_scorefxn,
	core::scoring::ScoreFunctionCOP pack_scorefxn ) : Mover() {
	user_defined_ = true;
	dock_scorefxn_ = dock_scorefxn->clone();
	pack_scorefxn_ = pack_scorefxn->clone();

	init(loops_in, false);
}

LHRepulsiveRampLegacy::LHRepulsiveRampLegacy( AntibodyInfoOP antibody_in,
	core::scoring::ScoreFunctionCOP dock_scorefxn,
	core::scoring::ScoreFunctionCOP pack_scorefxn ) : Mover() {
	user_defined_ = true;
	dock_scorefxn_ = dock_scorefxn->clone();
	pack_scorefxn_ = pack_scorefxn->clone();

	init(*(antibody_in->get_AllCDRs_in_loopsop()),false);
}


// default destructor
LHRepulsiveRampLegacy::~LHRepulsiveRampLegacy() {}

//clone
protocols::moves::MoverOP LHRepulsiveRampLegacy::clone() const {
	return( protocols::moves::MoverOP( new LHRepulsiveRampLegacy() ) );
}


void LHRepulsiveRampLegacy::init(loops::Loops loops_in, bool camelid ) {
	set_default();

	is_camelid_ = camelid;
	all_loops_ = loops_in;

	tf_ = pack::task::TaskFactoryOP( new pack::task::TaskFactory );

}


void LHRepulsiveRampLegacy::set_default() {
	benchmark_       = false;

	rep_ramp_cycles_ = 3 ;
	rot_mag_         = 2.0 ;
	trans_mag_       = 0.1 ;
	temperature_     = 0.8;
	min_threshold_   = 15.0;
	num_repeats_     = 4;
	min_type_        = "lbfgs_armijo_nonmonotone";

	if ( !user_defined_ ) {
		dock_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		dock_scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
		dock_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
		pack_scorefxn_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

}


std::string LHRepulsiveRampLegacy::get_name() const {
	return "LHRepulsiveRampLegacy";
}


void LHRepulsiveRampLegacy::finalize_setup(pose::Pose & pose ) {
	TR<<"   start finalize_setup function ..."<<std::endl;

	tf_= setup_packer_task(pose);


	( *dock_scorefxn_ )( pose );

	//setting MoveMap
	cdr_dock_map_ = kinematics::MoveMapOP( new kinematics::MoveMap() );
	cdr_dock_map_->clear();
	cdr_dock_map_->set_chi( false );
	cdr_dock_map_->set_bb( false );
	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
	utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );

	select_loop_residues( pose, all_loops_, false/*include_neighbors*/, bb_is_flexible);
	cdr_dock_map_->set_bb( bb_is_flexible );
	select_loop_residues( pose, all_loops_, true/*include_neighbors*/, sc_is_flexible);
	cdr_dock_map_->set_chi( sc_is_flexible );
	cdr_dock_map_->set_jump( 1, true );
	for ( Size ii = 2; ii <= all_loops_.num_loop() + 1; ii++ ) {
		cdr_dock_map_->set_jump( ii, false );
	}


	//set up sidechain movers for rigid body jump and loop & neighbors
	utility::vector1_size rb_jump;
	rb_jump.push_back( 1 );
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	// selecting movable c-terminal residues
	ObjexxFCL::FArray1D_bool loop_residues( pose.total_residue(), false );
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		loop_residues(i) = sc_is_flexible[i];
	} // check mapping

	using namespace protocols::toolbox::task_operations;
	tf_->push_back( TaskOperationCOP( new RestrictToInterface( rb_jump, loop_residues ) ) );


	TR<<"   finish finalize_setup function !!!"<<std::endl;

}


///////////////////////////////////////////////////////////////////////////
///
/// @brief ramping up the fullatom repulsive weight slowly to allow the
///        partners to relieve clashes and make way for each other
///
/// @details This routine is specially targetted to the coupled
///           optimization of docking partners and the loop region.  The
///           loop modelling & all previous  steps  involve mainly
///           centroid  mode .On switching  on fullatom mode, one is bound
///           to end up with clashes.To relieve the clashes, it is
///           essential to slowly  dial up the  repulsive weight instead of
///           turning it on to the maximum  value in one single step
///
/// @param[in] input pose which is assumed to have a docking fold tree
///
/// @global_read fa_rep : fullatom repulsive weight
///
/// @global_write fa_rep ( It is reset to the original value at the end )
///
/// @remarks A particular portion is  commented out,which can be
///          uncommented if one  uses a  low  resolution  homology  model.
///          Check details in the beginning of the commented out region
///
/// @references
///
/// @author Aroop 07/13/2010
///
///////////////////////////////////////////////////////////////////////////

void LHRepulsiveRampLegacy::apply( pose::Pose & pose ) {
	TR<<"start apply function ..."<<std::endl;


	finalize_setup(pose );

	// remove cutpoints variants for all cdrs
	// "true" forces removal of variants even from non-cutpoints
	loops::remove_cutpoint_variants( pose, true );
	using namespace core::chemical;
	for ( loops::Loops::const_iterator it = all_loops_.begin(),
			it_end = all_loops_.end(); it != it_end; ++it ) {
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, it->cut() );
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER,it->cut()+1);
	}

	// add scores to map
	( *dock_scorefxn_ )( pose );

	// dampen fa_rep weight
	core::Real rep_weight_max = dock_scorefxn_->get_weight( core::scoring::fa_rep );

	if ( benchmark_ ) {
		rep_ramp_cycles_ = 1;
		num_repeats_ = 1;
		min_threshold_ = 150.0;
	}


	core::Real rep_ramp_step = (rep_weight_max - 0.02) / core::Real(rep_ramp_cycles_-1);
	core::scoring::ScoreFunctionOP temp_scorefxn = dock_scorefxn_->clone();

	for ( Size i = 1; i <= rep_ramp_cycles_; i++ ) {
		core::Real rep_weight = 0.02 + rep_ramp_step * Real(i-1);
		TR<<"   repulsive ramp cycle "<<i<<":     rep_weight = "<<rep_weight<<std::endl;
		temp_scorefxn->set_weight( core::scoring::fa_rep, rep_weight );

		snugfit_MC_min(pose, temp_scorefxn);

	}

	TR<<"finish apply function !!!"<<std::endl;

}


//JQX: since the scorefxn needs to be changed, you have to
//     completely re-build the "simple_mcm_repeat" mover
void LHRepulsiveRampLegacy::snugfit_MC_min(pose::Pose & pose, core::scoring::ScoreFunctionOP temp_scorefxn) {
	using namespace moves;

	simple_moves::MinMoverOP
		min_mover( new simple_moves::MinMover( cdr_dock_map_, temp_scorefxn,  min_type_, min_threshold_, true/*nb_list*/ ) );

	rigid::RigidBodyPerturbMoverOP
		rb_perturb( new rigid::RigidBodyPerturbMover(pose, *cdr_dock_map_, rot_mag_, trans_mag_, rigid::partner_downstream, true ) );

	simple_moves::RotamerTrialsMoverOP
		pack_rottrial( new simple_moves::RotamerTrialsMover( pack_scorefxn_, tf_ ) );

	SequenceMoverOP rb_mover( new SequenceMover );
	rb_mover->add_mover( rb_perturb );
	rb_mover->add_mover( pack_rottrial );

	JumpOutMoverOP
		rb_mover_min( new JumpOutMover( rb_mover, min_mover, temp_scorefxn, min_threshold_ ) );


	MonteCarloOP
		mc( new MonteCarlo( pose, *temp_scorefxn, temperature_ ) );

	TrialMoverOP
		rb_mover_min_trial( new TrialMover( rb_mover_min, mc) );

	RepeatMoverOP
		simple_mcm_repeat( new RepeatMover( rb_mover_min_trial, num_repeats_ ) );


	simple_mcm_repeat->apply( pose );
}


void LHRepulsiveRampLegacy:: set_task_factory(pack::task::TaskFactoryCOP tf) {
	tf_ = pack::task::TaskFactoryOP( new pack::task::TaskFactory(*tf) );
}


} // namespace antibody
} // namespace protocols


