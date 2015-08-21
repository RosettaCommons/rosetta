// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FullatomRelaxMover.cc
/// @author Robin A Thottungal  (rathottungal@gmail.com)
/// @author Michael Pacella (mpacella88@gmail.com)

//Unit header
#include <protocols/surface_docking/FullatomRelaxMover.hh>

// Package header
#include <protocols/surface_docking/SurfaceOrientMover.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>

// Project headers
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RB_geometry.hh>

// Numeric Headers
#include <numeric/random/random.hh>

//Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/Tracer.hh>
// Utility Headers
#include <utility/exit.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.SurfaceDocking.FullatomRelaxMover" );

namespace protocols {
namespace surface_docking {

using namespace core;
using protocols::jd2::JobDistributor;
using namespace protocols::moves;

//constructor
FullatomRelaxMover::FullatomRelaxMover() : Mover()
{
	setup_defaults();
}

FullatomRelaxMover::FullatomRelaxMover( FullatomRelaxMover const & src) : Mover(src)
{
	copy_data(*this, src);
}

//destructor
FullatomRelaxMover::~FullatomRelaxMover() {}

protocols::moves::MoverOP
FullatomRelaxMover::clone() const
{
	return protocols::moves::MoverOP( new FullatomRelaxMover(*this) );
}

protocols::moves::MoverOP
FullatomRelaxMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new FullatomRelaxMover() );
}

void FullatomRelaxMover::setup_defaults()
{
	score_high_res_ = scoring::get_score_function();
	score_high_res_->set_weight( core::scoring::fa_elec, 0.25 );
	small_min_type_="linmin";
	shear_min_type_="dfpmin";
	nmoves_ = 6;
	kT_ = 0.5;
	Mover::type("FullatomRelaxMover");

	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
		encounter_cycle_ = 1;
	} else {
		encounter_cycle_= numeric::random::rg().random_range(1,5);
	}

	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
		outer_loop_cycles_ = 1;
		inner_loop_cycles_ = 1;
	} else {
		outer_loop_cycles_ = 5;
		inner_loop_cycles_ = 5;
	}
	TR << "Random Number Generated"<< std::endl;
	TR<< "Adsorption will occur after "<<encounter_cycle_<<" cycles..."<<std::endl;
}

void FullatomRelaxMover::setup_movers( const core::pose::Pose & pose )
{
	// Setting Common Parameters for movers
	Size const first_protein_residue = pose.num_jump() + 1;

	move_map_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
	move_map_->set_bb_true_range(first_protein_residue , pose.total_residue());
	move_map_->set_chi_true_range( first_protein_residue , pose.total_residue() );
	move_map_->set_jump(pose.num_jump(),false);

	Real tolerance=0.001;
	small_mover_ = simple_moves::SmallMoverOP( new simple_moves::SmallMover(move_map_,kT_,nmoves_) );
	small_min_mover_ = simple_moves::MinMoverOP( new simple_moves::MinMover(move_map_,
		score_high_res_,small_min_type_,tolerance,true,false,false) );

	small_sequence_mover_ = moves::SequenceMoverOP( new moves::SequenceMover() );
	small_sequence_mover_->add_mover(small_mover_);
	small_sequence_mover_->add_mover(small_min_mover_);


	shear_mover_ = simple_moves::ShearMoverOP( new simple_moves::ShearMover(move_map_,kT_,nmoves_/2) );

	shear_min_mover_ = simple_moves::MinMoverOP( new simple_moves::MinMover(move_map_,score_high_res_,
		shear_min_type_,tolerance,true,false,false) );

	shear_sequence_mover_ = moves::SequenceMoverOP( new moves::SequenceMover() );
	shear_sequence_mover_->add_mover(shear_mover_);
	shear_sequence_mover_->add_mover(shear_min_mover_);

	monte_carlo_ = moves::MonteCarloOP( new moves::MonteCarlo(pose,*score_high_res_,kT_) );
	//smallTrialMover
	small_trial_min_mover_ = moves::TrialMoverOP( new moves::TrialMover(small_sequence_mover_,monte_carlo_) );
	//shearTrialMover
	shear_trial_min_mover_ = moves::TrialMoverOP( new moves::TrialMover(shear_sequence_mover_,monte_carlo_) );
	dock_mcm_ = docking::DockMCMProtocolOP( new docking::DockMCMProtocol(pose.num_jump()) );
}

void FullatomRelaxMover::set_smallmovesize(Size scale){
	small_mover_->angle_max('H',scale/3);
	small_mover_->angle_max('E',scale/2);
	small_mover_->angle_max('L',scale);
}

void FullatomRelaxMover::set_ljrepulsion_weight(Real weight_scale){
	score_high_res_->set_weight(core::scoring::fa_rep,weight_scale);
}

void FullatomRelaxMover::set_ecounter(Size ecount){
	encounter_=ecount;
}

void FullatomRelaxMover::set_surface_contact_mover(protocols::docking::FaDockingSlideIntoContactOP surface_contact_mover)
{
	surface_contact_mover_ = surface_contact_mover;
}

void FullatomRelaxMover::set_surface_orient_mover( SurfaceOrientMoverOP surface_orient )
{
	surface_orient_ = surface_orient;
}

void FullatomRelaxMover::inner_loop_refinement( core::pose::Pose & pose )
{
	//k<=5 original value;k=2 test
	small_min_mover_->min_type("linmin");
	shear_min_mover_->min_type("linmin");
	set_secondary_struct(pose);
	small_trial_min_mover_->apply(pose); //linmin
	set_secondary_struct(pose);
	shear_trial_min_mover_->apply(pose);  //linmin
	//crank mover is missing
}

void FullatomRelaxMover::output_solution_state( core::pose::Pose & pose )
{
	TR<<"Solution State Structure Found"<<std::endl;
	//Write the solutionState Structure
	//getting lowest energy pose
	monte_carlo_->recover_low(pose);
	monte_carlo_->show_state();
	monte_carlo_->show_counters();
	TR<<"Calculating Secondary Structure of Solution State..."<<std::endl;
	calc_secondary_struct(pose);
	// Creating a job compatible with JD2
	protocols::jd2::JobOP job2
		=jd2::JobDistributor::get_instance()->current_job();
	//std::string job_name (JobDistributor::get_instance()->
	//        job_outputter()->output_name( job2 ) );
	job2->add_string_real_pair("Total weighted score: ", pose.energies().total_energy());
	job2->add_string_string_pair("SolState_SecondaryStructure:",sec_struct_);
	JobDistributor::get_instance()->job_outputter()->
		other_pose( job2,pose, "SolState_");
	sol_sec_struct_ = sec_struct_;

}

void FullatomRelaxMover::reorient_and_slide_into_surface( core::pose::Pose & pose )
{
	// Random Orient the Partner (make sure this is the peptide)
	TR<<"Preparing to dock protein to surface"<<std::endl;
	TR<<"Randomizing orientation..."<<std::endl;
	Size rb_jump_=pose.num_jump(); //default value
	rigid::RigidBodyRandomizeMover rmover( pose, rb_jump_, rigid::partner_upstream );
	rmover.apply( pose );
	//Axis Spin
	rigid::RigidBodySpinMover smover( rb_jump_ );
	smover.apply( pose );
	TR<<"Repositioning protein above surface..."<<std::endl;
	reposition_above_surface(pose);
	// SurfaceOrient Mover
	TR<<"Reorienting protein to central unit cell..."<<std::endl;
	surface_orient_->apply(pose);
	TR<<"Sliding protein into contact with the surface..."<<std::endl;
	surface_contact_mover_->apply( pose );
	monte_carlo_->reset(pose);
	TR<<"Protein adsorbed onto the Surface..."<<std::endl;

}

void FullatomRelaxMover::dock_mcm_on_surface(core::pose::Pose & pose)
{
	TR<<"Starting high-resolution docking on surface..."<<std::endl;
	//move_map_->set_chi(true);

	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
		dock_mcm_->set_first_cycle( 1 );
		dock_mcm_->set_second_cycle( 1 );
	}
	dock_mcm_->apply(pose);
}

void FullatomRelaxMover::outer_loop_refinement_solution(core::pose::Pose & pose)
{
	small_min_mover_->min_type("dfpmin");
	shear_min_mover_->min_type("dfpmin");
	set_secondary_struct(pose);
	small_trial_min_mover_->apply(pose); // dfpmin
	set_secondary_struct(pose);
	shear_trial_min_mover_->apply(pose);  //dfpmin
	for ( Size k = 1; k <= inner_loop_cycles_; ++k ) {
		inner_loop_refinement( pose );
	}

}

void FullatomRelaxMover::outer_loop_refinement_adsorbed(core::pose::Pose & pose)
{
	surface_orient_->apply(pose);
	TR<<"outer_loop_refinement: sliding protein into contact"<<std::endl;
	surface_contact_mover_->apply( pose );
	//reposition_above_surface( pose );
	outer_loop_refinement_solution(pose);
}

void FullatomRelaxMover::reposition_above_surface(core::pose::Pose & pose)
{
	Vector protein_centroid, surf_centroid;
	Vector slide_into, slide_away;
	Size const rb_jump=pose.num_jump();
	protocols::geometry::centroids_by_jump (pose, rb_jump, surf_centroid, protein_centroid);

	slide_away = surface_parameters_->slide_axis().negated();

	//getting a point 30 angstroms above surface centroid
	Vector point_above = surf_centroid+slide_away.normalized()*80;
	//vector needed to move protein centroid to point above
	Vector position_above = point_above-protein_centroid;
	protocols::rigid::RigidBodyTransMover position_above_surface = protocols::rigid::RigidBodyTransMover(position_above, pose.num_jump());
	position_above_surface.step_size(position_above.magnitude());
	position_above_surface.apply(pose);
}

void FullatomRelaxMover::refinement_cycle(pose::Pose & pose)
{
	monte_carlo_->reset(pose);
	calc_secondary_struct(pose);
	TR<<"Fullatom refinement cycle: "<<encounter_<<std::endl;
	TR<<"Adsoprtion occuring after: "<<encounter_cycle_<<" cycles"<<std::endl;
	// object for slide into contact
	for ( Size j = 1; j <=outer_loop_cycles_; ++j ) {
		if ( encounter_ <= encounter_cycle_ ) {
			outer_loop_refinement_solution(pose);
		} else {
			outer_loop_refinement_adsorbed(pose);
		}
	}
	if ( encounter_ == encounter_cycle_ ) {
		output_solution_state(pose);
		move_map_->set_jump(pose.num_jump(),true);
		reorient_and_slide_into_surface(pose);
	}
	if ( encounter_ > encounter_cycle_ ) { //I think this should be greater than to match the algorithm description (addressed)
		dock_mcm_on_surface(pose);
	}
}

void FullatomRelaxMover::apply(core::pose::Pose & pose)
{
	setup_movers(pose);
	core::Size lj_ramp_cycle=5;
	core::Size total_cycles=7;
	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
		lj_ramp_cycle = 2;
		total_cycles = 2;
	}
	core::Real lj_increment = ( .44 - 0.02 )/ (lj_ramp_cycle-1);
	for ( Size i=1; i<=total_cycles; ++i ) {
		set_smallmovesize(30/i);
		if ( i<=5 ) {
			set_ljrepulsion_weight((0.02+(i-1)*lj_increment));
		}
		set_ecounter(i); // when value of i matches the random number
		// in the allatomrelax, protein is slide into the surface
		refinement_cycle(pose);
	}
	// Final Side Chain re-packing using rtmin; adopted from Dave's protocol
	calc_secondary_struct(pose);
	ads_sec_struct_ = sec_struct_;
}

void FullatomRelaxMover::set_nmoves( core::Size const nmoves_in ){
	nmoves_ = nmoves_in;
}

void FullatomRelaxMover::set_surface_parameters(protocols::surface_docking::SurfaceParametersOP surface_parameters)
{
	surface_parameters_ = surface_parameters;
}
std::string FullatomRelaxMover::get_name() const {
	return "FullatomRelaxMover";
}

void FullatomRelaxMover::calc_secondary_struct(core::pose::Pose & pose){
	// creating a new pose and splitting!
	//split the pose into separate chains
	sec_struct_ = "";
	pose::Pose pose_tmp= pose::Pose(pose);
	utility::vector1< pose::PoseOP > singlechain_poses;
	singlechain_poses = pose_tmp.split_by_chain();
	core::scoring::dssp::Dssp dssp( *singlechain_poses[2] );
	Size const last_surface_residue( pose.num_jump());
	for ( Size ii = 1; ii <= singlechain_poses[2]->total_residue(); ++ii ) {
		if ( dssp.get_dssp_secstruct(ii) == ' ' ) {
			pose.set_secstruct(last_surface_residue + ii, 'L');

		} else {
			pose.set_secstruct(last_surface_residue + ii, dssp.get_dssp_secstruct(ii));
		}
		sec_struct_+=dssp.get_dssp_secstruct(ii);
	}
}

void FullatomRelaxMover::set_secondary_struct(core::pose::Pose & pose)
{
	Size const last_surface_residue( pose.num_jump());
	for ( Size ii = 1; ii <= sec_struct_.length(); ++ii ) {
		if ( sec_struct_[ii] == ' ' ) {
			pose.set_secstruct(last_surface_residue + ii, 'L');

		} else {
			pose.set_secstruct(last_surface_residue + ii, sec_struct_[ii]);
		}
	}
}

std::string FullatomRelaxMover::get_sol_secondary_struct(){
	return sol_sec_struct_;
}

std::string FullatomRelaxMover::get_ads_secondary_struct(){
	return ads_sec_struct_;
}

void FullatomRelaxMover::copy_data(FullatomRelaxMover object_to_copy_to, FullatomRelaxMover object_to_copy_from)
{
	object_to_copy_to.kT_ = object_to_copy_from.kT_;
	object_to_copy_to.nmoves_ = object_to_copy_from.nmoves_;
	object_to_copy_to.encounter_ = object_to_copy_from.encounter_;
	object_to_copy_to.encounter_cycle_ = object_to_copy_from.encounter_cycle_;
	object_to_copy_to.angle_max_ = object_to_copy_from.angle_max_;
	object_to_copy_to.score_high_res_ = object_to_copy_from.score_high_res_;
	object_to_copy_to.small_min_type_ = object_to_copy_from.small_min_type_;
	object_to_copy_to.shear_min_type_ = object_to_copy_from.shear_min_type_;
	object_to_copy_to.move_map_ = object_to_copy_from.move_map_;
	object_to_copy_to.monte_carlo_ = object_to_copy_from.monte_carlo_;
	object_to_copy_to.small_mover_ = object_to_copy_from.small_mover_;
	object_to_copy_to.small_min_mover_ = object_to_copy_from.small_min_mover_;
	object_to_copy_to.small_sequence_mover_ = object_to_copy_from.small_sequence_mover_;
	object_to_copy_to.small_trial_min_mover_ = object_to_copy_from.small_trial_min_mover_;
	object_to_copy_to.shear_mover_ = object_to_copy_from.shear_mover_;
	object_to_copy_to.shear_min_mover_ = object_to_copy_from.shear_min_mover_;
	object_to_copy_to.shear_sequence_mover_ = object_to_copy_from.shear_sequence_mover_;
	object_to_copy_to.shear_trial_min_mover_ = object_to_copy_from.shear_trial_min_mover_;
	object_to_copy_to.sol_sec_struct_ = object_to_copy_from.sol_sec_struct_;
	object_to_copy_to.ads_sec_struct_ = object_to_copy_from.ads_sec_struct_;
	object_to_copy_to.sec_struct_ = object_to_copy_from.sec_struct_;
	object_to_copy_to.surface_contact_mover_ = object_to_copy_from.surface_contact_mover_;
	object_to_copy_to.surface_orient_ = object_to_copy_from.surface_orient_;
	object_to_copy_to.dock_mcm_ = object_to_copy_from.dock_mcm_;
	object_to_copy_to.outer_loop_cycles_ = object_to_copy_from.outer_loop_cycles_;
	object_to_copy_to.inner_loop_cycles_ = object_to_copy_from.inner_loop_cycles_;
	object_to_copy_to.surface_parameters_ = object_to_copy_from.surface_parameters_;
}


} //surface_docking

} //protocol
