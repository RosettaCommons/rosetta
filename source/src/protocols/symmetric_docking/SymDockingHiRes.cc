// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SymDockingHiRes
/// @brief protocols that are specific to high resolution docking
/// @details
///		This contains the functions that create initial positions for docking
/// 	Also contains docking mcm protocol
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Sid Chaudhury
/// @author Modified by Jacob Corn
/// @author ingemar André. Based on the standard docking protocol

#include <protocols/symmetric_docking/SymDockingHiRes.hh>
#include <protocols/symmetric_docking/SymSidechainMinMover.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>


#include <basic/options/option.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/Residue.hh> // for design() flag
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymDof.hh>

#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/RepeatMover.hh>
//for resfile reading
#include <basic/options/keys/packing.OptionKeys.gen.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>


#include <basic/Tracer.hh>
using basic::T;

// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.symmetric_docking.SymDockingHiRes" );

using namespace core;

namespace protocols {
namespace symmetric_docking {

// default constructor
SymDockingHiRes::SymDockingHiRes() : Mover()
{
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	scorefxn_pack_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	moves::Mover::type( "SymDockingHiRes" );
	init_task_factory_=NULL;
	design_ = false;
}

// constructor with arguments
SymDockingHiRes::SymDockingHiRes(
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::scoring::ScoreFunctionOP scorefxn_pack_in
) : Mover(), scorefxn_(scorefxn_in), scorefxn_pack_(scorefxn_pack_in)
{
	moves::Mover::type( "SymDockingHRes" );
	init_task_factory_=NULL;
	design_ = false;
}

//destructor
SymDockingHiRes::~SymDockingHiRes() {}

//clone
protocols::moves::MoverOP SymDockingHiRes::clone() const {
	return protocols::moves::MoverOP( new SymDockingHiRes(*this) );
}
// what type of minimization is used?
void SymDockingHiRes::set_min_type( std::string min_type_in ) { min_type_ = min_type_in;}

// repack sidechains during docking?
void SymDockingHiRes::set_repack( bool repack_switch){ repack_switch_ = repack_switch;}

// pointer to mc object
moves::MonteCarloOP SymDockingHiRes::get_mc() { return mc_; }

// set the packer task
void
SymDockingHiRes::task_factory( core::pack::task::TaskFactoryOP task )
{
	init_task_factory_ = task;
}

core::pack::task::TaskFactoryOP &
SymDockingHiRes::task_factory(){
	return init_task_factory_;
}

// Use design during docking. Currently not tested...
void SymDockingHiRes::design( bool const des ) {
	design_ = des;
}

// Are we designing during the docking procedure
bool SymDockingHiRes::design() const { return design_; }

// Set default values for the docking protocol
void SymDockingHiRes::set_default( core::pose::Pose & pose ) {
		using namespace basic::options;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		// SETS UP the stuff in pose
		(*scorefxn_)( pose );

		// Get the docking translation and rotation magnitutes from command line. If
		// none is given the default values of 0.1 A and 5 degrees are used
		trans_magnitude_ = option[ OptionKeys::docking::dock_mcm_trans_magnitude ]();
		rot_magnitude_ = option[ OptionKeys::docking::dock_mcm_rot_magnitude ]();
		// use rtmin_ and sc_min. rtmin for symmetry is not implemented?! Turned off 		// for now. The flag is still here which adds to the confusion.

		rtmin_ = option[ OptionKeys::docking::dock_rtmin ]();
		scmin_ = option[ OptionKeys::docking::sc_min ]();

		temperature_ = 0.8;
		repack_switch_ = true;
		repack_period_ = 8;

		//sets up MC object
		mc_ = moves::MonteCarloOP( new moves::MonteCarlo( pose, *scorefxn_, temperature_ ) );

		//sets up default movemap
		bb_ = false;
		chi_ = false;
		movemap_ = core::kinematics::MoveMapOP( new kinematics::MoveMap() );
		movemap_->set_chi( chi_ );
		movemap_->set_bb( bb_ );
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );

		//sets up minimization parameters
		min_tolerance_ = 0.01;
		min_type_ = std::string( "dfpmin_armijo_nonmonotone" );
		nb_list_ = true;

		setup_packing( pose );

		//set up docking protocol based on options
		set_protocol( pose );
	}

// Set the movemap
void SymDockingHiRes::set_move_map(core::kinematics::MoveMapOP movemap_in){ movemap_ = movemap_in; }

// define which protocol is used...
void SymDockingHiRes::set_protocol( pose::Pose & pose ){
	using namespace basic::options;

	if ( option[ OptionKeys::docking::dock_min ]() ) {
		set_dock_min_protocol();
		}
		else if ( option[ OptionKeys::docking::dock_ppk ]() ){
			set_dock_ppk_protocol( pose );
		}

	else {
	  set_dock_mcm_protocol( pose );
		}
}

/*void SymDockingHiRes::define_loops( pose::Pose const & pose, loops::Loops & loop_set, Real & interface_dist ) {
	//runtime_assert( movable_jumps_.size() == 1 ); // CURRENTLY ONLY SUPPORTED WITH SIMPLE DOCKING
	//core::Size const rb_jump = movable_jumps_[1];

	loop_set.clear();

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	//RestrictTaskForDockingOP rtfd = new RestrictTaskForDocking( scorefxn_, rb_jump_, true, interface_dist );
	pack::task::TaskFactory tf;
	//tf.push_back( rtfd );
	//tf.push_back( new RestrictTaskForDocking( scorefxn_, rb_jump_, true, interface_dist ) );
	tf.push_back( new RestrictToInterface( movable_jumps_, interface_dist ) );
	pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );

	// extend one residue beyond borders of repackable regions, don't allow 1-residue loops
	core::Size const nres = pose.total_residue();
	utility::vector1<bool> flexible_region( nres, false );
	for ( Size i=2; i < nres; ++i ) {
		int num_flexible(0);
		if ( pose.pdb_info()->chain(i-1) == pose.pdb_info()->chain(i+1) ) {
			if ( task->pack_residue(i-1) ) ++num_flexible;
			if ( task->pack_residue(i) ) ++num_flexible;
			if ( task->pack_residue(i+1) ) ++num_flexible;
		}
		if ( num_flexible > 1 ) {
			flexible_region.at(i-1) = true;
			flexible_region.at(i) = true;
			flexible_region.at(i+1) = true;
		}
	}

	// if we have a single fixed residue between two loops, make this flexible
	for ( Size i=2; i < nres; ++i ) {
		if ( flexible_region.at(i-1) && flexible_region.at(i+1) ) flexible_region.at(i) = true;
	}

	// jk For now, don't let the first or last two residues of a chain be flexible
	flexible_region.at(1) = false;
	flexible_region.at(2) = false;
	for ( Size i=3; i < (nres-1); ++i ) {
		if ( pose.pdb_info()->chain(i-1) != pose.pdb_info()->chain(i) ) {
			flexible_region.at(i-2) = false;
			flexible_region.at(i-1) = false;
			flexible_region.at(i) = false;
			flexible_region.at(i+1) = false;
		}
	}
	flexible_region.at(nres-1) = false;
	flexible_region.at(nres) = false;

	// disallow one-residue loops
	for ( Size i=2; i < nres; ++i ) {
		if ( ( ! flexible_region.at(i-1)) && ( ! flexible_region.at(i+1)) ) flexible_region.at(i) = false;
	}
	// disallow two-residue loops
	for ( Size i=3; i < nres; ++i ) {
		if ( ( ! flexible_region.at(i-2)) && ( ! flexible_region.at(i+1)) ) {
			flexible_region.at(i-1) = false;
			flexible_region.at(i) = false;
		}
	}

	// setup loops
	core::Size loop_start=0;
	core::Size loop_stop=0;
	for ( Size i=1; i < nres; ++i ) {
		if ( flexible_region.at(i) ) {
			loop_start = i;
			loop_stop = i;
			for ( Size j=i+1; j <= nres; ++j ) {
				// if j is on a different chain than i, break
				if ( pose.pdb_info()->chain(j) != pose.pdb_info()->chain(i) ) {
					break;
				}
				// if j is not flexible, break
				if ( ( ! flexible_region.at(j) ) ) {
					break;
				}
				loop_stop = j;
			}
			loop_set.add_loop( loop_start, loop_stop, 0 );
			i = loop_stop;
		}
	}

	loop_set.choose_cutpoints( pose );

	return;
} */


////////////////////////////////////////////////////////////////////////////////
/// @brief
/// @details
///		decides what to call according to options
void SymDockingHiRes::apply( core::pose::Pose & pose )
{
	using namespace scoring;
	using namespace basic::options;

	TR << "in SymDockingHiRes.apply" << std::endl;

	// jec sanity check to avoid overwriting newly-set minimizers on every apply
	if ( !mc_ ) {
		set_default( pose );
	}

	mc_->reset( pose );

	docking_highres_protocol_mover_->apply( pose );
	mc_->recover_low( pose );
}

std::string
SymDockingHiRes::get_name() const {
	return "SymDockingHiRes";
}

///////////////////////////////////////////////////////////////////////////////////
///
/// @brief main entrance for normal rigid-body minimization
/// @details
///		retrieve the structure in the low array and do the normal minimization
///		by calling using a min_mover to optimize the score accourding to the
///		scorefunction that has been set
///
/// @remarks
///
/// @references docking_minimize_trial from docking_minimize.cc
///				pose_docking_minimize_trial from pose_docking.cc
///
/// @author Monica Berrondo June 14 2007, modified for symmetric docking by
///					 Ingemar Andre
///
/////////////////////////////////////////////////////////////////////////////////
void SymDockingHiRes::set_dock_min_protocol() {
	using namespace moves;

	TR << "::::::::::::::::::DOCK_MIN:::::::::::::::::::" << std::endl;

	protocols::simple_moves::MinMoverOP min_mover( new simple_moves::symmetry::SymMinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	TrialMoverOP minimize_trial( new TrialMover( min_mover, mc_ ) );
	docking_highres_protocol_mover_ = moves::SequenceMoverOP( new SequenceMover );
	docking_highres_protocol_mover_->add_mover( minimize_trial );
}

///////////////////////////////////////////////////////////////////////////////////
///
/// @brief main entrance to do monte carlo minimization
/// @details
///			a total of 50 cycles of monte-carlo minimization will be
///			carried out if the minimized structure can pass the filter
///			after the first and fifth cycle.  Then it is rigid-body minimized
///			to a stringent tolerance.
///
/// @remarks
///
/// @references docking_mcm_protocol from docking_minimize.cc
///				pose_docking_monte_carlo_minimize from pose_docking.cc
///
/// @author Sid Chaudhury May 28 2009, modified for symmetric docking
///					 by Ingemar Andre
///
/////////////////////////////////////////////////////////////////////////////////
void SymDockingHiRes::set_dock_mcm_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::conformation::symmetry;
	using namespace protocols::toolbox::task_operations;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	//set up rigid body movers
	rigid::RigidBodyDofSeqPerturbMoverOP rb_perturb( new rigid::RigidBodyDofSeqPerturbMover( dofs , rot_magnitude_, trans_magnitude_ ) );

	//set up minimizer movers
	protocols::simple_moves::MinMoverOP min_mover( new simple_moves::symmetry::SymMinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ ) );

	//set up sidechain movers for each movable jump
	//tf_->push_back( new RestrictToInterface( 1 ) );

	//fpd -- repack any interface formed by a movable jump
	utility::vector1<int> movable_jumps;
	for (std::map<Size,SymDof>::iterator i=dofs.begin(),i_end=dofs.end(); i!=i_end; ++i) {
		movable_jumps.push_back( i->first );
	}
	tf_->push_back( TaskOperationCOP( new RestrictToInterface( movable_jumps ) ) );

	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new simple_moves::symmetry::SymRotamerTrialsMover( scorefxn_pack_, tf_ ) );

	SequenceMoverOP interface_repack_and_move_loops( new moves::SequenceMover );

	std::string const flex_bb_docking_type = option[ OptionKeys::docking::flexible_bb_docking ]();
	if ( flex_bb_docking_type == "fixedbb" ) {
		// Call pack_rotamers, no backbone movement
		protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new simple_moves::symmetry::SymPackRotamersMover( scorefxn_pack_ ) );
		pack_interface_repack->task_factory(tf_);
		interface_repack_and_move_loops->add_mover( pack_interface_repack );
	} else {

		// Call pack_rotamer before and after loop movement
		protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new simple_moves::symmetry::SymPackRotamersMover( scorefxn_pack_ ) );
		pack_interface_repack->task_factory(tf_);
		interface_repack_and_move_loops->add_mover( pack_interface_repack );

		// This does not work for symmetry yet. Anyone is free to implement it...
		if ( ( flex_bb_docking_type == "ccd" ) || ( flex_bb_docking_type == "kic" ) ||  ( flex_bb_docking_type == "backrub" ) ) {
			TR << "[ ERROR ] flexible_bb_docking is not implemented for symmetric docking yet..." << std::endl;
			exit(1);
		} else {
			TR << "[ ERROR ] Unknown flexible_bb_docking type: " << flex_bb_docking_type << std::endl;
			exit(1);
		}

		interface_repack_and_move_loops->add_mover( pack_interface_repack );
	}


	TrialMoverOP pack_interface_and_move_loops_trial( new TrialMover( interface_repack_and_move_loops, mc_ ) );

//	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin = new protocols::simple_moves::RotamerTrialsMinMover( scorefxn_pack_, tf_ );
//	TrialMoverOP rtmin_trial = new TrialMover( rtmin, mc_ );

	//InterfaceSidechainMinMoverOP scmin_mover = new InterfaceSidechainMinMover(rb_jump_, scorefxn_pack_ );
	SymSidechainMinMoverOP scmin_mover( new SymSidechainMinMover(scorefxn_pack_, tf_ ) );
	TrialMoverOP scmin_trial( new TrialMover( scmin_mover, mc_ ) );

	// the standard mcm cycle : rb perturbation->rotamer trials->minimization->MC accept
	SequenceMoverOP rb_mover( new SequenceMover );
	rb_mover->add_mover( rb_perturb );
	if ( repack_switch_ ) rb_mover->add_mover( pack_rottrial );

	core::Real minimization_threshold = 15.0;
	JumpOutMoverOP rb_mover_min( new JumpOutMover( rb_mover, min_mover, scorefxn_, minimization_threshold ) );
	TrialMoverOP rb_mover_min_trial( new TrialMover( rb_mover_min, mc_ ) );

	//every step (rb_mover_min_trial): the standard mcm cycle
	//every 8th step (repack_step): standard mcm cycle + repacking
	//try moving loops too, if desired

	SequenceMoverOP repack_step( new SequenceMover );
	repack_step->add_mover(rb_mover_min_trial);

	if ( repack_switch_ ){
		repack_step->add_mover( pack_interface_and_move_loops_trial );
	//	if (rtmin_) repack_step->add_mover(rtmin_trial);
		if (scmin_) repack_step->add_mover(scmin_trial);
		}

	CycleMoverOP rb_mover_min_trial_repack( new CycleMover );
	for ( Size i=1; i<repack_period_; ++i ) rb_mover_min_trial_repack->add_mover( rb_mover_min_trial );
	rb_mover_min_trial_repack->add_mover( repack_step );

	//set up initial repack mover
	SequenceMoverOP initial_repack( new SequenceMover );
	initial_repack->add_mover(pack_interface_and_move_loops_trial);
	//if (rtmin_) initial_repack->add_mover(rtmin_trial);
	if (scmin_) initial_repack->add_mover(scmin_trial);

	//set up initial and final min_trial movers for docking
	TrialMoverOP minimize_trial( new TrialMover( min_mover, mc_ ) );

	//set up mcm cycles and mcm_repack cycles
	RepeatMoverOP mcm_four_cycles( new RepeatMover( rb_mover_min_trial, 4 ) );
	RepeatMoverOP mcm_fortyfive_cycles( new RepeatMover( rb_mover_min_trial_repack, 45 ) );
	//set up protocol mover
	TR << "::::::::::::::::::DOCK_MCM:::::::::::::::::::" << std::endl;

	docking_highres_protocol_mover_ = moves::SequenceMoverOP( new SequenceMover );
	if (repack_switch_) docking_highres_protocol_mover_->add_mover( initial_repack );
	docking_highres_protocol_mover_->add_mover( minimize_trial );
	docking_highres_protocol_mover_->add_mover( mcm_four_cycles );
	docking_highres_protocol_mover_->add_mover( mcm_fortyfive_cycles );
	docking_highres_protocol_mover_->add_mover( minimize_trial );

}

// noncost pose for load_unboundrot csts
// detup the options for packing including setup up the packer task
void SymDockingHiRes::setup_packing( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	//set upconstructor packer options
	tf_ = core::pack::task::TaskFactoryOP( new TaskFactory );
	if( init_task_factory_ )
	{
		TR << "Using user-defined TaskFactory." << std::endl;
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory( *init_task_factory_ ) );
	}
	if( design_ ) {
		TR << "Designing during docking" << std::endl;
	}
	else { // default case -- restrict everything to repacking.
		tf_->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	}
//	tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
	tf_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf_->push_back( TaskOperationCOP( new IncludeCurrent ) );
	tf_->push_back( TaskOperationCOP( new NoRepackDisulfides ) );
	if( option[OptionKeys::packing::resfile].user() ) tf_->push_back( TaskOperationCOP( new ReadResfile ) );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot( new core::pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation( new operation::AppendRotamerSet( unboundrot ) );
	tf_->push_back( unboundrot_operation );
	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	// note that RestrictToInterfaceOperation is added during set_dock_mcm_protocol
}

/// define the prepacking protocol
void SymDockingHiRes::set_dock_ppk_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::conformation::symmetry;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	TR << "::::::::::::::::::DOCK_PPK:::::::::::::::::::" << std::endl;

	//set up translate-by-axis movers
	rigid::RigidBodyDofSeqTransMoverOP translate_away( new rigid::RigidBodyDofSeqTransMover( dofs ) );
	rigid::RigidBodyDofSeqTransMoverOP translate_back( new rigid::RigidBodyDofSeqTransMover( dofs ) );
	translate_away->step_size( 1000 );
	translate_back->step_size( -1000 );


	PackerTaskOP task = tf_->create_task_and_apply_taskoperations( pose ); // does not include restrict to interface

	protocols::simple_moves::PackRotamersMoverOP prepack_full_repack( new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn_pack_, task ) );
	//RotamerTrialsMinMoverOP rtmin_mover = new symmetry::SymRotamerTrialsMinMover( scorefxn_pack_, *task );
	SymSidechainMinMoverOP scmin_mover( new SymSidechainMinMover(scorefxn_pack_, task) );

	// set up protocol
	docking_highres_protocol_mover_ = moves::SequenceMoverOP( new SequenceMover );
	if (scmin_) docking_highres_protocol_mover_->add_mover( scmin_mover );
	docking_highres_protocol_mover_->add_mover( translate_away );
	docking_highres_protocol_mover_->add_mover( prepack_full_repack );
	//if (rtmin_) docking_highres_protocol_mover_->add_mover( rtmin_mover );
	if (scmin_) docking_highres_protocol_mover_->add_mover( scmin_mover );
	docking_highres_protocol_mover_->add_mover( translate_back );
}


} // namespace docking
} // namespace protocols
