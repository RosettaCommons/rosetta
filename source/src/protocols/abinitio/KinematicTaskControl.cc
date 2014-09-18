// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KinematicTaskCenter
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @detailed responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/KinematicTaskControl.hh>

// Package Headers
#include <protocols/abinitio/ResolutionSwitcher.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/kinematics/util.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>

// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>


// AUTO-REMOVED #include <basic/options/option.hh>

#include <core/conformation/util.hh> //idealize

#include <protocols/loops/Loop.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/util.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>


//// C++ headers
// AUTO-REMOVED #include <fstream>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.general_abinitio", basic::t_info );

namespace protocols {
namespace abinitio {

using namespace core;

void
KinematicTaskControl::init( core::pose::Pose const & /* pose */ )  {
}

KinematicTaskControl::~KinematicTaskControl() {}

//@brief basic apply for the generalized protocol:
void KinematicTaskControl::apply( pose::Pose &pose ) {
	core::kinematics::simple_visualize_fold_tree( pose.fold_tree(), tr.Debug );
	tr.Debug << "KinematicTaskControl settings: "
					 << " return_full_atom " <<( !return_centroid() ? "yes" : "no" )
					 << " add side-chains " << ( sampling_protocol_->return_centroid() ? "yes" : "no" )
					 << std::endl;
	//resoltuion switch: if pose is fullatom ---> make it centroid but keep full-atom copy for the retrieval of side-chains later on.
	res_switch_ = new	ResolutionSwitcher(
			pose,
			b_input_is_fullatom_,
			sampling_protocol_->start_from_centroid(),
			sampling_protocol_->return_centroid()
	);
	//needs the full-atom scorefxn for repacking later
	res_switch().set_scorefxn( fullatom_scorefxn() ); //needed for repacking

	//-------------------------------------------------
	// the actual work
  bool success( true );
	// start pose ( centroid or fullatom based on the initialization of res_switch )
	pose = res_switch().start_pose();
	if ( tr.Debug.visible() ) output_debug_structure( pose, "start_pose" );
	//* do sampling:

	sampling_protocol_->set_current_job( get_current_job() );
	sampling_protocol_->set_current_tag( get_current_tag() );

	success = inner_loop( pose );

	set_current_tag( sampling_protocol_->get_current_tag() );

	if ( !success ) tr.Warning << "[WARNING] no success in sampler... could be the loop-closing " << std::endl;
	// apply res_switch to get back to a full-atom pose ( copies side-chains in unmoved regions )
	if ( tr.Debug.visible() ) output_debug_structure( pose, "before resolution switch ");

	if ( !return_centroid() ) res_switch().apply( pose );
	tr.Debug << "return from KinematicTaskControl" << std::endl;
	if (! (success && (res_switch().get_last_move_status() == moves::MS_SUCCESS) ) ) set_last_move_status( moves::FAIL_RETRY );
}

std::string
KinematicTaskControl::get_name() const {
	return "KinematicTaskControl";
}

//@brief sampling: simple version: get new kinematics (movemap+jumps) and call sampling protocol.
	//overwrite this guy if you want to do more stuff ... i.e., extend loops if things didn't work out in the first place.
bool KinematicTaskControl::inner_loop( core::pose::Pose& pose ) {
	bool success( false );

	Size fail( 0 );
	current_kinematics_ = NULL;
	while ( fail++ <= 10 && !current_kinematics_ ) {// get new setup
		//this may add constraints to the pose ...or should this be handled via the KinematicControl object?!
		current_kinematics_ = new_kinematics( pose );
	}

	//debug output
	if ( current_kinematics_ && tr.Info.visible() ) {
		tr.Info << "kinematic choice:\n";
		core::kinematics::simple_visualize_fold_tree_and_movemap(
				current_kinematics_->sampling_fold_tree(),
				current_kinematics_->movemap(),
				tr.Info );
		tr.Info << "\n" << jumping::JumpSample( current_kinematics_->sampling_fold_tree() );
		tr.Info << "\nfinal_fold-tree:\n";
		core::kinematics::simple_visualize_fold_tree( current_kinematics_->final_fold_tree(), tr.Info );
	}

	// if setup valid...
  if ( current_kinematics_ ) {
		// sample with this setup
		sampling_protocol_->set_kinematics( current_kinematics() );
		sampling_protocol_->apply( pose );
		success = ( sampling_protocol_->get_last_move_status() == moves::MS_SUCCESS );
	}
	// done...
	return success;
}

void KinematicTaskControl::set_extended_torsions_and_idealize_loops( core::pose::Pose& pose, loops::Loops loops ) const {

	// if no loops, we want to have extended structure everywhere.
	// it is a by-value parameter -- as intenden the change is kept local

	if ( loops.size() == 0 ) loops.add_loop( 1, pose.total_residue(), 0 );
	tr.Debug << "extend structure for " << loops << std::endl;
	for ( loops::Loops::const_iterator it = loops.begin(), eit = loops.end(); it != eit; ++it ) {
		Size const end_extended ( std::min( (int) it->stop(), (int) pose.total_residue() ) );
		for ( Size pos = std::max( 1, (int) it->start()); pos<=end_extended; pos++ ) {
			core::conformation::idealize_position( pos, pose.conformation() );
		}
 		Real const init_phi  ( -150.0 );
		Real const init_psi  (  150.0 );
		Real const init_omega(  180.0 );
		//special thing for residue 1 since idealize has a bug for residue 1
	// 	if ( it->start() == 1 ) {
// 			core::conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( pose.residue_type( 1 ) ) );
// 			pose.replace_residue( 1, *new_rsd , false /*orient backbone*/ );
// 		}

		for ( Size pos = it->start(); pos <= end_extended; pos++ ) {
			if( pos != it->start() )	pose.set_phi( pos,  init_phi );
			if( pos != end_extended ) pose.set_psi( pos,  init_psi );
			if( ( pos != it->start() ) && ( pos != end_extended ) ) pose.set_omega( pos,  init_omega );
		}
	}
}


}
}
