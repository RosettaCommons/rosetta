// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/RigidSearchMover.cc
///
/// @brief
/// @author Ian W. Davis

// IS THIS DEPRECATED? IS IT USED ANYWHERE?


#include <protocols/ligand_docking/RigidSearchMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <numeric/random/random.hh>

#include <cmath>

#include <utility/vector1.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

static THREAD_LOCAL basic::Tracer TR( "protocols.ligand_docking.RigidSearchMover", basic::t_debug );

namespace protocols {
namespace ligand_docking {


RigidSearchMover::RigidSearchMover(int jump_id, int num_trials, core::scoring::ScoreFunctionCOP scorefxn):
	Mover(),
	jump_id_(jump_id),
	num_trials_(num_trials),
	scorefxn_(scorefxn),
	temperature_(2),
	rotate_deg_(3),
	translate_Ang_(0.1),
	rotate_rsd_(0),
	rotate_atom_(0),
	recover_low_(true)
{
	Mover::type( "RigidSearch" );
}


RigidSearchMover::RigidSearchMover(RigidSearchMover const & that):
	//utility::pointer::ReferenceCount(),
	Mover(),
	jump_id_( that.jump_id_ ),
	num_trials_( that.num_trials_ ),
	scorefxn_( that.scorefxn_ ),
	temperature_( that.temperature_ ),
	rotate_deg_( that.rotate_deg_ ),
	translate_Ang_( that.translate_Ang_ ),
	rotate_rsd_( that.rotate_rsd_ ),
	rotate_atom_( that.rotate_atom_ ),
	recover_low_( that.recover_low_ )
{
}


RigidSearchMover::~RigidSearchMover()
{
}


void RigidSearchMover::apply(core::pose::Pose & pose)
{
	clock_t start_time = clock();
	using namespace core;
	using kinematics::Jump;
	using kinematics::Stub;
	Jump best_jump_so_far( pose.jump(jump_id_) );
	Jump last_accepted_jump( best_jump_so_far );
	Real best_score_so_far( (*scorefxn_)( pose ) );
	Real last_accepted_score( best_score_so_far );
	int num_accepts(0), num_improves(0);
	int const report_interval = std::min( 50, num_trials_ / 5 );
	TR << "Starting score " << best_score_so_far << std::endl;
	for ( int i = 1; i <= num_trials_; ++i ) {
		if ( i % report_interval == 0 ) {
			TR << "Cycle " << i << ", " << num_accepts << " accepts, " << num_improves << " improves, score = " << (*scorefxn_)( pose ) << std::endl;
		}
		// Do move (copied from RigidBodyPerturbMover)
		// Want to update our center of rotation every time we take a step.
		// Can either rotate around downstream centroid (default) or a specific atom.
		Vector dummy_up, rot_center;
		if ( rotate_rsd_ <= 0 || rotate_atom_ <= 0 ) protocols::geometry::centroids_by_jump(pose, jump_id_, dummy_up, rot_center);
		else rot_center = pose.residue(rotate_rsd_).xyz(rotate_atom_);
		Jump curr_jump = pose.jump( jump_id_ );
		// comments explain which stub to use when...
		Stub downstream_stub = pose.conformation().downstream_jump_stub( jump_id_ );
		curr_jump.set_rb_center( 1 /*n2c*/, downstream_stub, rot_center );
		curr_jump.gaussian_move( 1 /*n2c*/, translate_Ang_, rotate_deg_ );
		pose.set_jump( jump_id_, curr_jump );
		// score and do Boltzmann test
		Real const curr_score = (*scorefxn_)( pose );
		Real const deltaE = last_accepted_score - curr_score;
		if ( deltaE < 0 ) {
			// copied from MonteCarlo::boltzmann()
			Real const boltz = std::max(-40.0, deltaE / temperature_);
			if ( numeric::random::rg().uniform() >= std::exp(boltz) ) { // rejected!
				pose.set_jump( jump_id_, last_accepted_jump );
				continue;
			}
		}
		// accepted, thermally or otherwise
		num_accepts += 1;
		last_accepted_jump = pose.jump( jump_id_ );
		last_accepted_score = curr_score;
		if ( last_accepted_score < best_score_so_far ) {
			num_improves += 1;
			best_jump_so_far = last_accepted_jump;
			best_score_so_far = last_accepted_score;
		}
	}
	// recover lowest energy pose
	if ( recover_low_ ) pose.set_jump( jump_id_, best_jump_so_far );
	TR << "Best score " << best_score_so_far << ", end score " << (*scorefxn_)( pose ) << std::endl; // should be same!
	clock_t end_time = clock();
	TR << "Speed: " << num_trials_ / (double(end_time - start_time) / CLOCKS_PER_SEC) << " cycles per second" << std::endl;
}

std::string
RigidSearchMover::get_name() const {
	return "RigidSearchMover";
}


} // namespace ligand_docking
} // namespace protocols
