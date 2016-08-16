// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/mmt_msd/MMTMinPackingJob.hh
/// @brief  declaration for class MMTMinPackingJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <devel/mmt_msd/MMTMinPackingJob.hh>

#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationGraph.hh>

#include <core/pack/min_pack.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>
#include <core/pack/scmin/SidechainStateAssignment.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/optimization/MinimizerOptions.hh>

#include <utility/excn/Exceptions.hh>

namespace devel {
namespace mmt_msd {

MMTMinPackingJob::MMTMinPackingJob() :
	MMTPackingJob(),
	cartesian_( false ),
	nonideal_( false )
{}

MMTMinPackingJob::~MMTMinPackingJob() {}

void MMTMinPackingJob::setup()
{
	assert( has_pose() );
	assert( has_sfxn() );
	assert( has_task() );

	core::pack::min_pack_setup(
		pose(), sfxn(), task(), cartesian_, nonideal_,
		rotsets_, scminmap_, mingraph_,
		atc_, min_options_ );

}

void MMTMinPackingJob::optimize()
{

	best_assignment_ = core::pack::scmin::SidechainStateAssignmentOP( new core::pack::scmin::SidechainStateAssignment( rotsets_->nmoltenres() ) );

	core::pack::min_pack_optimize(
		pose(), sfxn(), task(), rotsets_, scminmap_,
		mingraph_, atc_, *min_options_, *best_assignment_ );

	// it'd be better if the minpacker tracked the total energy
	update_pose( pose() );
	best_assignment_->assign_energy( sfxn()( pose() ) );

}

/// @throws If the input pose / task do not match the saved best_state,
/// then this function throws an exception as it cannot apply the previously
/// saved state in a consistent way.
void MMTMinPackingJob::update_pose( core::pose::Pose & final_pose )
{
	// error checking
	std::string errormsg = "";
	if ( ! has_pose() ) { errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; no internal pose set\n"; }
	if ( ! has_sfxn() ) { errormsg += "MMTMinPackingJob::update_pose() --  Could not update pose; no sfxn present\n"; }
	if ( ! has_task() ) { errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; no task present\n"; }
	if ( ! atc_ ) { errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; atom-tree collection not yet initialized\n"; }
	if ( ! rotsets_ ) { errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; rotamer sets not yet initialized\n"; }
	if ( ! best_assignment_ ) { errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; no best state assignment present\n"; }

	if ( errormsg != "" ) { throw utility::excn::EXCN_Msg_Exception( errormsg ); }

	if ( rotsets_->nmoltenres() != best_assignment_->nmoltenres() ) {
		errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; nmoltenres disagreement between rotsets_ and best_assignment_\n";
	}
	if ( pose().total_residue() != final_pose.total_residue() ) {
		errormsg += "MMTMinPackingJob::update_pose() -- Could not update pose; internal pose and final_pose have different numbers of residues\n";
	}
	if ( errormsg != "" ) { throw utility::excn::EXCN_Msg_Exception( errormsg ); }

	// ok, if we made it this far, we're golden.
	core::pack::min_pack_place_opt_rotamers_on_pose( final_pose, sfxn(), rotsets_, atc_, *best_assignment_, 0.0 );

}

bool MMTMinPackingJob::best_assignment_exists() const
{
	return best_assignment_ != 0;
}


MMTMinPackingJob::SidechainStateAssignment const &
MMTMinPackingJob::get_best_assignment() const
{
	return *best_assignment_;
}

core::Real
MMTMinPackingJob::final_energy() const
{
	return best_assignment_->energy();
}


void MMTMinPackingJob::clean_up() {
	rotsets_.reset();
	scminmap_.reset();
	mingraph_.reset();
	atc_.reset();
	min_options_.reset();
	MMTPackingJob::clean_up();
}


}
}
