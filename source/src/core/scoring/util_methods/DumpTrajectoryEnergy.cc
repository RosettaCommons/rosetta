// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/util_methods/DumpTrajectoryEnergy.cc
/// @brief An EnergyMethod that dumps trajectories to file.
/// @details Dumps trajectories of the minimizer and packer to file when the dump_trajectory
/// ScoreType is enable. Output may be controlled through the dump_trajectory:* flags.
/// @author Brian Coventry (bcov@uw.edu). - dump minimizer

// Unit headers
#include <core/scoring/util_methods/DumpTrajectoryEnergy.hh>
#include <core/scoring/util_methods/DumpTrajectoryEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/io/pdb/pdb_writer.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/dump_trajectory.OptionKeys.gen.hh>

// Utility headers
#include <utility/sys_util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace util_methods {

static basic::Tracer TR("core.scoring.util_methods.DumpTrajectoryEnergy");

/// @brief This must return a fresh instance of the DumpTrajectoryEnergy class, never an instance already in use.
///
core::scoring::methods::EnergyMethodOP
DumpTrajectoryEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & options ) const
{
	return core::scoring::methods::EnergyMethodOP( new DumpTrajectoryEnergy( options ) );
}

/// @brief Defines the score types that this energy method calculates.
///
ScoreTypes
DumpTrajectoryEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( dump_trajectory );
	return sts;
}

/// @brief Options constructor.
///
DumpTrajectoryEnergy::DumpTrajectoryEnergy ( core::scoring::methods::EnergyMethodOptions const & ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new DumpTrajectoryEnergyCreator ) ),
	parent2( ),
	state_( IDLE ),
	dump_filename_( "" ),
	dumped_frames_( 0 ),
	dump_prefix_( basic::options::option[ basic::options::OptionKeys::dump_trajectory::prefix ]() ),
	dump_as_gz_( basic::options::option[ basic::options::OptionKeys::dump_trajectory::gz ]() )
{}

/// @brief Copy constructor.
/// @details Must move to IDLE before we duplicate.
DumpTrajectoryEnergy::DumpTrajectoryEnergy( DumpTrajectoryEnergy const &src ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new DumpTrajectoryEnergyCreator ) ),
	parent2( src ) {

	// Move to IDLE first. This class abuses mutable so we must do this first.
	src.try_move_to_idle();
	assert( src.state_ == IDLE );

	state_ = src.state_;
	dump_filename_ = src.dump_filename_;
	dumped_frames_ = src.dumped_frames_;
}

/// @brief Default destructor.
///
DumpTrajectoryEnergy::~DumpTrajectoryEnergy() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP DumpTrajectoryEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new DumpTrajectoryEnergy(*this) );
}

/// @brief DumpTrajectoryEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void DumpTrajectoryEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief DumpTrajectoryEnergy is version 1.0 right now.
///
core::Size DumpTrajectoryEnergy::version() const
{
	return 1; // Initial versioning
}


/********************************************************************
Minimizing
********************************************************************/

/// @brief This is where the state_ moves to MINIMIZING.
///
void
DumpTrajectoryEnergy::setup_for_minimizing( pose::Pose & , ScoreFunction const & , kinematics::MinimizerMapBase const &) const {
	start_new_activity( MINIMIZING );

}

/// @brief This is where the minimizer frames are dumped.
/// @details Assumption: All minimizers must evaluate derivatives repeatedly while minimizing.
/// Although these won't be equally spaced, dumping them is the best we can do without modifying the minimizers.
void
DumpTrajectoryEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {
	dump_pose( pose );
}


/********************************************************************
Packing
********************************************************************/


/// @brief This is where the state_ moves to PACKING.
///
void
DumpTrajectoryEnergy::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose &,
	core::pack::rotamer_set::RotamerSets const &,
	core::scoring::ScoreFunction const &/*sfxn*/
) {
	start_new_activity( PACKING );
}

/// @brief Should be possible to dump the packer from here. Not implemented yet though
///
core::Real
DumpTrajectoryEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &,
	core::Size const /*substitution_position = 0*/ ) const {
	return 0;
}


/********************************************************************
Dumping
********************************************************************/

/// @brief Ensure that the previous activity has finished before starting a new one
/// @details This function is to be called when starting a new activity to catch invalid states.
/// This throws an exception when an invalid state is detected.
void
DumpTrajectoryEnergy::try_move_to_idle() const {
	switch ( state_ ) {
	case IDLE : {
		break;
	}
	case MINIMIZING : {
		// We don't get an ending signal on the minimizer. So we can always switch.
		break;
	}
	case PACKING : {
		// We don't get an ending signal on the packer. So we can always switch.
		break;
	}
	}
	// no one threw an exception. So switching to IDLE must be ok.
	state_ = IDLE;
}


/// @brief Perform initial setup related to an activity
/// @details This ensure the switch is valid based on the current state_.
void
DumpTrajectoryEnergy::start_new_activity( DumpState new_state ) const {
	try_move_to_idle();

	state_ = new_state;

	std::stringstream filename;
	filename << dump_prefix_;
	filename << "_" << utility::timestamp_millis_short();

	switch ( state_ ) {
	case IDLE : {
		break;
	}
	case MINIMIZING : {
		if ( TR.visible() ) {
			TR << "Dumping minimization trajectory to: ";
		}
		filename << "min";
		break;
	}
	case PACKING : {
		if ( TR.visible() ) {
			TR << "Dumping packer trajectory to: ";
		}
		filename << "pack";
	}
	}

	filename << ".pdb";
	if ( dump_as_gz_ ) {
		filename << ".gz";
	}
	dump_filename_ = filename.str();
	dumped_frames_ = 0;

	if ( state_ != IDLE ) {
		TR << dump_filename_ << std::endl;
	}
}


void
DumpTrajectoryEnergy::dump_pose( core::pose::Pose const & pose ) const {
	if ( state_ == IDLE ) {
		utility_exit_with_message("Tried to dump pose from IDLE state!!! Email bcov@uw.edu");
	}

	dumped_frames_++;
	std::stringstream model_tag;
	model_tag << std::setfill('0') << std::setw( 6 ) << dumped_frames_;

	core::io::StructFileRepOptionsOP options( new core::io::StructFileRepOptions() );
	options->set_output_pose_energies_table( false );         // this causes problems during the minimizer

	io::pdb::add_to_multimodel_pdb( pose, dump_filename_, model_tag.str(), options );
}




} // util_methods
} // scoring
} // core
