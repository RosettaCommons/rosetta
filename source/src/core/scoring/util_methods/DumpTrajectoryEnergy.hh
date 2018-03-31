// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/util_methods/DumpTrajectoryEnergy.hh
/// @brief Headers for an EnergyMethod that dumps trajectories to file.
/// @details Dumps trajectories of the minimizer and packer to file when the dump_trajectory
/// ScoreType is enable. Output may be controlled through the dump_trajectory:* flags.
/// @author Brian Coventry (bcov@uw.edu). - dump minimizer


// Goals for this class
// 1. Enabling the dump_trajectory ScoreType should be sufficient to dump trajectories.
// 2. Other classes should not need to cooperate. This class must be robust.
// 3. Different activities get different filenames (minimizing vs. packing)
// 4. Each call to an activity produces exactly 1 file (i.e. each call to MinMover produces 1 file)
// 5. For a given trajectory, filenames should be alphabetically organized chronologically
// 6. Low overhead. Write the poses out as they come in. Don't store them.
// 7. The activity being dumped should give the same results as when not dumped.
// 8. All configuration for this class should be done through the dump_trajectory::* options


#ifndef INCLUDED_core_scoring_util_methods_DumpTrajectoryEnergy_hh
#define INCLUDED_core_scoring_util_methods_DumpTrajectoryEnergy_hh

// Unit headers
#include <core/scoring/util_methods/DumpTrajectoryEnergy.fwd.hh>

// Package headers
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace util_methods {

/// @brief DumpTrajectoryEnergy, a util_method that allows one to dump
/// rosetta trajectories through the use of the dump_trajectory ScoreType.
class DumpTrajectoryEnergy : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;

public:

	/// @brief Options constructor.
	///
	DumpTrajectoryEnergy( core::scoring::methods::EnergyMethodOptions const &options );

	/// @brief Copy constructor.
	///
	DumpTrajectoryEnergy( DumpTrajectoryEnergy const &src );

	/// @brief Default destructor.
	///
	~DumpTrajectoryEnergy() override;

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief DumpTrajectoryEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief DumpTrajectoryEnergy is version 1.0 right now.
	///
	core::Size version() const override;


	/********************************************************************
	Minimizing
	********************************************************************/

	/// @brief This is where the state_ moves to MINIMIZING.
	///
	void setup_for_minimizing( pose::Pose & , ScoreFunction const & , kinematics::MinimizerMapBase const &) const override;

	/// @brief This is where the minimizer frames are dumped.
	/// @details Assumption: All minimizers must evaluate derivatives repeatedly while minimizing.
	/// Although these won't be equally spaced, dumping them is the best we can do without modifying the minimizers.
	void setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const override;


	/********************************************************************
	Packing
	********************************************************************/

	/// @brief This is where the state_ moves to PACKING.
	///
	void set_up_residuearrayannealableenergy_for_packing ( core::pose::Pose const &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, core::scoring::ScoreFunction const &sfxn) override;

	/// @brief Should be possible to dump the packer from here. Not implemented yet though
	///
	core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, core::Size const substitution_position = 0 ) const override;


private:

	/******************
	Private constants:
	******************/

	enum DumpState {
		IDLE,
		MINIMIZING,  // has no ending state, moving to IDLE is always allowable
		PACKING   // not implemented yet
	};


	/******************
	Private functions:
	******************/

	/// @brief Ensure that the previous activity has finished before starting a new one
	/// @details This function is to be called when starting a new activity to catch invalid states.
	/// This throws an exception when an invalid state is detected.
	void try_move_to_idle() const;


	/// @brief Perform initial setup related to an activity
	/// @details This ensure the switch is valid based on the current state_.
	void start_new_activity( DumpState new_state ) const;


	void dump_pose( core::pose::Pose const & pose ) const;

	/******************
	Private variables:
	******************/

	/// @brief The activity we are currently dumping
	/// @details One may set this manually. But moves to new activities should be done through start_new_activity()
	mutable DumpState state_;

	/// @brief The filename where frames are being dumped to
	/// @details
	mutable std::string dump_filename_;

	/// @brief Number of dumped frames so far
	/// @details This is used to give the models numbers in the pdb
	mutable Size dumped_frames_;


	std::string dump_prefix_;
	bool dump_as_gz_;

};

} // util_methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
