// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MonteCarloInterface.hh
/// @brief A MonteCarlo object for optimizing the interface dG as defined using InterfaceAnalyzer. The dG and Total energy can be weighted.  This is so that the interface energy itself can be optimized through a protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_MonteCarloInterface_hh
#define INCLUDED_protocols_simple_moves_MonteCarloInterface_hh

#include <protocols/simple_moves/MonteCarloInterface.fwd.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace simple_moves {

///@brief A MonteCarlo object for optimizing the interface dG as defined using InterfaceAnalyzer. The dG and Total energy can be weighted.  This is so that the interface energy itself can be optimized through a protocol.
class MonteCarloInterface : public protocols::moves::MonteCarlo {

public:

	/// @brief Constructs a useable MonteCarlo object
	///
	/// mc = MonteCarlo( init_pose , scorefxn , temp )
	///
	/// Pose           init_pose   /manipulated during the simulation
	/// ScoreFunction  scorefxn    /evaluates pose scores
	/// Real (float)   temp        /used in the Metropolis Criterion
	MonteCarloInterface(
		Pose const & init_pose, // PoseCOP init_pose,
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		Real const temperature,
		std::string interface
	);

	/// @brief Constructor without Pose -- call reset(pose) before first use
	MonteCarloInterface(
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		Real const temperature,
		std::string interface
	);

	MonteCarloInterface(MonteCarloInterface const & src);

	~MonteCarloInterface() override;

	protocols::moves::MonteCarloOP
	clone() override;

	/// @brief Applies the Metropolis Criterion on pose based on
	/// the Interface dG and Total Score according to Set weights.
	///
	/// Interface dG is calculated using the InterfaceAnalyzer Application.
	///
	/// pose. This method evaluates the change in score, compares
	/// the trial pose to the last accepted pose, and updates the
	/// pose structure and simulation statistics appropriately
	///
	/// example(s):
	///     mc.boltzmann( pose )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	bool
	boltzmann(
		Pose & pose,//PoseOP pose,
		std::string const & move_type = "unk",
		core::Real const proposal_density_ratio = 1,
		core::Real const inner_score_delta_over_temperature = 0
	) override;

	/// @brief Sets lowest score pose and last accepted pose to
	/// the score of  <pose>
	/// @note (does not reset counters)
	///
	/// example(s):
	///     mc.reset(pose)
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	///     MonteCarlo.lowest_scored_pose
	void
	reset( Pose const & pose ) override;

public:

	///Set the interface that we will be using to calculate interface energy.
	void
	set_interface( std::string const & interface );

	///@brief Should we pack the interface while separated?  Default TRUE.
	void
	set_pack_interface_separated( bool pack_separated );

	///@brief Set the weight on the dG component of the energy (Default = 1.0)
	void
	set_dG_weight( core::Real dG_weight );

	///@brief Set the weight on the total component of the enegy (Default = 0.0)
	void
	set_total_weight( core::Real total_weight );

	void
	score_function(ScoreFunction const & scorefxn ) override;

	///@brief Calculate the score given the pose and the weights.
	core::Real
	calculate_score( Pose & pose );

	///@brief Should we repack only the interface residues of a cloned pose on calculating the interface energy?
	/// Default TRUE.
	void
	set_repack_separated( bool repack_separated );

private:

	void
	init_interface_analyzer();

	void
	read_cmd_line_options();

private:

	std::string interface_definition_ = "";
	protocols::analysis::InterfaceAnalyzerMoverOP analyzer_;

	core::Real dG_weight_ = 1.0;
	core::Real total_weight_ = 0.0;

	core::Real last_dG_ = 0.0;
	core::Real last_score_ = 0.0;

	bool repack_separated_ = true;


};


} //protocols
} //simple_moves



#endif //INCLUDED_protocols_simple_moves_MonteCarloInterface_hh





