// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMCMover.hh
/// @brief definition of SidechainMCMover class and functions
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_SidechainMCMover_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_SidechainMCMover_hh

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.fwd.hh>

// Protocols Headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>


// Core Headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>


#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>


#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#include <core/conformation/Residue.hh>
#endif


namespace protocols {
namespace simple_moves {
namespace sidechain_moves {


/// @brief class for non-discrete side chain sampling using Dunbrack rotamer probabilities/distributions
class SidechainMCMover : public protocols::simple_moves::sidechain_moves::SidechainMover {

public:
//
	/// @brief default constructor
	SidechainMCMover();

	/// @brief constructor with user provided rotamer library
	SidechainMCMover(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	~SidechainMCMover();

	void
	show_counters( std::ostream & out );

	void
	setup( core::scoring::ScoreFunctionCOP sfxn );

	//parser stuff
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	/// @brief apply a sidechain move to a Pose object
	void
	apply(
		core::pose::Pose & pose
	);

	virtual std::string get_name() const;

	void
	set_ntrials( core::Size ntrial ){
		ntrials_ = ntrial;
	};

	core::Size
	ntrials(){
		return 	ntrials_;
	}

	void
	set_temperature( core::Real temp ){
		temperature_ = temp;
	}

	core::Real
	temperature(){
		return temperature_;
	}

	void
	set_inherit_scorefxn_temperature( bool inherit_scorefxn_temperature )
	{
		inherit_scorefxn_temperature_ = inherit_scorefxn_temperature;
	}

	bool
	inherit_scorefxn_temperature() const
	{
		return inherit_scorefxn_temperature_;
	}

	void
	set_scorefunction( core::scoring::ScoreFunction const & sfxn ){
		sfxn_ = sfxn.clone();
	}

	core::scoring::ScoreFunctionCOP
	scorefunction(){
		return sfxn_;
	}

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual
	core::Real
	last_proposal_density_ratio()
	{
		return 1;
	}

	virtual
	bool
	is_multi_trial()
	{
		return true;
	}

	virtual
	core::Real
	last_inner_score_delta_over_temperature()
	{
		return (score_post_apply_-score_pre_apply_)/temperature_;
	}

	virtual
	protocols::canonical_sampling::MetropolisHastingsMoverAP
	metropolis_hastings_mover()
	{
		return metropolis_hastings_mover_;
	}

	virtual
	void
	set_metropolis_hastings_mover(
		protocols::canonical_sampling::MetropolisHastingsMoverAP metropolis_hastings_mover
	)
	{
		metropolis_hastings_mover_ = metropolis_hastings_mover;
	}


private:

	bool pass_metropolis( core::Real delta_energy , core::Real last_proposal_density_ratio );

	void
	perturb_chi(numeric::random::RandomGenerator Rand,
							core::Real max_deviation,
							utility::vector1<core::Real> & current_chi,
							utility::vector1<core::Real> & new_chi
	);

	//ek for fast sidechain sampling and internal mc trials
	utility::vector1< core::conformation::ResidueOP > current_;
	utility::vector1< core::conformation::ResidueOP > previous_;
	utility::vector1< core::conformation::ResidueOP > best_;
	core::Real temperature_;
	core::Size ntrials_;
	core::Real best_energy_;
	core::Real current_energy_;
	core::scoring::ScoreFunctionOP sfxn_;
	bool inherit_scorefxn_temperature_;
	core::pack::interaction_graph::SimpleInteractionGraphOP ig_;
	core::Real accepts_;
	core::Real current_ntrial_;
	core::Real score_pre_apply_;
	core::Real score_post_apply_;
	protocols::canonical_sampling::MetropolisHastingsMoverAP metropolis_hastings_mover_;

}; //SidechainMCMover


} // sidechain_moves
} // simple_moves
} // protocols

#endif
