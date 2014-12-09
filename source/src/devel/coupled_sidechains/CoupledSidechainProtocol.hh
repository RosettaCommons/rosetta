// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   TemperedCoupled_Sidechains.hh
///
/// @brief allows low-resolution coupled_sidechains using simulated or parallel tempering
/// @author Oliver Lange

#ifndef INCLUDED_devel_coupled_sidechains_CoupledSidechainProtocol_hh
#define INCLUDED_devel_coupled_sidechains_CoupledSidechainProtocol_hh

// Unit Headers
#include <devel/coupled_sidechains/CoupledSidechainProtocol.fwd.hh>

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

#include <protocols/canonical_sampling/SimulatedTemperingObserver.hh>
#include <protocols/canonical_sampling/MultiTemperatureTrialCounter.hh>

#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

// Protocols Headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>


// Core Headers
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

#include <numeric/random/random.hh>

#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/vector0.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

//Auto Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1_bool.hh>

#include <core/scoring/ScoreFunction.hh>

namespace devel {
namespace coupled_sidechains {


/// @brief class for non-discrete side chain sampling using Dunbrack rotamer probabilities/distributions
class CoupledSidechainProtocol : public protocols::simple_moves::sidechain_moves::SidechainMover {

public:
//
	static void register_options();
	/// @brief default constructor
	CoupledSidechainProtocol();

	CoupledSidechainProtocol( CoupledSidechainProtocol const & ) {};

	/// @brief constructor with user provided rotamer library
	CoupledSidechainProtocol(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	~CoupledSidechainProtocol();

	void
	show_counters( std::ostream & out );


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

	void set_defaults();
	void init_from_options();
	void setup_objects();

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
	set_scorefunction( core::scoring::ScoreFunctionOP sfxn );

	core::scoring::ScoreFunction const& scorefunction() const {
		return *sfxn_;
	}

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle //non-zero if trajectory is restarted
	);

private:

	bool pass_metropolis( core::Real delta_energy , core::Real last_proposal_density_ratio );

	void
	perturb_chi(numeric::random::RandomGenerator Rand,
							core::Real max_deviation,
							utility::vector1<core::Real> & current_chi,
							utility::vector1<core::Real> & new_chi
	);

	core::Size output_count( core::Size ct ) {
		if ( ct % score_stride_ == 0 ) {
			return ct / score_stride_;
		}	else return 0;
	}

	void observe_rotamers( core::Size ct, std::string const& traj_file_tag );
	void update_rotamers( core::Size resid );

	//ek for fast sidechain sampling and internal mc trials
	utility::vector1< core::conformation::ResidueOP > current_;
	utility::vector1< core::conformation::ResidueOP > previous_;
	utility::vector1< core::conformation::ResidueOP > best_;
	utility::vector1< core::pack::dunbrack::ChiVector > chi_vectors_;
	utility::vector1< core::pack::dunbrack::RotVector > rot_vectors_;

	core::Size score_stride_;
	core::Size traj_stride_;
	core::Size rotamer_stride_;

	core::Real temperature_;
	core::Size ntrials_;
	core::Real best_energy_;
	core::Real current_energy_;
	core::scoring::ScoreFunctionOP sfxn_;

	core::pack::interaction_graph::SimpleInteractionGraphOP ig_;

	core::Real accepts_;
	core::Real current_ntrial_;
	core::Real score_pre_apply_;
	core::Real score_post_apply_;

	utility::io::ozstream rotamer_stream_;

	protocols::canonical_sampling::SimulatedTemperingObserverOP tempering_;
	protocols::canonical_sampling::MultiTemperatureTrialCounter counters_;

	static bool options_registered_;
}; //CoupledSidechainProtocol


} // moves
} // protocols

#endif
