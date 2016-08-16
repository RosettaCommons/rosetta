// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   TemperedCoupled_Sidechains.hh
///
/// @brief allows low-resolution coupled_sidechains using simulated or parallel tempering
/// @author Oliver Lange

#ifndef INCLUDED_devel_coupled_sidechains_CoupledSidechainProtocol2_hh
#define INCLUDED_devel_coupled_sidechains_CoupledSidechainProtocol2_hh

// Unit Headers
#include <devel/coupled_sidechains/CoupledSidechainProtocol.fwd.hh>

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

#include <protocols/canonical_sampling/SimulatedTempering.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

// Protocols Headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>


// Core Headers
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

//Auto Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1_bool.hh>

#include <core/scoring/ScoreFunction.hh>

namespace devel {
namespace coupled_sidechains {


/// @brief class for non-discrete side chain sampling using Dunbrack rotamer probabilities/distributions
class CoupledSidechainProtocol : public protocols::moves::Mover {

public:
	//
	static void register_options();
	/// @brief default constructor
	CoupledSidechainProtocol();

	CoupledSidechainProtocol( CoupledSidechainProtocol const & ) : protocols::moves::Mover() {}

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

private:
	core::Real prob_jump_rot_;
	core::Real prob_withinrot_;
	core::Real prob_pert_chi_;
	core::Real pert_magnitude_;
	core::Size ntrials_;
	core::Size stride_;
	protocols::canonical_sampling::MetropolisHastingsMoverOP sampler_;
	static bool options_registered_;
}; //CoupledSidechainProtocol


} // moves
} // protocols

#endif
