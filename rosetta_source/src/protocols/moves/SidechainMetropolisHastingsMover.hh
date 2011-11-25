// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/MetropolisHastingsMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_SidechainMetropolisHastingsMover_hh
#define INCLUDED_protocols_moves_SidechainMetropolisHastingsMover_hh

// Unit Headers
#include <protocols/moves/SidechainMetropolisHastingsMover.fwd.hh>
#include <protocols/moves/MetropolisHastingsMover.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/ThermodynamicMover.hh>
#include <protocols/moves/ThermodynamicObserver.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace moves {

///@details
class SidechainMetropolisHastingsMover : public MetropolisHastingsMover {

public:
	typedef MetropolisHastingsMover Parent;

	SidechainMetropolisHastingsMover();

	SidechainMetropolisHastingsMover( core::Size stride );

	SidechainMetropolisHastingsMover(
		SidechainMetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	~SidechainMetropolisHastingsMover();

	virtual
	void
	apply( core::pose::Pose & pose );

	virtual
	std::string
	get_name() const;

	MoverOP
	clone() const;

	virtual
	MoverOP
	fresh_instance() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		DataMap & data,
		protocols::filters::Filters_map const & filters,
		Movers_map const & movers,
		core::pose::Pose const & pose
	);

	bool pass_metropolis( core::Real delta_energy , core::Real last_proposal_density_ratio ) const;

	core::Size output_count( core::Size ct ) const;

	void set_stride( core::Size setting ) { stride_ = setting; };

private:
	core::Size stride_;
}; //end SidechainMetropolisHastingsMover

} //namespace moves
} //namespace protocols

#endif //INCLUDED_protocols_moves_SidechainMetropolisHastingsMover_HH
