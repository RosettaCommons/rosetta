// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/canonical_sampling/MetropolisHastingsMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_SidechainMetropolisHastingsMover_hh
#define INCLUDED_protocols_canonical_sampling_SidechainMetropolisHastingsMover_hh

// Unit Headers
#include <protocols/canonical_sampling/SidechainMetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Run a sidechain-only canonical Monte Carlo simulation.
///
/// @details The Monte Carlo algorithm present in apply() has been optimized
/// for the case where only the sidechains can moves.  This makes it possible
/// speed up score function evaluation by either precalculating and/or caching
/// residue pair energies.  In this specific case, however, I'm not exactly
/// sure how the algorithm is doing its optimization.
///
/// @warning Although this class inherits from MetropolisHastingsMover, it
/// doesn't support all of its parent's interface.  In particular, since the
/// algorithm is customized for a particular sidechain move, movers added via
/// add_mover() or its related methods are ignored.  However, observers added
/// via add_observer() can still be used to report on the state of the
/// simulation.

class SidechainMetropolisHastingsMover : public protocols::canonical_sampling::MetropolisHastingsMover {

public:
	typedef MetropolisHastingsMover Parent;

	/// @brief Default constructor.
	SidechainMetropolisHastingsMover();

	/// @brief Constructor with stride parameter.
	SidechainMetropolisHastingsMover( core::Size stride );

	/// @brief Copy constructor.
	SidechainMetropolisHastingsMover(
		SidechainMetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief Default destructor
	virtual
	~SidechainMetropolisHastingsMover();

	/// @brief Run the sidechain-only simulation.
	virtual
	void
	apply( core::pose::Pose & pose );

	virtual
	std::string
	get_name() const;

	protocols::moves::MoverOP
	clone() const;

	virtual
	protocols::moves::MoverOP
	fresh_instance() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief Return true if a move should be accepted, given @a delta_energy
	/// and @a proposal_density_ratio.
	bool pass_metropolis( core::Real delta_energy , core::Real last_proposal_density_ratio ) const;

	/// @brief Return non-zero if the observers should be invoked on this
	/// iteration.
	/// @see set_stride()
	core::Size output_count( core::Size ct ) const;

	/// @brief Set the frequency with which the observers should be invoked.
	/// @see output_count()
	void set_stride( core::Size setting ) { stride_ = setting; };

private:
	core::Size stride_;
}; //end SidechainMetropolisHastingsMover

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_SidechainMetropolisHastingsMover_HH
