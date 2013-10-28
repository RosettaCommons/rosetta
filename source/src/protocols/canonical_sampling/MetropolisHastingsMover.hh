// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/MetropolisHastingsMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_hh
#define INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_hh

// Unit Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.fwd.hh>

#ifdef WIN32
	#include <protocols/canonical_sampling/ThermodynamicMover.hh>
	#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#endif


namespace protocols {
namespace canonical_sampling {

///@details
class MetropolisHastingsMover : public protocols::moves::Mover {

public:

	MetropolisHastingsMover();

	MetropolisHastingsMover(
		MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	~MetropolisHastingsMover();

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
	bool
	reinitialize_for_each_job() const;

	virtual
	bool
	reinitialize_for_new_input() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	protocols::moves::MonteCarloCOP
	monte_carlo() const;

	void
	set_monte_carlo(
		protocols::moves::MonteCarloOP monte_carlo
	);

	void
	set_tempering(
		TemperatureControllerOP
	);

	TemperatureControllerCOP
	tempering() const;

	core::Size
	ntrials() const { return ntrials_; }

	void
	set_ntrials(
		core::Size ntrials
	);

	std::string const &
	output_name() const;

	void
	set_output_name(
		std::string const & output_name
	);

	std::string
	output_file_name(
		std::string const & suffix,
		bool cumulate_jobs = false,
		bool cumulate_replicas = false
	) const;

	bool
	finished() const;

	virtual
	ThermodynamicMoverOP
	random_mover();

	virtual
	void
	add_mover(
		ThermodynamicMoverOP mover,
		core::Real weight,
		utility::tag::TagCOP const& subtag
	);

	virtual void
	add_mover(
		ThermodynamicMoverOP mover,
		core::Real weight
	);

	void
	add_backrub_mover(
		core::Real weight
	);

	void
	add_small_mover(
		core::Real weight
	);

	void
	add_shear_mover(
		core::Real weight
	);

	void
	add_sidechain_mover(
		core::Real weight,
		core::Real prob_uniform,
		core::Real prob_withinrot,
		bool preserve_cbeta
	);

	void
	add_sidechain_mc_mover(
		core::Real weight,
		core::Real prob_uniform,
		core::Real prob_withinrot,
		bool preserve_cbeta,
		core::Size ntrials
	);

	void
	add_observer(
		ThermodynamicObserverOP observer
	);

	ThermodynamicMover const&
	last_move() const;

	bool
	last_accepted() const {
		return last_accepted_;
	}

protected:

	TemperatureControllerOP const&
	tempering() { return tempering_; }

	ThermodynamicMoverOP
	mover_by_index(numeric::Size idx) const
	{
		return movers_[idx];
	}

	void wind_down_simulation( core::pose::Pose& pose);

	//initialize all movers
	//return cycle number if a restart, 0 otherwise
	core::Size prepare_simulation( core::pose::Pose& pose);

	void set_last_accepted( bool setting ) {
		last_accepted_ = setting;
	}

	void set_last_move( ThermodynamicMoverOP setting );

	typedef utility::vector1< ThermodynamicObserverOP > ObserverList;
	ObserverList const& observers() { return observers_; }

	protocols::moves::MonteCarlo& nonconst_monte_carlo();
private:

	///configurables...
	protocols::moves::MonteCarloOP monte_carlo_;
	core::Size ntrials_;
	std::string output_name_;
	utility::vector1< ThermodynamicMoverOP > movers_;
	utility::vector1< ThermodynamicObserverOP > observers_;
	TemperatureControllerOP tempering_;

	//helper
	numeric::random::WeightedSampler weighted_sampler_;

	///some status is necessary for the observers
	ThermodynamicMoverOP last_move_;
	bool last_accepted_;
	core::Size trial_;

	//internal book keeping
	bool output_name_from_job_distributor_;

}; //end MetropolisHastingsMover

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_HH
