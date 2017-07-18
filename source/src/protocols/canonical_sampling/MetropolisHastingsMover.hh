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

#ifndef INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_hh
#define INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_hh

// Unit Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <vector>
#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.fwd.hh>

#ifdef WIN32
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#endif


namespace protocols {
namespace canonical_sampling {

/// @brief Manage the main loop of a canonical Monte Carlo simulation.
///
/// @details To make the simulation flexible, most aspects of the algorithm
/// have been delegated to other classes.  Use the add_mover() methods to
/// control which moves are used during the simulation.  Use the
/// set_tempering() method to control how the temperature changes during the
/// simulation.  This can be used to setup simulated annealing or parallel
/// tempering runs.  Management of the score function is delegated to the
/// underlying MonteCarlo object, so use set_monte_carlo() to specify a score
/// function.  Use add_observer() to keep track of statistics and to record the
/// trajectory.

class MetropolisHastingsMover : public protocols::moves::Mover {

public:

	/// @brief Default constructor.
	MetropolisHastingsMover();

	/// @brief Copy constructor.
	MetropolisHastingsMover(
		MetropolisHastingsMover const & metropolis_hastings_mover
	);

	/// @brief Destructor.

	~MetropolisHastingsMover() override;

	/// @brief Run the Metropolis-Hastings simulation.

	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Return the name of this mover.

	// XRW TEMP  std::string
	// XRW TEMP  get_name() const override;

	/// @brief Return a copy of this mover.
	protocols::moves::MoverOP
	clone() const override;

	/// @brief Return a newly instantiated mover.

	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief Return false.  This mover does not need to be reinitialized for
	/// each job.

	bool
	reinitialize_for_each_job() const override;

	/// @brief Return false.  This mover does not need to be reinitialized for
	/// new input.

	bool
	reinitialize_for_new_input() const override;

	/// @brief Use a RosettaScripts tag to configure this mover.

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	/// @brief Return the MonteCarlo object being used by this simulation.
	protocols::moves::MonteCarloCOP
	monte_carlo() const;

	/// @brief Provide a MonteCarlo object to use for this simulation.
	void
	set_monte_carlo(
		protocols::moves::MonteCarloOP monte_carlo
	);

	/// @brief Return the TemperatureController being used by this simulation.
	TemperatureControllerCOP
	tempering() const;

	/// @brief Provide a TemperatureController to use for this simulation.
	void
	set_tempering(
		TemperatureControllerOP
	);

	/// @brief Return the number of iterations used by this simulation.
	core::Size
	ntrials() const { return ntrials_; }

	/// @brief Set the number of iterations to use for this simulation.
	void
	set_ntrials(
		core::Size ntrials
	);

	/// @brief Return the iteration currently being processed by the simulation.
	core::Size
	current_trial() const { return current_trial_; }

	/// @brief Return the file name used by some of the observers to output data.
	std::string const &
	output_name() const;

	/// @brief Set the file name used by some of the observers to output data.
	void
	set_output_name(
		std::string const & output_name
	);

	/// @brief Return a file name that is consistent with the given options.
	/// @details If @a cumulate_jobs is true, the same filename will be returned
	/// for different jobs, so that the jobs all get cumulated in the same file.
	/// Likewise, if @a cumulate_replicas is true, the same filename will be
	/// returned for all replicas.  If either of these options are set, MPI must
	/// be enabled.
	std::string
	output_file_name(
		std::string const & suffix,
		bool cumulate_jobs = false,
		bool cumulate_replicas = false
	) const;

	/// @brief Return true if the simulation has been completed.
	// Undefined, commenting out to fix PyRosetta build  bool finished() const;

	/// @brief Return a randomly chosen mover to use in the next iteration.
	virtual
	ThermodynamicMoverOP
	random_mover() const;

	/// @brief Add the given mover to the simulation.
	virtual
	void
	add_mover(
		ThermodynamicMoverOP mover,
		core::Real weight,
		utility::tag::TagCOP const& subtag
	);

	/// @brief Add the given mover to the simulation.
	virtual void
	add_mover(
		ThermodynamicMoverOP mover,
		core::Real weight
	);

	/// @brief Convenience method to add a backrub move to the simulation.
	void
	add_backrub_mover(
		core::Real weight
	);

	/// @brief Convenience method to add a kinematic closure move to the
	/// simulation.
	void
	add_kic_mover(
		core::Real weight,
		protocols::loops::Loop const & loop
	);

	/// @brief Convenience method to add a small move to the simulation.
	void
	add_small_mover(
		core::Real weight
	);

	/// @brief Convenience method to add a shear move to the simulation.
	void
	add_shear_mover(
		core::Real weight
	);

	/// @brief Convenience method to add a sidechain move to the simulation.
	void
	add_sidechain_mover(
		core::Real weight,
		core::Real prob_uniform,
		core::Real prob_withinrot,
		bool preserve_cbeta
	);

	/// @brief Convenience method to add a Monte Carlo sidechain move to the
	/// simulation.  This move uses an internal Monte Carlo loop to generate a
	/// whole new set of sidechain conformations.
	void
	add_sidechain_mc_mover(
		core::Real weight,
		core::Real prob_uniform,
		core::Real prob_withinrot,
		bool preserve_cbeta,
		core::Size ntrials
	);

	/// @brief Add the given observer to this simulation.
	void
	add_observer(
		ThermodynamicObserverOP observer
	);

	/// @brief Return the most recently used ThermodynamicMover.
	ThermodynamicMover const&
	last_move() const;

	/// @brief Return true if the last attempted move was accepted.
	bool
	last_accepted() const {
		return last_accepted_;
	}

	std::string
	get_last_checkpoint() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	static
	std::string
	mover_or_add_group_name();

	static
	std::string
	add_ct_name( std::string );


protected:

	/// @brief Protected non-const access to the TemperatureController.
	TemperatureControllerOP const&
	tempering() { return tempering_; }

	/// @brief Return the mover that was added at the given index.
	ThermodynamicMoverOP
	mover_by_index(numeric::Size idx) const
	{
		return movers_[idx];
	}

	/// @brief write checkpoint snapshots for restarting
	void
	write_checkpoint( core::pose::Pose const & pose );

	/// @brief get checkpoint_id for restarting
	bool
	get_checkpoints();

	/// @brief Finalize all the movers and observers used in this simulation, and
	/// write some debrief statistics to the tracer.
	void wind_down_simulation( core::pose::Pose& pose);

	/// @brief Initialize all the movers and observers to be used in this
	/// simulation.
	core::Size prepare_simulation( core::pose::Pose& pose);

	/// @brief Indicate whether or not the last attempted move was accepted.
	void set_last_accepted( bool setting ) {
		last_accepted_ = setting;
	}

	/// @brief Indicate what type of move was last attempted.
	void set_last_move( ThermodynamicMoverOP setting );

	/// @brief Return a list of all observers attached to this simulation.
	utility::vector1< ThermodynamicObserverOP > const& observers() {
		return observers_;
	}

	/// @brief Protected non-const access to the MonteCarlo object.
	protocols::moves::MonteCarlo& nonconst_monte_carlo();
private:

	// Configurable
	protocols::moves::MonteCarloOP monte_carlo_;
	core::Size ntrials_;
	core::Size current_trial_;
	std::string output_name_;
	utility::vector1< ThermodynamicMoverOP > movers_;
	utility::vector1< ThermodynamicObserverOP > observers_;
	TemperatureControllerOP tempering_;

	// Helper
	numeric::random::WeightedSampler weighted_sampler_;

	// Some status is necessary for the observers
	ThermodynamicMoverOP last_move_;
	bool last_accepted_;

	// Internal book keeping
	bool output_name_from_job_distributor_;

	core::Size checkpoint_count_;
	std::vector< std::string > checkpoint_ids_;

}; //end MetropolisHastingsMover

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_HH
