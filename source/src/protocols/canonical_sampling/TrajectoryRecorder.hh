// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/TrajectoryRecorder.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_TrajectoryRecorder_hh
#define INCLUDED_protocols_canonical_sampling_TrajectoryRecorder_hh

// Project forward headers
#include <protocols/canonical_sampling/TrajectoryRecorder.fwd.hh>

// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace canonical_sampling {

/// @brief Base class for recording a simulation trajectory.
///
/// @details This class seems a little too geared towards file IO, which will 
/// make it awkward to create the database trajectory subclass that I want.  
/// But I'll get it to work one way or another.

class TrajectoryRecorder : public ThermodynamicObserver {
public:

	/// @brief Associate relevant options with the TemperedDocking class.
	static void register_options();

	/// @brief Default constructor.
	TrajectoryRecorder();

	/// @brief Destructor.
	~TrajectoryRecorder();

	/// @brief Copy constructor.
	TrajectoryRecorder( TrajectoryRecorder const & );

private:

	/// @brief Assignment not allowed.
	TrajectoryRecorder&
	operator=( TrajectoryRecorder const & );

public:

	/// @brief Return the name of this mover.
	virtual	std::string	get_name() const;

	/// @brief Configure this mover from a RosettaScripts tag.
	virtual	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief Return the file name for the trajectory.
	std::string const &	file_name() const {
		return file_name_;
	}

	/// @brief Set the file name for the trajectory.
	void file_name( std::string const & file_name )	{
		file_name_ = file_name;
	}

	/// @brief Return the number of models that have been saved so far.
	core::Size model_count() const {
		return model_count_;
	}

	/// @brief Return the number of iterations that have occurred so far.
	core::Size step_count() const {
		return step_count_;
	}

	/// @brief Return how often models should be written to the trajectory.
	core::Size stride() const {
		return stride_;
	}
	
	/// @brief Set how often models should be written to the trajectory.
	/// @details This option can also be specified on the command line using the 
	/// <tt> -trajectory:stride </tt> flag.
	void stride( core::Size stride ) {
		stride_ = stride;
	}

	/// @brief Return the number of poses that can be cached.
	core::Size cache_limit() {
		return cache_limit_;
	}

	/// @brief Specify the maximum number of poses that can be cached.
	/// @details This option can also be specified on the command line using the 
	/// <tt> -trajectory:cache_limit </tt> flag.  Note that some recorders don't 
	/// use a cache at all, and will therefore ignore this option.  
	void cache_limit( core::Size limit ) {
		cache_limit_ = limit;
	}

	/// @brief Return true if poses from different jobs will be written to the 
	/// same trajectory file.
	/// @details I suspect this is only meant to be used in the context of jd2.  
	/// This option can only be set from the command line using the <tt>
	/// -trajectory:cumulate_jobs </tt> flag.
	bool cumulate_jobs() const {
		return cumulate_jobs_;
	}

	/// @brief Return true if poses from different replicas will be written to 
	/// the same trajectory file.
	/// @details I suspect this is only meant to be used in the context of jd2.
	/// This option can only be set from the command line using the <tt>
	/// -trajectory:cumulate_replicas </tt> flag.
	bool cumulate_replicas() const {
		return cumulate_replicas_;
	}

	/// @brief Callback executed whenever the simulation is initialized or reset.
	virtual void reset(
		protocols::moves::MonteCarlo const& mc,
		MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	/// @brief Save the given pose to the trajectory.
	/// @details I suspect this method is provided so that the trajectory classes 
	/// can be used outside the canonical_sampling namespace.  That explains why 
	/// it's public, when it doesn't do anything but provide a simpler way for 
	/// children to respond to observe_after_metropolis().  It also explains why 
	/// it has so many different signatures and default arguments.  Note that 
	/// this function expects a <tt>MetropolisHastingsMoverCAP</tt> instead of 
	/// the more usual <tt>MetropolisHastingsMover const &</tt>.  This is because 
	/// only pointers are compatible with null default arguments.
	/// 
	/// That said, it seems to me that the @a metropolis_hastings_mover argument 
	/// is basically useless because you can't ever be sure it'll actually be 
	/// provided.  It's also surprising that observe_after_metropolis() is still 
	/// virtual, because I would think that this method should supersede it.
	void update_after_boltzmann(
		core::pose::Pose const & pose,
		MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	/// @copydoc update_after_boltzmann(core::pose::Pose const &, MetropolisHastingsMoverCAP)
	void update_after_boltzmann(
		protocols::moves::MonteCarlo const& mc
	);

	/// @copydoc ThermodynamicObserver::apply
	virtual void apply( core::pose::Pose& pose );

	/// @copydoc ThermodynamicObserver::initialize_simulation
	virtual void initialize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	/// @copydoc ThermodynamicObserver::observe_after_metropolis
	virtual	void observe_after_metropolis(
		MetropolisHastingsMover const & metropolis_hastings_mover
	);

protected:

	/// @brief Pure virtual method called to write a model to the output file.
	virtual void 	write_model(
		core::pose::Pose const & pose,
		MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	) = 0;

private:
	core::Size stride_;
	core::Size model_count_;
	core::Size step_count_;
	core::Size cache_limit_;
	std::string file_name_;
	bool cumulate_jobs_;
	bool cumulate_replicas_;

	static bool options_registered_;
}; // TrajectoryRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_TrajectoryRecorder_HH
