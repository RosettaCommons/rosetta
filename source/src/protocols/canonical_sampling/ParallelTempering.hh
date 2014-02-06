// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/MetropolisHastingsMover.hh
/// @brief Implement replica exchange in the MetropolisHastingsMover Framework.
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_ParallelTempering_hh
#define INCLUDED_protocols_canonical_sampling_ParallelTempering_hh

#ifdef USEMPI
#include <mpi.h> // Keep first...don't know why...
#endif

// Unit Headers
#include <protocols/canonical_sampling/ParallelTempering.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

// C Headers
#include <time.h>

namespace protocols {
namespace canonical_sampling {


/// @brief Manage several replicas with periodic temperature swapping.
///
/// @details Parallel tempering is a standard Monte Carlo technique.  The idea 
/// is to simulate the same system at several different temperatures, and to 
/// periodically swap the temperatures between the simulations.  The 
/// simulations at higher temperatures sample more broadly, while the 
/// simulations at lower temperatures sample more deeply.  By swapping 
/// temperatures (such that detailed balance is maintained) the best of both 
/// worlds is obtained.  In order for a parallel tempering simulation to be 
/// successful, the replicas need to have good overlap.  That is to say, there 
/// has to be a reasonable chance of observing the system at temperature levels 
/// one above and one below the current.  If this is not the case, most of the 
/// temperature moves will be rejected and more levels need to be added.
///
/// This implementation requires multi-threading, as provided by the MPI 
/// library.  In addition, the number of replicas must exactly match the number 
/// of threads available to the job.  If there is a mismatch an error will be 
/// thrown.  On some systems it can be very challenging to compile Rosetta with 
/// support for MPI, but you can usually find help on the RosettaCommons forum.  
/// Also be aware that not every class in this module is thread-safe.  For 
/// example, the DbTrajectoryReporter will attempt to write to the database 
/// from several threads at once, which is forbidden when using SQLite.
/// 
/// Note that he only way to set the temperature range used for simulated 
/// annealing is to use the command line.  This is something I'd be interested 
/// in changed at some point, but for right now it's not a deal-breaker. The 
/// relevant options are:
///
/// @code
/// -tempering::temp::range <low> <high>
/// -tempering::temp::low <low> -tempering::temp::high <high>
/// -tempering::temp::levels <levels>
/// @endcode

class ParallelTempering : public TemperatureController {
	typedef TemperatureController Parent;

public:

	/// @brief Default constructor.
	ParallelTempering();

	/// @brief Non-trivial destructor used to free some C-style arrays.
	~ParallelTempering();

	/// @brief Copy constructor.
	ParallelTempering( ParallelTempering const& );

	/// @brief Assignment operator.
	ParallelTempering& operator=( ParallelTempering const& );

	virtual
	void apply( core::pose::Pose& ) {}

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

	virtual void
	initialize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & mover,
		core::Size cycle
	);

	/// @brief Execute a temperature move if neccesary
	/// @details This method blocks until it is reached by every thread.  Then 
	/// the last energies recorded in each replica are communicated to all the 
	/// other replicas, and the root thread coordinates the temperature swap.  
	/// Once all of this is done, the simulation continues.  Note that this can 
	/// produce a significant amount of dead time if called too often.  At the 
	/// end of the simulation, the amount of time spent waiting for MPI get 
	/// reported to the tracer.  Use this information to decide how often these 
	/// moves should be attempted.
	/// @see shuffle_temperatures()
	core::Real
	temperature_move(
			core::pose::Pose & pose,
			MetropolisHastingsMover & mover,
			core::Real score);

	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & metropolis_hastings_mover
	);

protected:
	
#ifdef USEMPI
	/// @brief Return a handle to the MPI communication group being used for this 
	/// job.
	MPI_Comm const& mpi_comm() const {
		return mpi_comm_;
	}

	/// @brief Keep a handle to the MPI communication group being used for this 
	/// job.
	void set_mpi_comm( MPI_Comm const& );

#endif

	/// @brief Return the rank of the MPI process being used to simulate this 
	/// specific temperature level.
	int rank() {
		return rank_;
	}

private:

	/// @brief Allocate the arrays used for MPI communication.
	/// @details Normally it's better to use vectors instead of C-style arrays, 
	/// because manually allocating and deallocating arrays is error-prone.  
	/// However, to use MPI there's no alternative.
	void allocate_buffers( core::Size );
	
	/// @brief Deallocate the arrays used for MPI communication.
	/// @details Normally it's better to use vectors instead of C-style arrays, 
	/// because manually allocating and deallocating arrays is error-prone.  
	/// However, to use MPI there's no alternative.
	void deallocate_buffers();

	/// @brief Setup the algorithm used to propose temperature swaps.
	/// @details This algorithm is surprisingly complex.  The gist of it goes 
	/// like this.  Half of the time, swaps between levels 1->2, 3->4, and so on 
	/// are considered.  The other half of the time, swaps between levels 2->3, 
	/// 4->5, and so on are considered.  In this way, every swap is considered 
	/// with the same frequency, and never no level will ever get swapped twice 
	/// in one move.  This method setups up a data structure to implement this 
	/// algorithm an populates it based on the number of replicas being used.
	/// @see shuffle_temperatures()
	void setup_exchange_schedule( Size nlevels );

	/// @brief Attempt swaps between half of the temperature level pairs.
	/// @details Swaps are accepted like any Monte Carlo move: using the 
	/// Metropolis criterion.
	/// @see setup_exchange_schedule()
	void shuffle_temperatures(
			MetropolisHastingsMover & mover, double *energies);

private:
#ifdef USEMPI
	MPI_Comm mpi_comm_;
#endif
	// rank within mpi_comm_
	int rank_;
	typedef utility::vector1<std::pair<int, int> > ExchangeSchedule;
	utility::vector0< ExchangeSchedule > exchange_schedules_;
	core::Size last_exchange_schedule_;
	// C-style arrays for communication in MPI
	double *last_energies_;
	int *rank2tlevel_;
	int *tlevel2rank_;
	std::map<std::pair<int, int>, core::Size> exchange_attempts_;
	std::map<std::pair<int, int>, core::Size> exchange_accepts_;
	clock_t start_time_;
	clock_t total_mpi_wait_time_;

}; //end ParallelTempering

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_ParallelTempering_HH
