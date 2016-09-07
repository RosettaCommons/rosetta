// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/canonical_sampling/MetropolisHastingsMover.hh
/// @brief Implement replica exchange in the MetropolisHastingsMover Framework.
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_ParallelTempering_hh
#define INCLUDED_protocols_canonical_sampling_ParallelTempering_hh

#ifdef USEMPI
#include <mpi.h> //keep first
#endif

// Unit Headers
#include <protocols/canonical_sampling/ParallelTempering.fwd.hh>
#include <protocols/canonical_sampling/TemperingBase.hh>
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
#include <ctime>

namespace protocols {
namespace canonical_sampling {


/// @brief
/// @details The only way to set the temperature range used for simulated
/// annealing is to use the command line.  The relevant options are:
///
/// @code
/// -tempering::temp::range <low> <high>
/// -tempering::temp::low <low> -tempering::temp::high <high>
/// @endcode
class ParallelTempering : public protocols::canonical_sampling::TemperingBase {
	typedef TemperingBase Parent;
public:

	/// @brief Default constructor.
	ParallelTempering();

	/// @brief Non-trivial destructor used to free some C-style arrays.
	~ParallelTempering() override;

	/// @brief Copy constructor.
	ParallelTempering( ParallelTempering const& );

	/// @brief Assignment operator.
	ParallelTempering& operator=( ParallelTempering const& );

	
	void apply( core::pose::Pose& ) override {};

	
	std::string
	get_name() const override;

	protocols::moves::MoverOP
	clone() const override;

	
	protocols::moves::MoverOP
	fresh_instance() const override;

	
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	core::Real
	temperature_move( core::Real score) override;

	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle  //non-zero if trajectory is restarted
	) override;

	
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	) override;

protected:
	void set_defaults();
	void init_from_options() override;

#ifdef USEMPI
	MPI_Comm const& mpi_comm() const {
		return mpi_comm_;
	}

	void set_mpi_comm( MPI_Comm const& );

#endif

	int rank() {
		return rank_;
	}

private:
	void deallocate_buffers();
	void allocate_buffers( core::Size );
	void setup_exchange_schedule( Size nlevels );
	void shuffle_temperatures( double *energies );

private:
	static bool options_registered_;

public:
	static void register_options();


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
