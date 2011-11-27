// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/MetropolisHastingsMover.hh
/// @brief Implement replica exchange in the MetropolisHastingsMover Framework.
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_moves_ParallelTempering_hh
#define INCLUDED_protocols_moves_ParallelTempering_hh

#ifdef USEMPI
#include <mpi.h> //keep first
#endif

// Unit Headers
#include <protocols/moves/ParallelTempering.fwd.hh>
#include <protocols/moves/TemperingBase.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/TemperatureController.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace moves {

///@details
class ParallelTempering : public TemperingBase {
	typedef TemperingBase Parent;
public:

	ParallelTempering();

	//important d'tor to delete some C-style arrays
	~ParallelTempering();

	ParallelTempering( ParallelTempering const& );

	ParallelTempering& operator=( ParallelTempering const& );

	virtual
	void apply( core::pose::Pose& ) {};

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

	/// @brief execute the temperatur move ( called by observer_after_metropolis )
	/// returns the current temperatur in kT.
	core::Real
	temperature_move( core::Real score);

	/// @brief callback executed before any Monte Carlo trials
	virtual void
	initialize_simulation();

	/// @brief callback executed after all Monte Carlo trials
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

protected:
	void set_defaults();
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

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
/// ------------------ register cmdline options ---------------------------

private:
	static bool options_registered_;

public:
	static void register_options();

/// ---------------- member variables --------------------------

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



}; //end ParallelTempering

} //namespace moves
} //namespace protocols

#endif //INCLUDED_protocols_moves_ParallelTempering_HH
