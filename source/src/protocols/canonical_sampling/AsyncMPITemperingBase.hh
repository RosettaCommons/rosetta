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
/// @author  Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_BiasEnergy_tmpl_hh
#define INCLUDED_protocols_canonical_sampling_BiasEnergy_tmpl_hh

#ifdef USEMPI
#include <mpi.h> //keep first
#endif

// Unit Headers
#include <protocols/canonical_sampling/TemperingBase.hh>

//#include <protocols/canonical_sampling/AsyncMPITemperingBase.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <protocols/canonical_sampling/MultiTemperatureTrialCounter.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

/// @details
class AsyncMPITemperingBase : public protocols::canonical_sampling::TemperingBase {
	typedef TemperingBase Parent;
	class ExchangeRequest {
	public:
		ExchangeRequest() :
			rank_( -1 ),
			level_( 0 ),
			requested_level_( 0 )
		{}

		ExchangeRequest( int rank, core::Size level, core::Size requested_level ) :
			rank_( rank ),
			level_( level ),
			requested_level_( requested_level )
		{}
		void show( std::ostream& ) const;
		int rank_;
		core::Size level_;
		core::Size requested_level_;
	};
	friend std::ostream& operator<< (std::ostream& os, AsyncMPITemperingBase::ExchangeRequest const& er );

	typedef std::map< core::Size, ExchangeRequest > ExchangeRequestMap;
public:

	AsyncMPITemperingBase();
	AsyncMPITemperingBase( AsyncMPITemperingBase const& );
	AsyncMPITemperingBase& operator=( AsyncMPITemperingBase const& );
	~AsyncMPITemperingBase() override;


	void apply( core::pose::Pose& ) override {};

	/// @brief callback executed before any Monte Carlo trials
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	) override;


	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size level,
		core::Real temperature,
		core::Size cycle
	) override;

	/// @brief callback executed after all Monte Carlo trials

	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	) override;


	bool finished_simulation( core::Size trials, core::Size ntrials) override;

protected:
#ifdef USEMPI
	MPI_Comm const& mpi_comm() const {
		return mpi_comm_;
	}

	void set_mpi_comm( MPI_Comm const& );

#endif


	bool time_for_temp_move() override;

	int rank() {
		return rank_;
	}
	//returns next exchange level or 0, if next schedule doesn't have exchange
	virtual
	core::Size next_exchange_level() const = 0;

	void find_exchange_partner( int& partner, bool& is_master );

private:
	typedef utility::vector1< std::pair< ExchangeRequest, ExchangeRequest > > MatchList;
	void process_match_list( MatchList const& list );
	void send_swap_invites( ExchangeRequest const& request1, ExchangeRequest const& request2 );
	void process_matching_exchange_requests();
	bool receive_exchange_requests();
	void remove_cycles_from_request_list();
	void start_listening();
	void listen_for_rank( int other_rank );
	void allocate_buffers();
	void deallocate_buffers();
	// rank within mpi_comm_
	int rank_;
	int comm_size_;
	bool ready_for_exchange_;
	bool finished_;
	int* recv_level_buffer_;
	int* send_level_buffer_;
	int* swap_invite_send_buffer_;
	int* swap_invite_recv_buffer_;
	int* finished_recv_buffer_;
	int* finished_send_buffer_;
#ifdef USEMPI
	MPI_Comm mpi_comm_;
	MPI_Request send_level_request_;
	MPI_Request* recv_level_requests_;
	MPI_Request* swap_invite_send_requests_;
	MPI_Request swap_invite_recv_request_;
	MPI_Request finished_send_request_;
	MPI_Request* finished_recv_requests_;
	MPI_Request finished_ack_recv_request_;
#endif
	ExchangeRequestMap exchange_requests_;
}; //end AsyncMPITemperingBase

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_AsyncMPITemperingBase_HH
