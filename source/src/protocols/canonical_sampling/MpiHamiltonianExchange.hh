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

#ifndef INCLUDED_protocols_canonical_sampling_MpiHamiltonianExchange_hh
#define INCLUDED_protocols_canonical_sampling_MpiHamiltonianExchange_hh

#ifdef USEMPI
#include <mpi.h> //keep first
#endif

// Unit Headers
#include <protocols/canonical_sampling/MpiHamiltonianExchange.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace canonical_sampling {

class MpiHamiltonianExchange : public TemperatureController {
	typedef TemperatureController Parent;
	typedef utility::vector1< core::Size > GridCoord;
	typedef utility::vector1< GridCoord > Grid;

public:

	MpiHamiltonianExchange();

	//important d'tor to delete some C-style arrays
	~MpiHamiltonianExchange();

	MpiHamiltonianExchange( MpiHamiltonianExchange const& );

	MpiHamiltonianExchange& operator=( MpiHamiltonianExchange const& );

	virtual
	void apply( core::pose::Pose& ) {};

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

	/// @brief callback executed before any Monte Carlo trials
	virtual void
	initialize_simulation(
  	core::pose::Pose & pose,
		MetropolisHastingsMover const & mover,
		core::Size cycle
	);

	core::Real
	temperature_move(
		core::pose::Pose & pose,
		MetropolisHastingsMover & mover,
		core::Real score
	);

	/// @brief callback executed after all Monte Carlo trials
	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		MetropolisHastingsMover const & metropolis_hastings_mover
	);

	void show( std::ostream& ) const;

	void next_exchange_schedule();

protected:
	void set_defaults();
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief return to uninitialized status
	void clear();

	/// @brief initialize temperatures and weights from file, return false if IO error occurrs
	virtual bool init_from_file( std::string const& filename );


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
	void setup_exchange_schedule();
	void find_exchange_partner( int& partner, bool& is_master );
	///@brief small helper function; compute unique key out of (z1, z2, ... , zN) excluding zD
	core::Size coord2key(
    GridCoord const& coord,
		GridCoord const& max_coord,
		Size exclude_dim = 0
	);
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

	/// the different score-funcitons
	utility::vector1< core::scoring::ScoreFunctionOP > hamiltonians_;

	// exchange schedule determines which nodes try to swap at a given trial phase
	typedef utility::vector1<std::pair<int, int> > ExchangeSchedule;
	utility::vector0< ExchangeSchedule > exchange_schedules_;
	core::Size current_exchange_schedule_;

	// exchange happens between neighbouring cells in a D dimensional grid
	// in temperature replica-exchange usually D=1.
	Grid exchange_grid_;
	core::Size exchange_grid_dimension_;
	bool successfully_initialized_;
}; //end MpiHamiltonianExchange

/// @brief Test IO operator for debug and Python bindings
std::ostream& operator << ( std::ostream & os, MpiHamiltonianExchange const& );

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_MpiHamiltonianExchange_HH
