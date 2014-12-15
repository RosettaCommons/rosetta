// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_BiasEnergy_hh
#define INCLUDED_protocols_canonical_sampling_BiasEnergy_hh


#ifdef USEMPI
#include <mpi.h> //keep first
#endif

// Unit Headers
#include <protocols/canonical_sampling/BiasEnergy.fwd.hh>

// Module Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobOutputterObserver.hh>

//for parse_my_tag
#include <utility/tag/Tag.hh>
#include <utility/pointer/ReferenceCount.hh>


// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace canonical_sampling {

class BiasEnergy : public ThermodynamicObserver, public protocols::jd2::JobOutputterObserver {

	template< typename T>
	class Histogram;

	template< typename T>
	class Histogram : public utility::pointer::ReferenceCount {
		typedef utility::pointer::ReferenceCount Parent;
		typedef T ValueType;

	public:
		typedef utility::pointer::shared_ptr< Histogram<T> > OP;
		typedef utility::pointer::shared_ptr< Histogram const > COP;

		Histogram( core::Real grid_min, core::Real grid_max, core::Size ngrid_cells_ = 100 );
		Histogram( Histogram const& );
		Histogram& operator=( Histogram const& );

		~Histogram();

#ifdef USEMPI
		void mpi_exchange( int partner, MPI_Comm const& mpi_comm );
#endif

		bool check_range( core::Real val ) const;
		core::Size cell_index( core::Real val ) const;
		core::Size size() const { return ngrid_cells_; };
		ValueType& at( core::Size index ) { return data_[ index - 1 ]; };
		ValueType at( core::Size index ) const { return data_[ index - 1 ]; };

		//switch to data in recv_buf which had previously been pre-staged with mpi_exchange
		void toggle();
		void reset();

		void write_to_stream( std::ostream& ) const;
		void write_to_string( std::string& ) const;

	private:
		void copy_data( Histogram const& );
		ValueType *data_; //use C-style array which allows direct sending and receiving with MPI
		ValueType *recv_buf_;
		core::Real grid_min_;
		core::Real grid_max_;
		core::Real delta_grid_;
		core::Size ngrid_cells_;
	};
///@details

 public:
	BiasEnergy();
	BiasEnergy( core::Size stride, core::Real omega, core::Real gamma );
	virtual ~BiasEnergy();
	virtual core::Real evaluate( core::pose::Pose const& ) const;

	//for trajectory output...
	virtual
	void add_values_to_job( core::pose::Pose const& pose, protocols::jd2::Job & ) const;

#ifdef USEMPI
	//in replica exchange
	virtual core::Real update_and_evaluate_replica(
    int exchange_partner,
		MPI_Comm const& mpi_comm,
		core::pose::Pose const& pose
	);
#endif

	virtual void update( core::pose::Pose const& );

	void swap_replicas();

	//write to string instead of stream to reduce buffer caused delay
	void write_to_string( std::string & str ) const;

	void set_temperature( core::Real setting );

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose &,
		MetropolisHastingsMover const &,
		core::Size //non-zero if trajectory is restarted
	);

	virtual
	void
	observe_after_metropolis(
		MetropolisHastingsMover const &
	) {};

	virtual
	void
	finalize_simulation(
		core::pose::Pose &,
		MetropolisHastingsMover const &
	);

	virtual
	bool
	restart_simulation(
	  core::pose::Pose &,
		MetropolisHastingsMover &,
		core::Size& cycle,
		core::Size& temp_level,
		core::Real& temperature
	);

 protected:
	virtual core::Real extract_collective_var( core::pose::Pose const& ) const = 0;

 private:
	//Histogram<float>::OP bias_grid_;
	//Histogram<int>::OP count_grid_;
	// To help with auto code rewriting:
	typedef Histogram<float>::OP HistogramFloatOP;
	typedef Histogram<int>::OP HistogramIntOP;
	HistogramFloatOP bias_grid_;
	HistogramIntOP count_grid_;

	core::Size step_counter_;

	core::Size stride_;
	core::Real omega_;
	core::Real gamma_;
	core::Real deltaT_;
	core::Real temperature_;

	core::Real default_bias_;

	std::string grid_output_file_;
};

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_BiasEnergy_HH
