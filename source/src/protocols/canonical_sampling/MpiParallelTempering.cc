// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/MpiParallelTemperingMover.cc
/// @brief MpiParallelTempering methods implemented
/// @author

// Unit headers
#include <protocols/canonical_sampling/MpiParallelTempering.hh>
#include <protocols/canonical_sampling/MpiParallelTemperingCreator.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

// Core headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/prof.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/string.functions.hh>
#include <numeric/random/random.hh>

// C++ Headers
#include <cmath>

using namespace std;
using namespace core;
using core::pose::Pose;
using protocols::moves::MoverOP;

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer tr( "protocols.canonical_sampling.MpiParallelTempering" );
static numeric::random::RandomGenerator RG(3227547);

namespace protocols {
namespace canonical_sampling {

string MpiParallelTemperingCreator::keyname() const { // {{{1
	return MpiParallelTemperingCreator::mover_name();
}

MoverOP MpiParallelTemperingCreator::create_mover() const { // {{{1
	return new MpiParallelTempering;
}

string MpiParallelTemperingCreator::mover_name() { // {{{1
	return "MpiParallelTempering";
}
// }}}1

MpiParallelTempering::MpiParallelTempering() // {{{1
	: rank_( -1 ),
	  last_energies_( NULL ),
	  rank2tlevel_( NULL ),
	  tlevel2rank_( NULL ),
	  start_time_(0),
	  total_mpi_wait_time_(0) {

#ifndef USEMPI
	utility_exit_with_message( "MpiParallelTempering requires MPI build" );
#endif
#ifdef USEMPI
	mpi_comm_ = MPI_COMM_NULL;
#endif
	set_defaults();
	type("parallel-tempering");
}

MpiParallelTempering::MpiParallelTempering(MpiParallelTempering const & other) // {{{1
	: Parent( other ),
	  exchange_schedules_( other.exchange_schedules_ ),
	  last_exchange_schedule_( other.last_exchange_schedule_ ),
	  last_energies_( NULL ),
	  rank2tlevel_( NULL ),
	  tlevel2rank_( NULL ),
	  start_time_( other.start_time_ ),
	  total_mpi_wait_time_( other.total_mpi_wait_time_ ) {

#ifndef USEMPI
	utility_exit_with_message( "MpiParallelTempering requires MPI build" );
#endif
#ifdef USEMPI
	set_mpi_comm( other.mpi_comm() );
#endif
}

MpiParallelTempering& MpiParallelTempering::operator=( // {{{1
		MpiParallelTempering const& other ) {

	if ( &other == this ) return *this;
	deallocate_buffers();
	Parent::operator=( other );
#ifdef USEMPI
	set_mpi_comm( other.mpi_comm() );
#endif
	Size const nlevels( n_temp_levels() );
	allocate_buffers( nlevels );
	runtime_assert( exchange_schedules_[ 0 ].size() == nlevels );
	return *this;
}

MpiParallelTempering::~MpiParallelTempering() { // {{{1
	deallocate_buffers();
}

void MpiParallelTempering::initialize_simulation( // {{{1
		pose::Pose & pose,
		MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle) {

	tr.Trace << "MpiParallelTempering::initialize_simul1... " << std::endl;
	Parent::initialize_simulation(pose, metropolis_hastings_mover,cycle);
	tr.Trace << "MpiParallelTempering::initialize_simul2... " << std::endl;
#ifdef USEMPI
	set_mpi_comm( jd2::current_mpi_comm() );
#endif
	tr.Trace << "MpiParallelTempering::initialize_simul3... " << std::endl;
	Size const nlevels( n_temp_levels() );
	allocate_buffers( nlevels );
	setup_exchange_schedule( nlevels );
	set_current_temp( rank2tlevel_[rank_] );

	last_exchange_schedule_ = 0;
	start_time_ = clock() / basic::SHRINK_FACTOR;
	total_mpi_wait_time_ = 0;
	tr.Trace << "Initialized MpiParallelTempering! " << std::endl;
}

Real MpiParallelTempering::temperature_move( // {{{1
		Pose&,
		MetropolisHastingsMover & MPI_ONLY(mover),
		Real MPI_ONLY(score)) {

	check_temp_consistency();
	if ( !time_for_temp_move() ) return temperature();

#ifdef USEMPI
	double last_energy = score;
	clock_t time_before_MPI = clock();

	// Whenever a node reaches MPI_Gather, it will immediately send its last
	// energy to the root node.  The root node will wait here until it has
	// received a message from every child node.  The child nodes will continue
	// on initially, but will wait once they reach MPI_Scatter.  Only after the
	// root has heard from all the children and sent out the new temperature
	// assignments will the simulation proceed.

	MPI_Gather(
			&last_energy, 1, MPI_DOUBLE,
			last_energies_, 1, MPI_DOUBLE,
			0, mpi_comm() );

	// Change the T_tag and T_rev at the root node.
	if ( rank()==0 ) shuffle_temperatures(mover, last_energies_ );

	int new_tlevel;

	MPI_Scatter(
			rank2tlevel_, 1, MPI_INT,
			&new_tlevel, 1, MPI_INT,
			0, mpi_comm() );

	total_mpi_wait_time_ += ( clock() - time_before_MPI ) / basic::SHRINK_FACTOR;
	set_current_temp( new_tlevel );
#endif
	return temperature();
}

void MpiParallelTempering::finalize_simulation( // {{{1
		pose::Pose& pose,
		MetropolisHastingsMover const & mhm) {

	deallocate_buffers();
	Parent::finalize_simulation( pose, mhm );

	if (rank() == 0) {
		tr << "Temperature Exchange Frequencies:" << std::endl;
		std::streamsize oldwidth( tr.width() );
		std::streamsize oldprecision( tr.precision() );
		std::ios_base::fmtflags oldflags( tr.flags() );
		tr.width(5);
		tr.precision(3);
		tr.setf(std::ios_base::fixed);
		for (core::Size i=0; i<n_temp_levels()-1; ++i) {
			std::pair<int, int> elem(i, i+1);
			//Original code (line below) fails on some versions of GCC (4.4.3 and 4.3.4 reported by users)
			//core::Real frequency(core::Real(exchange_accepts_[elem])/core::Real(exchange_attempts_[elem]));
			//this replacement code is reported to work for those compilers
			core::Real const a = exchange_accepts_[elem];
			core::Real const b = exchange_attempts_[elem];
			core::Real const frequency(a/b);
			tr << temperature(i+1) << " <-> " << temperature(i+2) << ": " << frequency
				 << " (" << exchange_accepts_[elem] << " of " << exchange_attempts_[elem] << ")" << std::endl;
		}
		tr.width(oldwidth);
		tr.precision(oldprecision);
		tr.flags(oldflags);
	}

	core::Real const clock_factor( ( (double) basic::SHRINK_FACTOR ) / CLOCKS_PER_SEC );
	clock_t total_time(clock()/basic::SHRINK_FACTOR - start_time_);
	core::Real fraction_waiting = total_mpi_wait_time_*clock_factor / ( total_time*clock_factor );
	tr << "Spent " << fraction_waiting*100 << "% time waiting for MPI temperature exchange ("
	   << total_mpi_wait_time_*clock_factor << " seconds out of " << total_time*clock_factor << " total)" << std::endl;
}

void MpiParallelTempering::allocate_buffers( core::Size nlevels ) { // {{{1
	tr.Trace << "MpiParallelTempering::allocate_buffers for " << nlevels << std::endl;
	deallocate_buffers();
	last_energies_ = new double[nlevels];
	rank2tlevel_ = new int[nlevels];
	tlevel2rank_ = new int[nlevels];
	for ( Size i=0; i<nlevels; ++i ) {
		rank2tlevel_[ i ] = i+1; //tlevels enumerate 1..N
		tlevel2rank_[ i ] = i;   //ranks enumerate 0...N-1
	}
}

void MpiParallelTempering::deallocate_buffers() { // {{{1
	if ( last_energies_ ) {
		delete [] last_energies_;
		delete [] rank2tlevel_;
		delete [] tlevel2rank_;
	}
	last_energies_ = NULL;
	rank2tlevel_ = NULL;
	tlevel2rank_ = NULL;
}

void MpiParallelTempering::setup_exchange_schedule( Size nlevels ) { // {{{1
	tr.Trace << "MpiParallelTempering::setup_exchange_schedule for " << nlevels << std::endl;
	exchange_schedules_.clear();
	ExchangeSchedule list;

	//0<->1, 2<->3...
	list.clear();
	for (int i=0; i < (int)nlevels-1; i+=2) {
		std::pair<int, int> elem(i, i+1);
		list.push_back(elem);
		if (rank() == 0) {
			exchange_attempts_[elem] = 0;
			exchange_accepts_[elem] = 0;
		}
	}
	exchange_schedules_.push_back(list);

	//1<->2, 3<->4...
	list.clear();
	for (int i=1; i<(int)nlevels-1; i+=2) {
		std::pair<int, int> elem(i, i+1);
		list.push_back(elem);
		if (rank() == 0) {
			exchange_attempts_[elem] = 0;
			exchange_accepts_[elem] = 0;
		}
	}
	exchange_schedules_.push_back(list);
}

void MpiParallelTempering::shuffle_temperatures( // {{{1
		MetropolisHastingsMover & mover, double *energies ) {

	last_exchange_schedule_ = ( last_exchange_schedule_ + 1 ) % 2;
	ExchangeSchedule const& ex( exchange_schedules_[ last_exchange_schedule_ ] );

	for ( Size i=1; i<=ex.size(); i++ ) {
		Size const rank1=tlevel2rank_[ex[i].first];
		Size const rank2=tlevel2rank_[ex[i].second];
		Real const invT1( 1.0 / temperature( rank2tlevel_[rank1] ) );
		Real const invT2( 1.0 / temperature( rank2tlevel_[rank2] ) );
		Real const deltaE( energies[rank2]-energies[rank1] );
		Real const delta( ( invT1 - invT2 ) * deltaE );

		++exchange_attempts_[ex[i]];
		mover.count_trial(type());

		if ( RG.uniform() < std::min( 1.0, std::exp( std::max(-40.0, -delta) ) ) ) {
			Size tmp;

			//Swap tlevel
			tmp=rank2tlevel_[rank1];
			rank2tlevel_[rank1]=rank2tlevel_[rank2];
			rank2tlevel_[rank2]=tmp;

			//Swap ranks
			tmp=tlevel2rank_[ex[i].first];
			tlevel2rank_[ex[i].first]=tlevel2rank_[ex[i].second];
			tlevel2rank_[ex[i].second]=tmp;

			++exchange_accepts_[ex[i]];
			mover.count_accepted(type());
			mover.count_energy_drop(type(), deltaE);
		}
	}
}

string MpiParallelTempering::get_name() const { // {{{1
	return "MpiParallelTempering";
}

MoverOP MpiParallelTempering::clone() const { // {{{1
	return new MpiParallelTempering(*this);
}

MoverOP MpiParallelTempering::fresh_instance() const { // {{{1
	return new MpiParallelTempering;
}

// }}}1

#ifdef USEMPI
void MpiParallelTempering::set_mpi_comm( MPI_Comm const& mpi_comm ) { // {{{1
	if ( mpi_comm != MPI_COMM_NULL ) {
		tr.Trace << "MpiParallelTempering::Duplicate mpi-communicator" << std::endl;
		MPI_Comm_dup( mpi_comm, &mpi_comm_ );
		MPI_Comm_rank( mpi_comm_, &rank_ );
		int size;
		MPI_Comm_size( mpi_comm_, &size );
		runtime_assert( size == (int) n_temp_levels() );
	} else {
		tr.Trace << "MpiParallelTempering::Duplicate mpi-communicator" << std::endl;
		mpi_comm_ = MPI_COMM_NULL;
	}
} // }}}1
#endif

} //moves
} //protocols

