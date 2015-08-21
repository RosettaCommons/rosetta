// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/ParallelTemperingMover.cc
/// @brief ParallelTempering methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/ParallelTempering.hh>
#include <protocols/canonical_sampling/ParallelTemperingCreator.hh>


// protocols headers
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <protocols/rosetta_scripts/util.hh>

//#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

// core headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/string.functions.hh>

// basic headers
#include <basic/prof.hh>

// C++ Headers
#include <cmath>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer tr( "protocols.canonical_sampling.ParallelTempering" );


bool protocols::canonical_sampling::ParallelTempering::options_registered_( false );

//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::canonical_sampling::ParallelTempering::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		Parent::register_options();
	}
}

namespace protocols {
namespace canonical_sampling {
using namespace core;

std::string
ParallelTemperingCreator::keyname() const {
	return ParallelTemperingCreator::mover_name();
}

protocols::moves::MoverOP
ParallelTemperingCreator::create_mover() const {
	return protocols::moves::MoverOP( new ParallelTempering );
}

std::string
ParallelTemperingCreator::mover_name() {
	return "ParallelTempering";
}

ParallelTempering::ParallelTempering() :
	rank_( -1 ),
	last_energies_( NULL ),
	rank2tlevel_( NULL ),
	tlevel2rank_( NULL ),
	start_time_(0),
	total_mpi_wait_time_(0)
{
#ifndef USEMPI
	utility_exit_with_message( "ParallelTempering requires MPI build" );
#endif
#ifdef USEMPI
	mpi_comm_ = MPI_COMM_NULL;
#endif
	set_defaults();
}

ParallelTempering::ParallelTempering( ParallelTempering const & other ) :
	Parent( other ),
	exchange_schedules_( other.exchange_schedules_ ),
	last_exchange_schedule_( other.last_exchange_schedule_ ),
	last_energies_( NULL ),
	rank2tlevel_( NULL ),
	tlevel2rank_( NULL ),
	start_time_( other.start_time_ ),
	total_mpi_wait_time_( other.total_mpi_wait_time_ )
{
#ifndef USEMPI
	utility_exit_with_message( "ParallelTempering requires MPI build" );
#endif
#ifdef USEMPI
	set_mpi_comm( other.mpi_comm() );
#endif
}

ParallelTempering& ParallelTempering::operator=( ParallelTempering const& other ) {
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

ParallelTempering::~ParallelTempering() {
	deallocate_buffers();
}


void
ParallelTempering::initialize_simulation(
	pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	tr.Trace << "ParallelTempering::initialize_simul1... " << std::endl;
	Parent::initialize_simulation(pose, metropolis_hastings_mover,cycle);
	tr.Trace << "ParallelTempering::initialize_simul2... " << std::endl;
#ifdef USEMPI
	set_mpi_comm( jd2::current_mpi_comm() );
#endif
	tr.Trace << "ParallelTempering::initialize_simul3... " << std::endl;
	Size const nlevels( n_temp_levels() );
	allocate_buffers( nlevels );
	setup_exchange_schedule( nlevels );
	set_current_temp( rank2tlevel_[rank_] );

	last_exchange_schedule_ = 0;
	start_time_ = clock() / basic::SHRINK_FACTOR;
	total_mpi_wait_time_ = 0;
	tr.Trace << "Initialized ParallelTempering! " << std::endl;
}

void
ParallelTempering::finalize_simulation(
	pose::Pose& pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & mhm
) {
	deallocate_buffers();
	Parent::finalize_simulation( pose, mhm );

	//if replica runs end at different times that leads to hanging in output for some reason
	// this seems to help.
#ifdef USEMPI
	MPI_Barrier( mpi_comm() );
#endif

	if ( rank() == 0 ) {
		tr << "Temperature Exchange Frequencies:" << std::endl;
		std::streamsize oldwidth( tr.width() );
		std::streamsize oldprecision( tr.precision() );
		std::ios_base::fmtflags oldflags( tr.flags() );
		tr.width(5);
		tr.precision(3);
		tr.setf(std::ios_base::fixed);
		for ( core::Size i=0; i<n_temp_levels()-1; ++i ) {
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

void ParallelTempering::setup_exchange_schedule( Size nlevels ) {
	tr.Trace << "ParallelTempering::setup_exchange_schedule for " << nlevels << std::endl;
	exchange_schedules_.clear();
	ExchangeSchedule list;

	//0<->1, 2<->3...
	list.clear();
	for ( int i=0; i < (int)nlevels-1; i+=2 ) {
		std::pair<int, int> elem(i, i+1);
		list.push_back(elem);
		if ( rank() == 0 ) {
			exchange_attempts_[elem] = 0;
			exchange_accepts_[elem] = 0;
		}
	}
	exchange_schedules_.push_back(list);

	//1<->2, 3<->4...
	list.clear();
	for ( int i=1; i<(int)nlevels-1; i+=2 ) {
		std::pair<int, int> elem(i, i+1);
		list.push_back(elem);
		if ( rank() == 0 ) {
			exchange_attempts_[elem] = 0;
			exchange_accepts_[elem] = 0;
		}
	}
	exchange_schedules_.push_back(list);
}

void ParallelTempering::allocate_buffers( core::Size nlevels ) {
	tr.Trace << "ParallelTempering::allocate_buffers for " << nlevels << std::endl;
	deallocate_buffers();
	last_energies_ = new double[nlevels];
	rank2tlevel_ = new int[nlevels];
	tlevel2rank_ = new int[nlevels];
	for ( Size i=0; i<nlevels; ++i ) {
		rank2tlevel_[ i ] = i+1; //tlevels enumerate 1..N
		tlevel2rank_[ i ] = i; //ranks enumerate 0...N-1
	}
}

void ParallelTempering::deallocate_buffers() {
	if ( last_energies_ ) {
		delete [] last_energies_;
		delete [] rank2tlevel_;
		delete [] tlevel2rank_;
	}
	last_energies_ = NULL;
	rank2tlevel_ = NULL;
	tlevel2rank_ = NULL;
}

core::Real
ParallelTempering::temperature_move( core::Real MPI_ONLY( score ) ) {
	check_temp_consistency();
	if ( !time_for_temp_move() ) return temperature();

	// Size const nlevels( n_temp_levels() );
#ifdef USEMPI
	//get infomation
	double last_energy = score;
	clock_t time_before_MPI = clock();
	MPI_Gather(&last_energy, 1, MPI_DOUBLE, last_energies_, 1, MPI_DOUBLE, 0, mpi_comm() );

	  //change the T_tag and T_rev at node0
	if ( rank()==0 ) shuffle_temperatures( last_energies_ );

	//public the new T_tag
	int new_tlevel;
	MPI_Scatter(rank2tlevel_, 1, MPI_INT, &new_tlevel, 1, MPI_INT, 0, mpi_comm() );
	total_mpi_wait_time_ += ( clock() - time_before_MPI ) / basic::SHRINK_FACTOR;
	set_current_temp( new_tlevel );
#endif
	return temperature();
}

void
ParallelTempering::shuffle_temperatures( double *energies ) {
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

		if ( numeric::random::rg().uniform() < std::min( 1.0, std::exp( std::max(-40.0, -delta) ) ) ) {
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
		}
	}
}

std::string
ParallelTempering::get_name() const
{
	return "ParallelTempering";
}

protocols::moves::MoverOP
ParallelTempering::clone() const
{
	return protocols::moves::MoverOP( new protocols::canonical_sampling::ParallelTempering(*this) );
}

protocols::moves::MoverOP
ParallelTempering::fresh_instance() const
{
	return protocols::moves::MoverOP( new ParallelTempering );
}

void
ParallelTempering::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	pose::Pose const & pose
) {
	Parent::parse_my_tag( tag, data, filters, movers, pose );
}


/// handling of options including command-line
void ParallelTempering::set_defaults() {
}

/// @brief Assigns user specified values to primitive members using command line options
void ParallelTempering::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	Parent::init_from_options();
}

#ifdef USEMPI
void ParallelTempering::set_mpi_comm( MPI_Comm const& mpi_comm ) {
	if ( mpi_comm != MPI_COMM_NULL ) {
		tr.Trace << "ParallelTempering::Duplicate mpi-communicator" << std::endl;
		MPI_Comm_dup( mpi_comm, &mpi_comm_ );
		MPI_Comm_rank( mpi_comm_, &rank_ );
		int size;
		MPI_Comm_size( mpi_comm_, &size );
		runtime_assert( size == (int) n_temp_levels() );
	} else {
		tr.Trace << "ParallelTempering::Duplicate mpi-communicator" << std::endl;
		mpi_comm_ = MPI_COMM_NULL;
	}
}
#endif

} //moves
} //protocols

