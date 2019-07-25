// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/random/init_random_generator.cc
/// @brief  Functions to initialize the random generator.
/// @author Originally by Sergey Lyskov.
/// @modified Moved to basic from core by Vikram K. Mulligan (vmulligan@flatironinstitute.org).  Added multithreading support
/// to ensure that all threads have unique random seeds.


#ifdef USEMPI
#include <mpi.h> // Must go first
#endif

//Prototypes
#include <basic/random/init_random_generator.hh>

// STL headers
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/random/RandomGeneratorSettings.hh>
#include <basic/options/option.hh>

// Numeric headers
#include <numeric/random/random.hh>


static basic::Tracer TR( "basic.random.init_random_generator" );

namespace basic {
namespace random {

void
init_random_number_generators(
#ifdef MULTI_THREADED
	platform::Size const curthread /*=0*/
#endif
) {
	basic::random::RandomGeneratorSettings random_generator_settings;
	random_generator_settings.initialize_from_options( basic::options::option );
#ifdef MULTI_THREADED
	int const real_seed ( determine_random_number_seed( random_generator_settings, curthread ) );
#else
	int const real_seed ( determine_random_number_seed( random_generator_settings ) );
#endif
	init_random_generators(real_seed, random_generator_settings.rng_type() );

	// Seed default random generator.  This will hopefully expose all code that use
	// non-approved random methods -- assuming that code is invoked in an integration
	// test
	srand( time(nullptr) );
}

int determine_random_number_seed(
	basic::random::RandomGeneratorSettings const & rgs
#ifdef MULTI_THREADED
	,
	platform::Size const curthread
#endif
) {
	using namespace numeric::random;

	int real_seed;
	int seed = rgs.seed();

	if ( rgs.const_seed() ) {
		real_seed = rgs.seed() + rgs.seed_offset();
#ifdef USEMPI
		{ // scope
			/// Give a different RNG seed to each processor
			int mpi_rank( 0 );
			MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
			real_seed += mpi_rank;
		}
#endif
#ifdef MULTI_THREADED
		adjust_seed_for_thread( real_seed, curthread );
#endif //MULTI_THREADED
		TR << utility::CSI_Red() << utility::CSI_Underline() << "Constant seed mode" << utility::CSI_Reset() << ", seed=" << seed << " seed_offset=" << rgs.seed_offset()
			<< " real_seed=" << real_seed;
#ifdef MULTI_THREADED
		TR << " thread_index=" << curthread;
#endif
		TR << std::endl;
	} else {
#if (defined WIN32) && (!defined WIN_PYROSETTA)
		bool const on_windows_platform = true;
#else
		bool const on_windows_platform = false;
#endif
		// attempt to open rng device, if failure then fallback to time
		std::ifstream random_device( rgs.random_device_name().c_str(), std::ios::in | std::ios::binary );
		if ( ( random_device.fail() && !on_windows_platform ) || rgs.use_time_as_seed() ) {
			if ( !rgs.use_time_as_seed() ) {
				// notify user that opening rng device has failed
				TR << "NOTICE: rng device failure, using time as seed" << std::endl;
			}
			random_device.close();

			//iwd  When using time-based seed on a cluster, seed_offset is usually from 0 to num_processes.
			//iwd  If jobs start a few seconds apart, a simple sum of seed and seed_offset can lead to collisions.
			//iwd  Thus we multiply the time by some number larger than the largest expected value of seed_offset.
			//iwd  If anyone's using this on more than 1000 machines at a time, we need to increase this value.
			//iwd  (Rosetta++ used a multiplier of 20, which helps some, but is nonetheless too small.)
			seed = time(nullptr);
			//seed = seed%10000; // PB-- USE THIS ON OUR CLUSTER TO GET UNIQUE RUNS
			//real_seed = seed + rgs.seed_offset();
			real_seed = 1000*seed + rgs.seed_offset();

#ifdef USEMPI
			// When we use MPI and time-based seeds on a cluster, adjust the RNG seed so that it is the seed of the head node
			// (the node with an mpi rank of zero) plus the rank of the processer. This or is garentees that each node will
			// have a unique seed.

			/// get the processor rank
			int mpi_rank( 0 );
			MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

			// To avoid deadlock, this should only be called if the calling code has not disabled the
			// mpi_bcast call -- this code will be invoked in the initial core::init(...) call, but
			// should not be called in a context in which the RNG needs to be set mid-run (e.g. in
			// a multi-threaded-MPI context.)
			if (
				rgs.mpi_bcast_seed_from_node0()
#ifdef MULTI_THREADED
				&& curthread == 0 //Only thread zero is allowed to do MPI communication
#endif //MULTI_THREADED
			) {
				// set the real_seed of each processor to the real seed of the head node
				MPI_Bcast( &real_seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
			}

			// adjust the real seed based on the rank
			real_seed += mpi_rank;
#endif //USEMPI

#ifdef MULTI_THREADED
			adjust_seed_for_thread( real_seed, curthread );
#endif //MULTI_THREADED

			TR << "'Time' seed mode, seed=" << seed << " seed_offset=" << rgs.seed_offset()
				<< " real_seed=" << real_seed;
#ifdef MULTI_THREADED
			TR << " thread_index=" << curthread;
#endif
			TR << std::endl;
		} else {
			// grab seeds from device
			uint32_t unsigned_32bit_seed;

			std::string random_device_name = rgs.random_device_name();

#if (defined WIN32) && (!defined WIN_PYROSETTA)
			// windows random device name
			random_device_name = "CryptGenRandom";

			// unique name for key container
			ostringstream_t key_container_name;
			key_container_name << "arzetta-" << GetCurrentProcessId() << "-" << GetCurrentThreadId();

			// init cryptographic provider, creates key container
			HCRYPTPROV hCryptProv = 0;

			if ( !CryptAcquireContext( &hCryptProv, key_container_name.str().c_str(), NULL, PROV_RSA_AES, CRYPT_NEWKEYSET ) ) {
				TR.Fatal << "CryptAcquireContext unable to acquire cryptographic provider!" << std::endl;
				TR.Fatal << std::hex << GetLastError() << std::endl;
				std::exit( 1 );
			}

			// grab the random number seed from CryptGenRandom
			BYTE pbData[ 4 ];
			if ( !CryptGenRandom( hCryptProv, 4, pbData ) ) {
				TR.Fatal << "Unable to obtain random number seed using CryptGenRandom!" << std::endl;
				TR.Fatal << std::hex << GetLastError() << std::endl;
				std::exit( 1 );
			}

			// store seed in 'unsigned_32bit_seed'
			unsigned_32bit_seed = pbData[ 0 ] + ( pbData[ 1 ] << 8 ) + ( pbData[ 2 ] << 16 ) + ( pbData[ 3 ] << 24 );

			// release cryptographic provider handle
			if ( !CryptReleaseContext( hCryptProv, 0 ) ) {
				TR.Fatal << "CryptReleaseContext failed to release cryptographic provider!" << std::endl;
				TR.Fatal << std::hex << GetLastError() << std::endl;
				std::exit( 1 );
			}

			// delete key container
			if ( !CryptAcquireContext( &hCryptProv, key_container_name.str().c_str(), NULL, PROV_RSA_AES, CRYPT_DELETEKEYSET ) ) {
				TR.Fatal << "CryptAcquireContext failed to delete key container!" << std::endl;
				TR.Fatal << std::hex << GetLastError() << std::endl;
				std::exit( 1 );
			}

#else
			random_device.read( reinterpret_cast< char * >( &unsigned_32bit_seed ), sizeof( uint32_t ));
#endif
			random_device.close();

			// calculate actual seeds
			seed = static_cast< int >( unsigned_32bit_seed );
			real_seed = seed + rgs.seed_offset();

#ifdef USEMPI
			// Although not as critical with device-based seeding as with time-based clusters, when we use MPI and
			// device-based seeds on a cluster, adjust the RNG seed so that it is the seed of the head node (the node with an
			// mpi rank of zero) plus the rank of the processer. This or is garentees that each node will have a unique
			// seed. Different OSs impliment their RNG differently and as we use increassingly large numbers of processors
			// this may become an issue (but probaly not).

			/// get the processor rank
			int mpi_rank( 0 );
			MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

			// To avoid deadlock, this should only be called if the calling code has not disabled the
			// mpi_bcast call -- this code will be invoked in the initial core::init(...) call, but
			// should not be called in a context in which the RNG needs to be set mid-run (e.g. in
			// a multi-threaded-MPI context.)
			if (
				rgs.mpi_bcast_seed_from_node0()
#ifdef MULTI_THREADED
				&& curthread == 0 //Only thread zero can participate in a broadcast.
#endif
			) {
				// set the real_seed of each processor to the real seed of the head node
				MPI_Bcast( &real_seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
			}

			// adjust the real seed based on the rank
			real_seed += mpi_rank;
#endif

#ifdef MULTI_THREADED
			adjust_seed_for_thread( real_seed, curthread );
#endif //MULTI_THREADED

			// log seeds
			TR << "'RNG device' seed mode, using '" << random_device_name << "', seed=" << seed << " seed_offset=" << rgs.seed_offset()
				<< " real_seed=" << real_seed;
#ifdef MULTI_THREADED
			TR << " thread_index=" << curthread;
#endif
			TR << std::endl;
		}

	}

	return real_seed;
}

#ifdef MULTI_THREADED
/// @brief Given a seed that has already been adjusted for MPI (if this is the MPI build),
/// adjust for the current thread index.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
adjust_seed_for_thread(
	int & real_seed,
	platform::Size const curthread
) {
	//Ensure that each thread also has a different random seed.
#ifdef USEMPI
	//If we're using MPI, figure out the number of MPI ranks, and add (thread_index)*(num_ranks) to the seed.
	int total_mpi_processes(0);
	MPI_Comm_size( MPI_COMM_WORLD, &total_mpi_processes );
	debug_assert( total_mpi_processes > 0 );
	real_seed += curthread * total_mpi_processes;
#else //!USEMPI
	//If we're not using MPI, just add the thread index to the seed.
	real_seed += curthread;
#endif //USEMPI
}
#endif //MULTI_THREADED

/// @brief Initialize random generator systems (and send debug io to tracer with seed/mode info).
void init_random_generators( int const start_seed, std::string const & RGtype )
{
	TR << "RandomGenerator:init: Normal mode, seed=" << start_seed <<
		" RG_type=" << RGtype << std::endl;

	numeric::random::rg().set_seed( RGtype, start_seed );
}

} // namespace random
} // namespace basic
