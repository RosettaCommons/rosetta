// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_multithreaded_packing_performance.cc
/// @brief An app to test the multithreaded performance of the packer.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// protocol headers

// core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetsFactory.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/thread_manager/RosettaThreadManager.hh>

// C++ headers
#include <chrono>

static basic::Tracer TR("test_multithreaded_packing_performance");

OPT_KEY( Integer, replicates )

/// @brief Indicate which options are relevant.
void
register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( multithreading::total_threads );
	option.add_relevant( in::file::s );
	NEW_OPT( replicates, "The number of replicates to perform for each packing attempt.  Default 10.", 10 );
}

#ifdef MULTI_THREADED
/// @brief Actually carry out the multithreading test.
void
do_test(
	core::pose::PoseCOP master_pose,
	core::Size const nthreads,
	core::Size const nreplicates
) {
	basic::thread_manager::RosettaThreadManager::get_instance()->launch_threads();

	core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function( true ) );
	core::pack::rotamer_set::RotamerSetsOP rotsets( core::pack::rotamer_set::RotamerSetsFactory::create_rotamer_sets( *master_pose ) );

	utility::vector1< utility::vector1< core::Size > > times_by_threadcount( nthreads );

	for ( core::Size ithread(nthreads); ithread >= 1; --ithread ) { // Do nthreads to 1 threads.
		TR << "Trying packing on " << ithread << " threads (of " << nthreads << ")." << std::endl;
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( *master_pose, ithread ) );
		utility::vector1< core::Size > times_for_current_threadcount(nreplicates);
		for ( core::Size ireplicate(1); ireplicate <= nreplicates; ++ireplicate ) { // Do 1 to nreplicates replicates.

			//Copy the pose to create a working copy:
			core::pose::PoseOP pose( master_pose->clone() );

			core::pack::interaction_graph::AnnealableGraphBaseOP intgraph(nullptr);

			std::chrono::time_point< std::chrono::high_resolution_clock > const start( std::chrono::high_resolution_clock::now() );
			core::pack::pack_rotamers_setup( *pose, *sfxn, task, rotsets, intgraph );
			std::chrono::time_point< std::chrono::high_resolution_clock > const end( std::chrono::high_resolution_clock::now() );

			TR << "Replicate " << ireplicate << " of " << nreplicates << " for " << ithread << " interaction graph setup threads completed in " << std::chrono::duration_cast< std::chrono::microseconds >( end-start ).count() << " microseconds." << std::endl;
			times_for_current_threadcount[ireplicate] = static_cast<core::Size>( std::chrono::duration_cast< std::chrono::microseconds >( end-start ).count() );
		}
		times_by_threadcount[ithread] = times_for_current_threadcount;
	}

	TR << "\n\tTHREADS\n";
	TR << "REPLICATE";
	for ( core::Size ithread(1); ithread<=nthreads; ++ithread ) {
		TR << "\t" << ithread;
	}
	TR << "\n";
	for ( core::Size ireplicate(1); ireplicate<=nreplicates; ++ireplicate ) {
		TR << ireplicate;
		for ( core::Size ithread(1); ithread<=nthreads; ++ithread ) {
			TR << "\t" << times_by_threadcount[ithread][ireplicate];
		}
		TR << "\n";
	}

	// Compute mean of inverse:
	TR << "Mean_of_inv(igs/sec)";
	utility::vector1< core::Real > inv_averages( nthreads );
	for ( core::Size ithread(1); ithread<=nthreads; ++ithread ) {
		core::Size sum(0);
		for ( core::Size ireplicate(1); ireplicate<=nreplicates; ++ireplicate ) {
			sum += times_by_threadcount[ithread][ireplicate];
		}
		inv_averages[ithread] = static_cast<core::Real>(nreplicates)*1.0e6/static_cast<core::Real>(sum);
		TR << "\t" << inv_averages[ithread];
	}
	TR << "\n";

	// Compute std. err of mean of inverse:
	TR << "STDERR[Mean_of_inv]";
	for ( core::Size ithread(1); ithread<=nthreads; ++ithread ) {
		core::Real accumulator(0.0);
		for ( core::Size ireplicate(1); ireplicate<=nreplicates; ++ireplicate ) {
			accumulator += std::pow(1.0e6/static_cast<core::Real>( times_by_threadcount[ithread][ireplicate] ) - inv_averages[ithread], 2);
		}
		accumulator /= static_cast<core::Real>(nreplicates);
		accumulator = std::sqrt( accumulator );
		TR << "\t" << accumulator;
	}

	TR << std::endl;
}
#endif //MULTI_THREADED

/// @brief Entry point for program execution.
int
#ifdef MULTI_THREADED
main( int argc, char * argv [] )
#else
main( int, char * [] )
#endif
{
	try {
#ifdef MULTI_THREADED
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		TR << "Starting test_multithreaded_packing_performance application." << std::endl;
		TR << "Pilot app created 20 May 2019 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitue.org)." << std::endl;

		runtime_assert_string_msg( option[ in::file::s ].user(), "Please specify an input PDB file with the -in:file:s option.");
		runtime_assert_string_msg( option[ in::file::s ]().size() == 1, "Please specify exactly one input PDB file with the -in:file:s option.");
		runtime_assert_string_msg( option[multithreading::total_threads].user(), "Please specify the number of threads to launch with the -multithreading:total_threads option." );

		core::Size const nthreads( option[multithreading::total_threads]() );
		core::Size const nreplicates( option[replicates]() );
		std::string const filename( option[in::file::s]()[1] );

		core::pose::PoseOP pose( core::import_pose::pose_from_file(filename) );
		runtime_assert_string_msg( pose != nullptr, "Problem constructing pose from data in \"" + filename + "\"." );

		do_test( pose, nthreads == 0 ? basic::thread_manager::RosettaThreadManager::total_threads() : nthreads, nreplicates );

		TR.flush();

#else
		utility_exit_with_message("This application cannot be run in the single-threaded compilation of Rosetta.  Please compile with the \"extras=cxx11thread\" option.");
#endif

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	TR << "Completed test_multithreaded_packing_performance.  Exiting with status 0 (no errors)." << std::endl;
	return 0;
}
