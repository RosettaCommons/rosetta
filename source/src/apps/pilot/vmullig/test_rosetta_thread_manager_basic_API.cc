// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_rosetta_thread_manager_basic_API.cc
/// @brief An application to test the RosettaThreadManager's default API, which accepts a vector of work.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/types.hh>

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>
#include <numeric/random/random.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>

// multithreading headers
#ifdef MULTI_THREADED
#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#endif //MULTI_THREADED

static basic::Tracer TR_main("test_rosetta_thread_manager_basic_API");
static basic::Tracer TR_level1("test_rosetta_thread_manager_basic_API_level1");
static basic::Tracer TR_level2("test_rosetta_thread_manager_basic_API_level2");
static basic::Tracer TR_level3("test_rosetta_thread_manager_basic_API_level3");
static basic::Tracer TR_output("test_rosetta_thread_manager_basic_API_output");

OPT_KEY( Integer, level1_threads )
OPT_KEY( Integer, level2_threads )
OPT_KEY( Integer, level3_threads )
OPT_KEY( Boolean, whole_row_mode )

/// @brief Indicate which options are relevant.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( level1_threads, "The number of level one threads to request.  The actual number used may be less than the number requested depending on thread availability.", 1 );
	NEW_OPT( level2_threads, "The number of level two threads to request.  The actual number used may be less than the number requested depending on thread availability.", 1 );
	NEW_OPT( level3_threads, "The number of level three threads to request.  The actual number used may be less than the number requested depending on thread availability.", 1 );
	NEW_OPT( whole_row_mode, "If true, a single thread computes a whole row of the table.  False by default.  Cannot be used with the -level3_threads option.", false );

	option.add_relevant( multithreading::total_threads );
}

#ifdef MULTI_THREADED

/// @brief Function for level three threads: computing a single entry in the array
void level3_function(
	utility::vector1< utility::vector1< core::Size > > & random_numbers,
	core::Size const row,
	core::Size const column,
	basic::thread_manager::RosettaThreadAssignmentInfo const & thread_assignments
) {
	using namespace basic::thread_manager;

	core::Size const curthread( RosettaThreadManager::get_instance()->get_rosetta_thread_index() );

	TR_level3 << "Level 3 reporting in from thread " << curthread << ".  " << thread_assignments.get_requested_thread_count() << " threads were requested, and " << thread_assignments.get_assigned_total_thread_count() << " were assigned.  Computing (row, col)=(" << row << "," << column << ")." << std::endl;

	random_numbers[row][column] += numeric::random::random_range( 1, 10 );

	TR_level3 << "Level 2 thread " << curthread << " signing off." << std::endl;
}

/// @brief Function for level two threads: computing a row in the array.
void level2_function(
	utility::vector1< utility::vector1< core::Size > > & random_numbers,
	core::Size const row,
	core::Size const level_3_threads,
	bool const whole_row,
	basic::thread_manager::RosettaThreadAssignmentInfo const & thread_assignments
) {
	using namespace basic::thread_manager;

	if ( whole_row ) {
		runtime_assert( level_3_threads == 1 );
	}

	core::Size const curthread( RosettaThreadManager::get_instance()->get_rosetta_thread_index() );

	TR_level2 << "Level 2 reporting in from thread " << curthread << ".  " << thread_assignments.get_requested_thread_count() << " threads were requested, and " << thread_assignments.get_assigned_total_thread_count() << " were assigned.  Computing row " << row << "." << std::endl;

	if ( whole_row ) {
		for ( core::Size column(1); column <= 100; ++column ) {
			random_numbers[row][column] += numeric::random::random_range( 1, 10 );
		}
	} else {
		//Create the vector of work to do:
		utility::vector1< RosettaThreadFunction > level3_thread_functions;
		RosettaThreadAssignmentInfo thread_assignments3( RosettaThreadRequestOriginatingLevel::CORE_GENERIC );
		for ( core::Size i(0); i < 100; ++i ) {
			level3_thread_functions.push_back( std::bind( &level3_function, std::ref(random_numbers), row, i+1, std::cref(thread_assignments3) ) );
		}
		RosettaThreadManager::get_instance()->do_work_vector_in_threads( level3_thread_functions, level_3_threads, thread_assignments3 );
	}

	TR_level2 << "Level 2 thread " << curthread << " signing off." << std::endl;
}

/// @brief Function for level one threads: computing a block of rows in the array.
void level1_function(
	utility::vector1< utility::vector1< core::Size > > & random_numbers,
	core::Size const row_start,
	core::Size const row_end,
	core::Size const level_2_threads,
	core::Size const level_3_threads,
	bool const whole_row,
	basic::thread_manager::RosettaThreadAssignmentInfo const & thread_assignments
) {
	using namespace basic::thread_manager;
	debug_assert( row_end >= row_start );

	core::Size const curthread( RosettaThreadManager::get_instance()->get_rosetta_thread_index() );

	TR_level1 << "Level 1 reporting in from thread " << curthread << ".  " << thread_assignments.get_requested_thread_count() << " threads were requested, and " << thread_assignments.get_assigned_total_thread_count() << " were assigned.  Computing rows " << row_start << " through " << row_end << "." << std::endl;

	//Create the vector of work to do:
	utility::vector1< RosettaThreadFunction > level2_thread_functions;
	RosettaThreadAssignmentInfo thread_assignments2( RosettaThreadRequestOriginatingLevel::PROTOCOLS_GENERIC );
	for ( core::Size i(0), imax(row_end - row_start + 1); i < imax; ++i ) {
		level2_thread_functions.push_back( std::bind( &level2_function, std::ref(random_numbers), row_start + i, level_3_threads, whole_row, std::cref(thread_assignments2) ) );
	}
	RosettaThreadManager::get_instance()->do_work_vector_in_threads( level2_thread_functions, level_2_threads, thread_assignments2 );

	TR_level1 << "Level 1 thread " << curthread << " signing off." << std::endl;
}

#endif //MULTI_THREADED

/// @brief Entry point for program execution.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		TR_main << "Starting test_rosetta_thread_manager_basic_API." << std::endl;
		TR_main << "Pilot app created 31 May 2019 by Vikram K. Mulligan, Flatiron Institute." << std::endl;
		TR_main << "For questions, please e-mail vmulligan@flatironinstitute.org." << std::endl;

#ifdef MULTI_THREADED

		using namespace basic::thread_manager;

		runtime_assert_string_msg(
			!(option[whole_row_mode].user() && option[whole_row_mode].value() && option[level3_threads].user()),
			"Error in test_rosetta_thread_manager_basic_API integration test: both -whole_row_mode and -level3_threads options were passed.  Please use one or the other."
		);

		//Global options -- local copies.
		bool const whole_row( option[whole_row_mode ].value() );
		core::Size const l1_threads( option[level1_threads]() );
		core::Size const l2_threads( option[level2_threads]() );
		core::Size const l3_threads( option[level3_threads]() );

		//We're going to construct a 2D array and fill it with random numbers.
		//Level one jobs will do big blocks of the array, level two jobs will do individual rows, and level three jobs will do individual entries.
		//This will not be as efficient as it can be, but the goal is testing, not efficiency.
		utility::vector1< utility::vector1 < core::Size > > random_numbers( 1000, utility::vector1< core::Size >( 100, 0 ) );

		//Create the vector of work to do:
		utility::vector1< RosettaThreadFunction > level1_thread_functions;
		RosettaThreadAssignmentInfo thread_assignments( RosettaThreadRequestOriginatingLevel::APPLICATIONS_OR_APPLICATION_PROTOCOLS );
		for ( core::Size i(0); i<50; ++i ) {
			level1_thread_functions.push_back( std::bind( &level1_function, std::ref(random_numbers), i*20+1, (i+1)*20, l2_threads, l3_threads, whole_row, std::cref(thread_assignments) ) );
		}

		//Pass it to the RosettaThreadManager:
		RosettaThreadManager::get_instance()->do_work_vector_in_threads( level1_thread_functions, l1_threads, thread_assignments );

		bool all_in_correct_range(true);
		for ( core::Size i(1), imax( random_numbers.size() ); i<=imax; ++i ) {
			for ( core::Size j(1), jmax( random_numbers[i].size() ); j<=jmax; ++j ) {
				TR_output << random_numbers[i][j] << " ";
				if ( random_numbers[i][j] == 0 || random_numbers[i][j] > 10 ) all_in_correct_range = false;
			}
			TR_output << std::endl;
		}

		runtime_assert_string_msg( all_in_correct_range, "Error!  Some entries were zero or greater than ten." );

#else //not MULTI_THREADED
		utility_exit_with_message( "This application cannot be run in the single-threaded build of Rosetta.  Please use the \"-extras=cxx11thread\" build option when compiling." );
#endif //MULTI_THREADED

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	TR_main << "Finished test_rosetta_thread_manager_basic_API.  Exiting with status 0 (no errors)." << std::endl;

	return 0;
}
