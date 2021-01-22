// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_tensorflow_graph_convolutional_nn.cc
/// @brief A pilot app that tests the ability to load and evaluate a graph convolutional
/// neural network made with Tensorflow within Rosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/types.hh>

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// basic headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>

static basic::Tracer TR("apps.pilot.vmullig.test_tensorflow_graph_convolutional_nn");


/// @brief Indicate which commandline flags are relevant to this application.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// option.add_relevant( in::file::s );
	// option.add_relevant( in::file::l );

}

OPT_KEY( String, expected_results_file )

#ifdef USE_TENSORFLOW
/// @brief Load expected results.  Throws on failure (e.g. empty file or file that doesn't exist).
/// @details Expects a file consisting of a column vector of floats.
utility::vector1< core::Real >
load_expected_results(
	std::string const & expected_results_filename
) {
	std::string const errmsg( "Error in load_expected_results() function in test_tensorflow_graph_convolutional_nn application: " );
	utility::vector1< core::Real > returnvec;
	utility::vector1< std::string > const lines( utility::split_by_newlines( utility::file_contents( expected_results_filename ) ) );
	for( std::string const & line : lines ) {
		utility::vector1< std::string > linesplit( utility::split_whitespace( line ) );
		if( linesplit.empty() ) continue;
		runtime_assert_string_msg( linesplit.size() == 1, errmsg + "Could not parse line \"" + line + "\"." );
		std::stringstream ss( linesplit[1] );
		core::Real curfloat;
		ss >> curfloat;
		runtime_assert_string_msg( !(ss.bad() || ss.fail() ), errmsg + "Could not interpret \"" + linesplit[1] + "\" as a float." );
		returnvec.push_back(curfloat);
	}
	runtime_assert_string_msg( !returnvec.empty(), errmsg + "No floating-point values found in file \"" + expected_results_filename + "\"." );
	return returnvec;
}

/// @brief Provide the test values for X, the node weight tensor.
void
initialize_X_values(
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & X
) {
	// From Python:
	// 	X = constant(
	//     [[
	//         [0.2,  0.5],
	//         [0.8,  -0.5],
	//         [0.4,  0.2],
	//         [0.7,  -0.8],
	//         [-0.1,  0.5],
	//         [0.0,  0.7],
	//         [0.1,  0.8]
	//     ]],
	//     dtype=float32
	//	)
	X( 1, 1, 1 ) = 0.2;
	X( 1, 1, 2 ) = 0.5;
	X( 1, 2, 1 ) = 0.8;
	X( 1, 2, 2 ) = -0.5;
	X( 1, 3, 1 ) = 0.4;
	X( 1, 3, 2 ) = 0.2;
	X( 1, 4, 1 ) = 0.7;
	X( 1, 4, 2 ) = -0.8;
	X( 1, 5, 1 ) = -0.1;
	X( 1, 5, 2 ) = 0.5;
	X( 1, 6, 1 ) = 0.0;
	X( 1, 6, 2 ) = 0.7;
	X( 1, 7, 1 ) = 0.1;
	X( 1, 7, 2 ) = 0.8;
}

/// @brief Provide the adjacency matrix defining the graph.
/// @details The test graph has 7 nodes and 8 edges.
void
initialize_A_values(
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & A
) {
	// From Python:
	// A = constant(
	//     [[
	//         [0, 1, 1, 0, 0, 0, 0],
	//         [1, 0, 0, 1, 1, 0, 0],
	//         [1, 0, 0, 1, 0, 0, 0],
	//         [0, 1, 1, 0, 1, 1, 0],
	//         [0, 1, 0, 1, 0, 0, 0],
	//         [0, 0, 0, 1, 0, 0, 1],
	//         [0, 0, 0, 0, 0, 1, 0]
	//     ]],
	//     dtype=float32
	// )

	utility::vector1< std::pair< unsigned short, unsigned short > > const edges{
		std::make_pair( 1, 2 ),
		std::make_pair( 1, 3 ),
		std::make_pair( 2, 4 ),
		std::make_pair( 3, 4 ),
		std::make_pair( 2, 5 ),
		std::make_pair( 4, 5 ),
		std::make_pair( 4, 6 ),
		std::make_pair( 6, 7 )
	};
	for( std::pair< unsigned short, unsigned short > const & p : edges ) {
		A(1, p.first, p.second ) = 1.0;
		A(1, p.second, p.first ) = 1.0;
	}
}

/// @brief Prepare the data on the edges:
void
initialize_E_values(
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & E
) {
	// From Python:
	// E = constant(
	// 	[[
	// 		[ [0,0], [7,3], [6,1], [0,0], [0,0], [0,0], [0,0] ],
	// 		[ [7,3], [0,0], [0,0], [8,9], [5,5], [0,0], [0,0] ],
	// 		[ [6,1], [0,0], [0,0], [6,6], [0,0], [0,0], [0,0] ],
	// 		[ [0,0], [8,9], [6,6], [0,0], [8,3], [7,2], [0,0] ],
	// 		[ [0,0], [5,5], [0,0], [8,3], [0,0], [0,0], [0,0] ],
	// 		[ [0,0], [0,0], [0,0], [7,2], [0,0], [0,0], [1,9] ],
	// 		[ [0,0], [0,0], [0,0], [0,0], [0,0], [1,9], [0,0] ]
	// 	]],
	// 	dtype=float32
	// )
	utility::vector1< std::tuple< unsigned short, unsigned short, float, float > > const edgedata {
		std::make_tuple( 1, 2, 7.0, 3.0 ),
		std::make_tuple( 1, 3, 6.0, 1.0 ),
		std::make_tuple( 2, 4, 8.0, 9.0 ),
		std::make_tuple( 2, 5, 5.0, 5.0 ),
		std::make_tuple( 3, 4, 6.0, 6.0 ),
		std::make_tuple( 4, 5, 8.0, 3.0 ),
		std::make_tuple( 4, 6, 7.0, 2.0 ),
		std::make_tuple( 6, 7, 1.0, 9.0 )
	};
	for( auto const & tp : edgedata ) {
		E( 1, std::get<0>(tp), std::get<1>(tp), 1 ) = std::get<2>(tp); 
		E( 1, std::get<0>(tp), std::get<1>(tp), 2 ) = std::get<3>(tp);
		E( 1, std::get<1>(tp), std::get<0>(tp), 1 ) = std::get<2>(tp); 
		E( 1, std::get<1>(tp), std::get<0>(tp), 2 ) = std::get<3>(tp);
	}
}

/// @brief Compare the expected and actual results.
void
compare_results (
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > const & actual,
	utility::vector1< core::Real > const & expected
) {
	std::string const errmsg( "Error in compare_results() function in test_tensorflow_graph_convolutional_nn application: " );
	bool success(true);
	runtime_assert_string_msg( expected.size() == 7, errmsg + "The expected results vector should have 7 entries, but I count " + std::to_string( expected.size() ) + "!" );

	TR << "NODE\tEXPECTED\tACTUAL\tPASS?" << std::endl;
	for( core::Size i(1); i<=7; ++i ) {
		TR << i << "\t" << expected[i] << "\t" << actual(1, i, 1) << "\t";
		if( std::abs( expected[i] - static_cast< core::Real >( actual(1, i, 1) ) ) < 1.0e-3 ) {
			TR << "YES" << std::endl;
		} else {
			TR << "NO" << std::endl;
			success = false;
		}
	}

	if( !success ) utility_exit_with_message( errmsg + "Comparison between expected and actual results failed!" );
}

/// @brief Carry out the test of loading and evaluating a graph convolutional neural network.
void
do_test( std::string const & expected_results_filename ) {
	using namespace basic::tensorflow_manager;

	//First, load expected results.  Throws on failure.
	utility::vector1< core::Real > expected_results( load_expected_results( expected_results_filename ) );

	//Next, load the GCN:
	std::string const gcn_fname( "protocol_data/tensorflow_graphs/gcn_test_model/" );
	RosettaTensorflowSessionContainerCOP sess( RosettaTensorflowManager::get_instance()->get_session( gcn_fname, "serve" ) );

	//Prepare the inputs:
	utility::vector1< RosettaTensorflowTensorContainer< float > > input_vect {
		RosettaTensorflowTensorContainer< float >( utility::vector1< int64_t >{ 1, 7, 2 }, 0.0 ),
		RosettaTensorflowTensorContainer< float >( utility::vector1< int64_t >{ 1, 7, 7 }, 0.0 ),
		RosettaTensorflowTensorContainer< float >( utility::vector1< int64_t >{ 1, 7, 7, 2 }, 0.0 )
	};
	initialize_X_values( input_vect[1] );
	initialize_A_values( input_vect[2] );
	initialize_E_values( input_vect[3] );

	//Run the session and collect output:
	RosettaTensorflowTensorContainer< float > output_tensor( utility::vector1< int64_t >{1, 7, 1}, 0.0 );
	sess->run_multiinput_session(
		utility::vector1< std::string >{"serving_default_X_in", "serving_default_A_in", "serving_default_E_in"},
		"StatefulPartitionedCall",
		input_vect,
		output_tensor
	);

	//Compare the results and return pass/fail.
	compare_results( output_tensor, expected_results );
}
#endif

/// @brief Program entry point.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT( expected_results_file, "A file listing the expected values from evaluation of the GCN.", "" );

		devel::init( argc, argv );
		register_options();

#ifdef USE_TENSORFLOW
		TR << "Starting test_tensorflow_graph_convolutional_nn application." << std::endl;
		TR << "Pilot application created 15 July 2020 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." << std::endl;
		TR << "This application is used to test Rosetta's ability to load and evaluate graph convolutional neural networks.  If this application exits with status 0 (no errors), then the test was successful." << std::endl;

		runtime_assert_string_msg(
			option[expected_results_file].user() && !option[expected_results_file]().empty(),
			"Error in test_tensorflow_graph_convolutional_nn application: the -expected_results_file flag must be provided, indicating a file listing expected results."
		);
		std::string const options_file( option[expected_results_file]() );
		debug_assert(!options_file.empty());

		do_test( options_file );
#else
		utility_exit_with_message("The test_tensorflow_graph_convolutional_nn application cannot be run unless Rosetta is compiled with the \"extras=tensorflow\" or \"extras=tensorflow_gpu\" option.");
#endif

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		TR.Error << "The test_tensorflow_graph_convolutional_nn application FAILED!  Check the output log for error messages." << std::endl;
		TR.Error.flush();
		return -1;
	}

#ifdef USE_TENSORFLOW
	TR << "The test_tensorflow_graph_convolutional_nn application PASSED.  Graph convolutional neural networks can be loaded and evaluated in Rosetta!" << std::endl;
#endif

	return 0;
}
