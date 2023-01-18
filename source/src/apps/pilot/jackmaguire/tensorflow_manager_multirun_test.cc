// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jackmaguire/tensorflow_manager_multirun_test.cc
/// @brief  See main/tests/integration/tests/tensorflow_manager/multirun/README.md
/// @details This code is scrapped together from google searches, so it does not meet our coding convention and I have no plans to correct it.
/// @author Jack Maguire

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>



#include <iostream>

#ifdef USE_TENSORFLOW
#include <core/select/residue_selector/ResidueSelector.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#include <iomanip>
#endif


static basic::Tracer TR( "apps.pilot.jackmaguire.tensorflow_manager_multirun_test" );

#ifdef USE_TENSORFLOW

OPT_1GRP_KEY( String, select, path_to_model )

#include <tensorflow/c/c_api.h>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::select::residue_selector;
using namespace basic::tensorflow_manager;

void test(
	double const expected,
	double const observed
){
  bool const pass( std::abs( observed - expected ) < 0.01 );
  TR << std::setprecision(5);
  TR << "Output from model: " << observed << "\tExpected: " << expected << "\t" << (pass ? "PASSED" : "FAILED") << std::endl;
	runtime_assert( pass );
}

void run_with_manager(
	std::string const & input_operation,
	std::string const & output_operation
){
	RosettaTensorflowSessionContainerCOP TF_session =
		RosettaTensorflowManager::get_instance()->get_session( option[ select::path_to_model ](), "serve" );

	utility::vector1< RosettaTensorflowTensorContainer< float > > input_tensors( 3 );
	input_tensors[ 1 ].initialize( TF_FLOAT, { 1, 3 } );
	input_tensors[ 1 ].value( 1, 1 ) = -1.0;
	input_tensors[ 1 ].value( 1, 2 ) =  0.0;
	input_tensors[ 1 ].value( 1, 3 ) =  1.0;

	input_tensors[ 2 ].initialize( TF_FLOAT, { 1, 3 } );
	input_tensors[ 2 ].value( 1, 1 ) =  1.0;
	input_tensors[ 2 ].value( 1, 2 ) =  2.0;
	input_tensors[ 2 ].value( 1, 3 ) =  3.0;

	input_tensors[ 3 ].initialize( TF_FLOAT, { 1, 3 } );
	input_tensors[ 3 ].value( 1, 1 ) =  0.0;
	input_tensors[ 3 ].value( 1, 2 ) =  0.0;
	input_tensors[ 3 ].value( 1, 3 ) =  0.0;

	utility::vector1< RosettaTensorflowTensorContainer< float > > output_tensors( 3 );
	output_tensors[ 1 ].initialize( TF_FLOAT, { 1, 1 } ); //expecting 0.63166666
	output_tensors[ 2 ].initialize( TF_FLOAT, { 1, 1 } ); //expecting 0.44829783
	output_tensors[ 3 ].initialize( TF_FLOAT, { 1, 1 } ); //expecting 0.4794323

	TF_session->multirun_session( input_operation, output_operation, input_tensors, output_tensors );

	test( output_tensors[ 1 ].value( 1, 1 ), 0.63166666 );
	test( output_tensors[ 2 ].value( 1, 1 ), 0.44829783 );
	test( output_tensors[ 3 ].value( 1, 1 ), 0.4794323  );

}

int main( int argc, char* argv[] ) {

	try {

		NEW_OPT( select::path_to_model, "Path to your directory of saved models. Check the integration test at main/tests/integration/tests/tensorflow_simple_model_load_and_evaluate for an example", "./saved_models" );

		devel::init( argc, argv );

		if ( ! option[ select::path_to_model ].user() ) {
			utility_exit_with_message( "Please use the -select:path_to_model flag" );
		}


		run_with_manager( "in1", "out/Sigmoid" );

		TR << "Success!" << std::endl; //if we got this far, everything worked

    // Ensure that the citation is written:
    basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info();

		return 0;
	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}

#else

int main( int argc, char* argv[] ){

	try {

		devel::init( argc, argv );

		TR << "This app is not useful without extras=tensorflow" << std::endl;

		return 0;

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}

#endif
