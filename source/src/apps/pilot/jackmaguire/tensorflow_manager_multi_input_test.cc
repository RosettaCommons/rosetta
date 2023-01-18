// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jackmaguire/tensorflow_manager_multi_input_test.cc
/// @brief  See main/tests/integration/tests/tensorflow_manager/multi_input_test/README.md
/// @details This code is scrapped together from google searches, so it does not meet our coding convention and I have no plans to correct it.
/// @author Jack Maguire

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>



#include <iostream>

#ifdef USE_TENSORFLOW
#include <core/select/residue_selector/ResidueSelector.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#include <basic/citation_manager/CitationManager.hh>
#endif


static basic::Tracer TR( "apps.pilot.jackmaguire.tensorflow_manager_multi_input_test" );

#ifdef USE_TENSORFLOW

OPT_1GRP_KEY( String, select, path_to_model )

#include <tensorflow/c/c_api.h>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::select::residue_selector;
using namespace basic::tensorflow_manager;

constexpr core::Real EXPECTED_SCORE = 0.500397;

void test(
	double const expected,
	double const observed
){
  bool const pass( std::abs( observed - expected ) < 0.01 );
  TR << std::setprecision(5);
  TR << "Output from model: " << observed << "\tExpected: " << expected << "\t" << (pass ? "PASSED" : "FAILED") << std::endl;
	runtime_assert( pass );
}

TF_Tensor *FloatTensor(const int64_t *dims, int num_dims, const float *values) {
  int64_t num_values = 1;

  for (int i = 0; i < num_dims; ++i) {
    num_values *= dims[i];
  }

  TF_Tensor *t =
    TF_AllocateTensor(TF_FLOAT, dims, num_dims, sizeof(float) * num_values);

  memcpy(TF_TensorData(t), values, sizeof(float) * num_values);

  return t;
}

core::Real run_with_no_manager(
  std::string const & input_operation1,
  std::string const & input_operation2,
  std::string const & output_operation
){

  auto status = TF_NewStatus();
  auto graph = TF_NewGraph();
  auto sess_opts = TF_NewSessionOptions();

  constexpr char kSavedModelTagServe[] = "serve";
  const char *tags[] = {kSavedModelTagServe};

  auto session = TF_LoadSessionFromSavedModel(
    sess_opts, nullptr, option[ select::path_to_model ]().c_str(), tags, 1, graph, nullptr, status);

  if (TF_GetCode(status) != TF_OK) {
    utility_exit_with_message( std::string( "Unable to load session: " ) + std::string( TF_Message(status) ) );
  }

  TF_DeleteSessionOptions(sess_opts);

  const int ninputs = 2;
  const int noutputs = 1;

  std::unique_ptr<TF_Output[]> inputs(new TF_Output[ninputs]);
  std::unique_ptr<TF_Tensor *[]> input_values(new TF_Tensor *[ninputs]);
  std::unique_ptr<TF_Output[]> outputs(new TF_Output[noutputs]);
  std::unique_ptr<TF_Tensor *[]> output_values(new TF_Tensor *[noutputs]);

  const float input1_row[ 3 ] = { -1, 0, 1 };
  const int64_t input1_dims[2] = {1, 3};
  auto input1_tensor = FloatTensor(input1_dims, 2, input1_row);

  const float input2_row[ 2 ] = { -2, 2 };
  const int64_t input2_dims[2] = {1, 2};
  auto input2_tensor = FloatTensor(input2_dims, 2, input2_row);

  {//input1
    TF_Operation * const op = TF_GraphOperationByName( graph, input_operation1.c_str() );
    runtime_assert( op );
    inputs.get()[0] =
      TF_Output{op, 0};
  }

  {//input2
    TF_Operation * const op = TF_GraphOperationByName( graph, input_operation2.c_str() );
    runtime_assert( op );
    inputs.get()[1] =
      TF_Output{op, 0};
  }

  if (TF_GetCode(status) != TF_OK) {
    utility_exit_with_message( std::string( "Unable to fetch input operation: " ) + std::string( TF_Message( status ) ) );
  }

  input_values.get()[0] = input1_tensor;
  input_values.get()[1] = input2_tensor;

  const float output_row[1] = { 0 };
  const int64_t output_dims[2] = {1, 1};
  auto output_tensor = FloatTensor( output_dims, 2, output_row );

  {
    TF_Operation * const out_op = TF_GraphOperationByName( graph, output_operation.c_str() );
    runtime_assert( out_op );
    outputs.get()[0] =
      TF_Output{ out_op, 0 };
  }

  output_values.get()[0] = output_tensor;

  TF_SessionRun(session, nullptr, inputs.get(), input_values.get(), ninputs,
    outputs.get(), output_values.get(), noutputs, nullptr, 0,
    nullptr, status);

  if (TF_GetCode(status) != TF_OK) {
    utility_exit_with_message( std::string( "Unable to run session: " ) + std::string( TF_Message(status) ) );
  }

  float *values = static_cast<float *>(TF_TensorData(output_values.get()[0]));
  float const result = values[ 0 ];

  TF_DeleteTensor(input1_tensor);
  TF_DeleteTensor(input2_tensor);
  TF_DeleteTensor(output_tensor);

  return result;
}

core::Real run_with_manager(
  std::string const & input_operation1,
  std::string const & input_operation2,
  std::string const & output_operation
){
  RosettaTensorflowSessionContainerCOP TF_session =
    RosettaTensorflowManager::get_instance()->get_session( option[ select::path_to_model ](), "serve" );

  utility::vector1< std::string > const input_names({ input_operation1, input_operation2 });

  utility::vector1< RosettaTensorflowTensorContainer< float > > input_tensors( 2 );

	utility::vector1< int64_t > dim1( { 1, 3 } );
	runtime_assert( dim1.size() == 2 );
  input_tensors[ 1 ].initialize( TF_FLOAT, dim1 );
  input_tensors[ 1 ]( 1, 1 ) = -1.0;
  input_tensors[ 1 ]( 1, 2 ) =  0.0;
  input_tensors[ 1 ]( 1, 3 ) =  1.0;

	utility::vector1< int64_t > dim2( { 1, 2 } );
	runtime_assert( dim2.size() == 2 );
  input_tensors[ 2 ].initialize( TF_FLOAT, dim2 );
  input_tensors[ 2 ]( 1, 1 ) = -2.0;
  input_tensors[ 2 ]( 1, 2 ) =  2.0;

  RosettaTensorflowTensorContainer< float > output_tensor;
	utility::vector1< int64_t > dim3( { 1, 1 } );
  output_tensor.initialize( TF_FLOAT, dim3 );

  TF_session->run_multiinput_session( input_names, output_operation, input_tensors, output_tensor );

  float const result = output_tensor( 1, 1 );
  return result;
}

void test_multirun(
  std::string const & input_operation1,
  std::string const & input_operation2,
  std::string const & output_operation
){
  RosettaTensorflowSessionContainerCOP TF_session =
    RosettaTensorflowManager::get_instance()->get_session( option[ select::path_to_model ](), "serve" );

  utility::vector1< std::string > const input_names({ input_operation1, input_operation2 });
	debug_assert( input_names.size() == 2 );

	utility::vector1< utility::vector1< RosettaTensorflowTensorContainer< float > > > input_tensors;
	utility::vector1< RosettaTensorflowTensorContainer< float > > output_tensors( 5 );

	input_tensors.resize( 2 );
	input_tensors[ 1 ].resize( 5 );
	input_tensors[ 2 ].resize( 5 );

	{
		core::Size tensor = 1;
		utility::vector1< int64_t > dim1( { 1, 3 } );
		input_tensors[ 1 ][ tensor ].initialize( TF_FLOAT, dim1 );
		input_tensors[ 1 ][ tensor ]( 1, 1 ) = -1.0;
		input_tensors[ 1 ][ tensor ]( 1, 2 ) =  0.0;
		input_tensors[ 1 ][ tensor ]( 1, 3 ) =  1.0;

		utility::vector1< int64_t > dim2( { 1, 2 } );
		input_tensors[ 2 ][ tensor ].initialize( TF_FLOAT, dim2 );
		input_tensors[ 2 ][ tensor ]( 1, 1 ) = -2.0;
		input_tensors[ 2 ][ tensor ]( 1, 2 ) =  2.0;

		utility::vector1< int64_t > dim3( { 1, 1 } );
		output_tensors[ tensor ].initialize( TF_FLOAT, dim3 );
	}

	{
		core::Size tensor = 2;
		utility::vector1< int64_t > dim1( { 1, 3 } );
		input_tensors[ 1 ][ tensor ].initialize( TF_FLOAT, dim1 );
		input_tensors[ 1 ][ tensor ]( 1, 1 ) =  1.0;
		input_tensors[ 1 ][ tensor ]( 1, 2 ) =  2.0;
		input_tensors[ 1 ][ tensor ]( 1, 3 ) =  3.0;

		utility::vector1< int64_t > dim2( { 1, 2 } );
		input_tensors[ 2 ][ tensor ].initialize( TF_FLOAT, dim2 );
		input_tensors[ 2 ][ tensor ]( 1, 1 ) = -4.0;
		input_tensors[ 2 ][ tensor ]( 1, 2 ) =  2.0;

		utility::vector1< int64_t > dim3( { 1, 1 } );
		output_tensors[ tensor ].initialize( TF_FLOAT, dim3 );
	}

	{
		core::Size tensor = 3;
		utility::vector1< int64_t > dim1( { 1, 3 } );
		input_tensors[ 1 ][ tensor ].initialize( TF_FLOAT, dim1 );
		input_tensors[ 1 ][ tensor ]( 1, 1 ) = -1.0;
		input_tensors[ 1 ][ tensor ]( 1, 2 ) =  5.0;
		input_tensors[ 1 ][ tensor ]( 1, 3 ) =  1.0;

		utility::vector1< int64_t > dim2( { 1, 2 } );
		input_tensors[ 2 ][ tensor ].initialize( TF_FLOAT, dim2 );
		input_tensors[ 2 ][ tensor ]( 1, 1 ) = -5.0;
		input_tensors[ 2 ][ tensor ]( 1, 2 ) =  2.0;

		utility::vector1< int64_t > dim3( { 1, 1 } );
		output_tensors[ tensor ].initialize( TF_FLOAT, dim3 );
	}

	{
		core::Size tensor = 4;
		utility::vector1< int64_t > dim1( { 1, 3 } );
		input_tensors[ 1 ][ tensor ].initialize( TF_FLOAT, dim1 );
		input_tensors[ 1 ][ tensor ]( 1, 1 ) = -1.0;
		input_tensors[ 1 ][ tensor ]( 1, 2 ) =  0.0;
		input_tensors[ 1 ][ tensor ]( 1, 3 ) = -1.0;

		utility::vector1< int64_t > dim2( { 1, 2 } );
		input_tensors[ 2 ][ tensor ].initialize( TF_FLOAT, dim2 );
		input_tensors[ 2 ][ tensor ]( 1, 1 ) = -2.0;
		input_tensors[ 2 ][ tensor ]( 1, 2 ) = -2.0;

		utility::vector1< int64_t > dim3( { 1, 1 } );
		output_tensors[ tensor ].initialize( TF_FLOAT, dim3 );
	}

	{
		core::Size tensor = 5;
		utility::vector1< int64_t > dim1( { 1, 3 } );
		input_tensors[ 1 ][ tensor ].initialize( TF_FLOAT, dim1 );
		input_tensors[ 1 ][ tensor ]( 1, 1 ) = -11.0;
		input_tensors[ 1 ][ tensor ]( 1, 2 ) =   0.0;
		input_tensors[ 1 ][ tensor ]( 1, 3 ) =  11.0;

		utility::vector1< int64_t > dim2( { 1, 2 } );
		input_tensors[ 2 ][ tensor ].initialize( TF_FLOAT, dim2 );
		input_tensors[ 2 ][ tensor ]( 1, 1 ) = -22.0;
		input_tensors[ 2 ][ tensor ]( 1, 2 ) =  22.0;

		utility::vector1< int64_t > dim3( { 1, 1 } );
		output_tensors[ tensor ].initialize( TF_FLOAT, dim3 );
	}

  TF_session->multirun_multiinput_session( input_names, output_operation, input_tensors, output_tensors );

	test( 0.500397,   output_tensors[ 1 ]( 1, 1 ) );
	test( 0.40970826, output_tensors[ 2 ]( 1, 1 ) );
	test( 0.40430167, output_tensors[ 3 ]( 1, 1 ) );
	test( 0.3356494,  output_tensors[ 4 ]( 1, 1 ) );
	test( 0.71121025, output_tensors[ 5 ]( 1, 1 ) );
}

int main( int argc, char* argv[] ) {

  try {

    NEW_OPT( select::path_to_model, "Path to your directory of saved models. Check the integration test at main/tests/integration/tests/tensorflow_simple_model_load_and_evaluate for an example", "./saved_models" );

    devel::init( argc, argv );

    if ( ! option[ select::path_to_model ].user() ) {
      utility_exit_with_message( "Please use the -select:path_to_model flag" );
    }


    core::Real const control = run_with_no_manager( "in1", "in2", "dense3/Sigmoid" );
    if( std::abs( EXPECTED_SCORE - control ) > 0.01 ){
      utility_exit_with_message( "Error with the test! Getting the wrong answer even without the manager!" );
    } else {
      TR << "Control is a success!" << std::endl;
    }

    core::Real const manager_score = run_with_manager( "in1", "in2", "dense3/Sigmoid" );
    if( std::abs( EXPECTED_SCORE - manager_score ) > 0.01 ){
      utility_exit_with_message( "Error! Manager returned " + std::to_string( manager_score ) + " instead of " + std::to_string( EXPECTED_SCORE ) );
    } else {
      TR << "Manager is a success!" << std::endl;
    }

		test_multirun( "in1", "in2", "dense3/Sigmoid" );

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

		TR << "This app is not useful with out extras=tensorflow" << std::endl;

		return 0;

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}

#endif
