// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jackmaguire/tensorflow_test2.cc
/// @brief  See main/tests/integration/tests/tensorflow_simple_model_load_and_evaluate/README.md
/// @details This code is scrapped together from google searches, so it does not meet our coding convention and I have no plans to correct it.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>



#include <iostream>

#ifdef USE_TENSORFLOW
#include <core/select/residue_selector/ResidueSelector.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <string>
#endif

static basic::Tracer TR( "apps.pilot.jackmaguire.tensorflow_test2" );

#ifdef USE_TENSORFLOW

OPT_1GRP_KEY( String, select, path_to_model )

#include <tensorflow/c/c_api.h>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::select::residue_selector;


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

void run( std::string const & input_operation, std::string const & output_operation ){

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

		const int ninputs = 1;
		const int noutputs = 1;

		std::unique_ptr<TF_Output[]> inputs(new TF_Output[ninputs]);
		std::unique_ptr<TF_Tensor *[]> input_values(new TF_Tensor *[ninputs]);
		std::unique_ptr<TF_Output[]> outputs(new TF_Output[noutputs]);
		std::unique_ptr<TF_Tensor *[]> output_values(new TF_Tensor *[noutputs]);

		const float input_row[ 10 ] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 0 };
		const int64_t input_dims[2] = {1, 10};
		auto input_tensor = FloatTensor(input_dims, 2, input_row);

		{//AAA
			size_t pos = 0;
			TF_Operation* oper;
			TR << "All Operations: " << std::endl;
			while ((oper = TF_GraphNextOperation(graph, &pos)) != nullptr) {
				TR << TF_OperationName( oper ) << std::endl;
			}
		}

		{
			TF_Operation * const firstop = TF_GraphOperationByName( graph, input_operation.c_str() );
			runtime_assert( firstop );
			inputs.get()[0] =
				TF_Output{firstop, 0};
		}

		if (TF_GetCode(status) != TF_OK) {
			utility_exit_with_message( std::string( "Unable to fetch input operation: " ) + std::string( TF_Message( status ) ) );
		}

		input_values.get()[0] = input_tensor;

		const float output_row[2] = { 0, 0 };
		const int64_t output_dims[2] = {1, 2};
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
		TR << "Values are " << values[ 0 ] << " and " << values[ 1 ] << ". We expect -0.40661177 and -0.5512837." << std::endl;
		//hoping for [[-0.40661177 -0.5512837 ]]
		if( values[ 0 ] < -0.4067 || values[ 0 ] > -0.4066 ){
			utility_exit_with_message( "First value should be close to -0.40661177 but is actually " + std::to_string( values[ 0 ] ) );
		}
		if( values[ 1 ] < -0.5513 || values[ 1 ] > -0.5512 ){
			utility_exit_with_message( "Second value should be close to -0.5512837 but is actually " + std::to_string( values[ 1 ] ) );
		}

		TF_DeleteTensor(input_tensor);
		TF_DeleteTensor(output_tensor);
}

int main( int argc, char* argv[] ) {

	try {

		NEW_OPT( select::path_to_model, "Path to your directory of saved models. Check the integration test at main/tests/integration/tests/tensorflow_simple_model_load_and_evaluate for an example", "./saved_models" );

		devel::init( argc, argv );

		if ( ! option[ select::path_to_model ].user() ) {
			utility_exit_with_message( "Please use the -select:path_to_model flag" );
		}

		//This was generated by the "AAA" scope in run
		utility::vector1< std::string > options = {
			"dense_input",//THIS IS INPUT (1)
			"dense/kernel/Initializer/random_uniform/shape",
			"dense/kernel/Initializer/random_uniform/min",
			"dense/kernel/Initializer/random_uniform/max",
			"dense/kernel/Initializer/random_uniform/RandomUniform",
			"dense/kernel/Initializer/random_uniform/sub",
			"dense/kernel/Initializer/random_uniform/mul",
			"dense/kernel/Initializer/random_uniform",
			"dense/kernel",
			"dense/kernel/IsInitialized/VarIsInitializedOp",
			"dense/kernel/Assign",
			"dense/kernel/Read/ReadVariableOp",
			"dense/bias/Initializer/zeros",
			"dense/bias",
			"dense/bias/IsInitialized/VarIsInitializedOp",
			"dense/bias/Assign",
			"dense/bias/Read/ReadVariableOp",
			"dense/MatMul/ReadVariableOp",
			"dense/MatMul",
			"dense/BiasAdd/ReadVariableOp",
			"dense/BiasAdd",
			"dense/Relu",
			"dense_1/kernel/Initializer/random_uniform/shape",
			"dense_1/kernel/Initializer/random_uniform/min",
			"dense_1/kernel/Initializer/random_uniform/max",
			"dense_1/kernel/Initializer/random_uniform/RandomUniform",
			"dense_1/kernel/Initializer/random_uniform/sub",
			"dense_1/kernel/Initializer/random_uniform/mul",
			"dense_1/kernel/Initializer/random_uniform",
			"dense_1/kernel",
			"dense_1/kernel/IsInitialized/VarIsInitializedOp",
			"dense_1/kernel/Assign",
			"dense_1/kernel/Read/ReadVariableOp",
			"dense_1/bias/Initializer/zeros",
			"dense_1/bias",
			"dense_1/bias/IsInitialized/VarIsInitializedOp",
			"dense_1/bias/Assign",
			"dense_1/bias/Read/ReadVariableOp",
			"dense_1/MatMul/ReadVariableOp",
			"dense_1/MatMul",//This is output? (40)
			"dense_1/BiasAdd/ReadVariableOp",
			"dense_1/BiasAdd",//This is the other option for output (42)
			"predict/group_deps",
			"Const",
			"Const_1",
			"Const_2",
			"Const_3",
			"Const_4",
			"Const_5",
			"Const_6",
			"Const_7",
			"Const_8",
			"RestoreV2/tensor_names",
			"RestoreV2/shape_and_slices",
			"RestoreV2",
			"Identity",
			"AssignVariableOp",
			"RestoreV2_1/tensor_names",
			"RestoreV2_1/shape_and_slices",
			"RestoreV2_1",
			"Identity_1",
			"AssignVariableOp_1",
			"RestoreV2_2/tensor_names",
			"RestoreV2_2/shape_and_slices",
			"RestoreV2_2",
			"Identity_2",
			"AssignVariableOp_2",
			"RestoreV2_3/tensor_names",
			"RestoreV2_3/shape_and_slices",
			"RestoreV2_3",
			"Identity_3",
			"AssignVariableOp_3",
			"VarIsInitializedOp",
			"VarIsInitializedOp_1",
			"VarIsInitializedOp_2",
			"VarIsInitializedOp_3",
			"init",
			"save/filename/input",
			"save/filename",
			"save/Const",
			"save/Const_1",
			"save/Const_2",
			"save/Const_3",
			"save/Const_4",
			"save/SaveV2/tensor_names",
			"save/SaveV2/shape_and_slices",
			"save/SaveV2",
			"save/control_dependency",
			"save/RestoreV2/tensor_names",
			"save/RestoreV2/shape_and_slices",
			"save/RestoreV2",
			"save/NoOp",
			"save/NoOp_1",
			"save/NoOp_2",
			"save/Identity",
			"save/AssignVariableOp",
			"save/Identity_1",
			"save/AssignVariableOp_1",
			"save/NoOp_3",
			"save/Identity_2",
			"save/AssignVariableOp_2",
			"save/Identity_3",
			"save/AssignVariableOp_3",
			"save/restore_all",
			"init_1"
		};

		run( options[ 1 ], options[ 40 ] );
		TR << "Success!" << std::endl; //if we got this far, everything worked
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
