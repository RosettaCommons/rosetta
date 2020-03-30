// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh
/// @brief A container for Rosetta Tensorflow sessions, allowing sessions to be loaded once and stored in the global Tensorflow session manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_tmpl_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_tmpl_hh

#ifdef USE_TENSORFLOW

#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>


namespace basic {
namespace tensorflow_manager {

/// @brief Given an input tensor and a place to stick the outputs, run the Tensorflow session.
/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
/// @param[in]	input_name		The name of the input tensor in the Tensorflow session graph.
/// @param[in]	output_name		The name of the output tensor in the Tensorflow session graph.
/// @param[in]	input_tensor	The tensor of inputs.
/// @param[out]	ouput_tensor	The place to stick the outputs, overwritten by this operation.
/// @param[out]	runtime			The time for the actual Tensorflow session evaluation, in microseconds.  The contents of runtime will be overwritten by tihs operation.
/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information, do something
/// with the runtime output variable.
template< typename T1, typename T2 >
void
RosettaTensorflowSessionContainer::run_session(
	std::string const & input_name,
	std::string const & output_name,
	RosettaTensorflowTensorContainer<T1> const & input_tensor,
	RosettaTensorflowTensorContainer<T2> & output_tensor,
	std::chrono::duration< double, std::micro > & runtime
) const {
	//Debug-mode sanity checks:
	debug_assert( status_ != nullptr );
	debug_assert( graph_ != nullptr );
	debug_assert( session_ != nullptr );

	//Get the locations of the input to the graph and the output from the graph:
	TF_Operation* graph_op_in( TF_GraphOperationByName( graph_, input_name.c_str() ) );
	TF_Operation* graph_op_out( TF_GraphOperationByName( graph_, output_name.c_str() ) );
	debug_assert( graph_op_in != nullptr );
	debug_assert( graph_op_out != nullptr );
	TF_Output graph_inputs( TF_Output{ graph_op_in, 0 } );
	TF_Output graph_outputs( TF_Output{ graph_op_out, 0 } );

	debug_assert( input_tensor.tensor_ != nullptr );

	//Actually run the network given the input data.  (Note that the RosettaTensorflowSessionContainer is the sole class allowed to access the private
	//member data of the RosettaTensorflowTensorContainer class.):
	std::chrono::time_point<std::chrono::system_clock> const starttime( ROSETTA_TENSORFLOW_CLOCK::now() );
	TF_SessionRun( session_, nullptr, &graph_inputs, &(input_tensor.tensor_), 1, &graph_outputs, &(output_tensor.tensor_), 1, nullptr, 0, nullptr, status_);
	std::chrono::time_point<std::chrono::system_clock> const endtime( ROSETTA_TENSORFLOW_CLOCK::now() );
	runtime = endtime - starttime;

	//Check that the run was okay:
	if ( TF_GetCode(status_) != TF_OK ) {
		utility_exit_with_message( "Unable to evaluate TensorFlow session: " + std::string( TF_Message(status_) ) );
	}

	//Ensure that the output tensor is in a good state for reading out the output:
	output_tensor.update_tensor_data_pointer();
}

} //namespace tensorflow_manager
} //namespace basic

#endif //USE_TENSORFLOW
#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_tmpl_hh
