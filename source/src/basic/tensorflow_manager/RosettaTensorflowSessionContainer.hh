// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh
/// @brief A container for Rosetta Tensorflow sessions, allowing sessions to be loaded once and stored in the global Tensorflow session manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_hh

#ifdef USE_TENSORFLOW

#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

//Tensorflow headers
#include <tensorflow/c/c_api.h>

//C++ headers:
#include <chrono>

namespace basic {
namespace tensorflow_manager {

/// @brief A container for Rosetta Tensorflow sessions, allowing sessions to be loaded once and stored in the global Tensorflow session manager.
/// @note Cannot be subclassed.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class RosettaTensorflowSessionContainer final : public utility::VirtualBase {

	friend class RosettaTensorflowManager;

private:

	/// @brief Default constructor -- explicitly deleted.
	RosettaTensorflowSessionContainer() = delete;

	/// @brief Initialization constructor.
	/// @details Takes name of stored session file.  Triggers read from disk.  Only accessible to TensorflowManager.
	/// @note If session options are passed in, then the RosettaTensorflowSessionContainer takes control of them and manages their
	/// lifetime and destruction.  They must be allocated on the heap with `new`!
	RosettaTensorflowSessionContainer( std::string const & filename, std::string const & tags, TF_SessionOptions* sess_options=nullptr );

public:

	/// @brief Copy constructor -- explicitly deleted.
	RosettaTensorflowSessionContainer(RosettaTensorflowSessionContainer const & ) = delete;

	/// @brief Assignment operator -- explicitly deleted.
	RosettaTensorflowSessionContainer operator=(RosettaTensorflowSessionContainer const & ) = delete;

	/// @brief Destructor.
	/// @details It is necessary to clean up Tensorflow objects here.
	~RosettaTensorflowSessionContainer() override;

	/// @brief Was this session initialized with custom Tensorflow session options?
	inline bool was_initialized_with_custom_options() const { return has_custom_options_; }

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
	void run_session(
		std::string const & input_name,
		std::string const & output_name,
		RosettaTensorflowTensorContainer<T1> const & input_tensor,
		RosettaTensorflowTensorContainer<T2> & output_tensor,
		std::chrono::duration< double, std::micro > & runtime
	) const;

	/// @brief Given a vector of input tensors and vector of outputs, run the Tensorflow session.
	/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
	/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
	/// @param[in]	input_name				The name of the input tensor in the Tensorflow session graph.
	/// @param[in]	output_name				The name of the output tensor in the Tensorflow session graph.
	/// @param[in]	input_tensor_vector		The vector of tensor of inputs.
	/// @param[out]	output_tensor_vector	The vector of places to stick the outputs, overwritten by this operation.
	/// @param[out]	runtime					The time for the actual Tensorflow session evaluation, in microseconds.  The contents of runtime will be overwritten by tihs operation.
	/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information, do something
	/// with the runtime output variable.  Also note that the tensor types must be the same in the input and output vectors.
	template< typename T1, typename T2 >
	void multirun_session(
		std::string const & input_name,
		std::string const & output_name,
		utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensor_vector,
		utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensor_vector,
		std::chrono::duration< double, std::micro > & runtime
	) const;

	/// @brief Overload of multirun_session for those who don't need runtime.
	/// @author Jack Maguire, Menten AI (jack@menten.ai).
	template< typename T1, typename T2 >
	void multirun_session(
		std::string const & input_name,
		std::string const & output_name,
		utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensor_vector,
		utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensor_vector
	) const {
		std::chrono::duration< double, std::micro > runtime;
		multirun_session( input_name, output_name, input_tensor_vector, output_tensor_vector, runtime );
	}

	/// @brief Equivalent of run_session() but for networks with multiple heads (multiple input tensors per run)
 	/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
 	/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
 	/// @param[in]	input_names		The name of the input tensors in the Tensorflow session graph.
 	/// @param[in]	output_name		The name of the output tensor in the Tensorflow session graph.
 	/// @param[in]	input_tensor	The tensors of inputs, must be in the same order as input_names.
 	/// @param[out]	ouput_tensor	The place to stick the outputs, overwritten by this operation.
 	/// @param[out]	runtime			The time for the actual Tensorflow session evaluation, in microseconds.  The contents of runtime will be overwritten by tihs operation.
 	/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information, do something
 	/// with the runtime output variable.
	/// @author Jack Maguire, Menten AI (jack@menten.ai).
 	template< typename T1, typename T2 >
 	void run_multiinput_session(
 		utility::vector1< std::string > const & input_names,
 		std::string const & output_name,
 		utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensors,
 		RosettaTensorflowTensorContainer<T2> & output_tensor,
 		std::chrono::duration< double, std::micro > & runtime
 	) const;

	/// @brief Overload of run_multiinput_session for those who don't need runtime.
	/// @author Jack Maguire, Menten AI (jack@menten.ai).
 	template< typename T1, typename T2 >
 	void run_multiinput_session(
 		utility::vector1< std::string > const & input_names,
 		std::string const & output_name,
 		utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensors,
 		RosettaTensorflowTensorContainer<T2> & output_tensor
 	) const {
 		std::chrono::duration< double, std::micro > runtime;
 		run_multiinput_session( input_names, output_name, input_tensors, output_tensor, runtime );
 	}


	/// @brief Run multiple passes of a network with multiple inputs
	/// @details for input_tensors, the INNER vector groups the different heads, the OUTER vector
	/// groups the different runs. input_tensors[x].size() == output_tensors.size(), if that helps.
	/// @author Jack Maguire, Menten AI (jack@menten.ai).
 	template< typename T1, typename T2 >
 	void multirun_multiinput_session(
 		utility::vector1< std::string > const & input_names,
 		std::string const & output_name,
 		utility::vector1< utility::vector1< RosettaTensorflowTensorContainer<T1> > > const & input_tensors,
 		utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensors,
 		std::chrono::duration< double, std::micro > & runtime
 	) const;

	/// @brief Overload of multirun_multiinput_session for those who don't need runtime.
	/// @author Jack Maguire, Menten AI (jack@menten.ai).
 	template< typename T1, typename T2 >
 	void multirun_multiinput_session(
 		utility::vector1< std::string > const & input_names,
 		std::string const & output_name,
 		utility::vector1< utility::vector1< RosettaTensorflowTensorContainer<T1> > > const & input_tensors,
 		utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensors
 	) const {
 		std::chrono::duration< double, std::micro > runtime;
 		multirun_multiinput_session( input_names, output_name, input_tensors, output_tensors, runtime );
 	}

	/// @brief List all of the operations available in a loaded model.
	void list_all_operations( std::ostream & outstream ) const;

private:

	/// @brief Status of the Tensorflow session.
	TF_Status* status_ = nullptr;

	/// @brief Graph for the Tensorflow session.
	TF_Graph* graph_ = nullptr;

	/// @brief Options for the Tensorflow session.
	TF_SessionOptions* sess_options_ = nullptr;

	/// @brief Was this initialized with custom TF_SessionOptions?
	bool has_custom_options_ = false;

	/// @brief The Tensorflow session.
	/// @details Note that this is stored by raw pointer.  The destructor for this container class destroys the session using
	/// Tensorflow's safe cleanup functions.
	TF_Session* session_ = nullptr;

};

} //tensorflow_manager
} //basic

#endif //USE_TENSORFLOW
#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_hh
