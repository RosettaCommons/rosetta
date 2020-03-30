// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh
/// @brief A container class for Tensorflow's TF_Tensor objects, which manages access, creation,
/// and destruction in a way that prevents memory leaks.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_tmpl_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_tmpl_hh

#ifdef USE_TENSORFLOW

// Project headers:
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>
#include <utility/vector0.hh>


namespace basic {
namespace tensorflow_manager {

/// @brief Initialization constructor.
/// @details Initializes this to be a tensor with N dimensions, given by the entries in a vector.
/// @note This does not fill the tensor with anything.
template < typename T >
RosettaTensorflowTensorContainer<T>::RosettaTensorflowTensorContainer(
	utility::vector1< int64_t > const & dimensions
) :
	utility::VirtualBase()
{
	initialize( get_datatype(), dimensions );
}

/// @brief Initialization constructor.
/// @details Initializes this to be a tensor with N dimensions, given by the entries in a vector.
/// @note This fills the tensor with the value passed as the third parameter.
template < typename T >
RosettaTensorflowTensorContainer<T>::RosettaTensorflowTensorContainer(
	utility::vector1< int64_t > const & dimensions,
	T const & init_value
) :
	utility::VirtualBase()
{
	initialize( get_datatype(), dimensions, init_value );
}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
template < typename T >
RosettaTensorflowTensorContainer<T>::RosettaTensorflowTensorContainer(
	RosettaTensorflowTensorContainer const & src
) :
	utility::VirtualBase(src)
{
	if( src.tensor_ == nullptr ) {
		tensor_ = nullptr;
	} else {
		tensor_ = clone_tensor( src.tensor_ );
	}
	update_tensor_data_pointer();
}

/// @brief Destructor.
template < typename T >
RosettaTensorflowTensorContainer<T>::~RosettaTensorflowTensorContainer() {
	if( tensor_ != nullptr ) {
		TF_DeleteTensor( tensor_ );
		//tensor_ = nullptr;
		//tensor_data_ = nullptr;
	}
}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
template < typename T >
RosettaTensorflowTensorContainerOP< T >
RosettaTensorflowTensorContainer<T>::clone() const {
	return utility::pointer::make_shared< RosettaTensorflowTensorContainer< T > >( *this );
}

/// @brief Initialize this to be a tensor with N dimensions, given by the entries in a vector.
/// @details If already initialized, the old tensor is first deleted.
/// @note This does not fill the tensor with anything.
template < typename T >
void
RosettaTensorflowTensorContainer<T>::initialize(
	TF_DataType const datatype,
	utility::vector1< int64_t > const & dimensions
) {
	runtime_assert_string_msg(!dimensions.empty(), "Error in RosettaTensorflowTensorContainer<T>::initialize(): The dimensions array cannot be empty.");
	platform::Size numentries(1);
	for( platform::Size const dim : dimensions ) {
		numentries *= dim;
	}
	if( tensor_ != nullptr ){
		TF_DeleteTensor( tensor_ );
	}
	tensor_ = TF_AllocateTensor(
		datatype,
		dimensions.data(),
		static_cast<int>( dimensions.size() ),
		static_cast< size_t >( numentries ) * TF_DataTypeSize(datatype)
	);
	runtime_assert_string_msg( tensor_ != nullptr, "Error in RosettaTensorflowTensorContainer<T>::initialize(): Failed to allocate tensor!");
	update_tensor_data_pointer();
	debug_assert( tensor_data_ != nullptr );
}

/// @brief Initialize this to be a tensor with N dimensions, given by the entries in a vector.
/// @details If already initialized, the old tensor is first deleted.
/// @note This fills the tensor with the value passed as the third parameter.
template < typename T >
void
RosettaTensorflowTensorContainer<T>::initialize(
	TF_DataType const datatype,
	utility::vector1< int64_t > const & dimensions,
	T const & init_value
) {
	initialize( datatype, dimensions );
	//Fill the tensor with the initialization value:
	debug_assert( tensor_ != nullptr && tensor_data_ != nullptr );
	for( platform::Size i(0), imax( num_tensor_elements() ); i<imax; ++i ) {
		tensor_data_[i] = init_value;
	}
}

/// @brief Is the TF_Tensor pointer pointing to something?
/// @details Returns false if tensor_ == nullptr, true otherwise.
template < typename T >
bool
RosettaTensorflowTensorContainer<T>::contains_data() const {
	return tensor_ != nullptr;
}

/// @brief Get the number of dimensions of this tensor.
/// @details Returns 0 if tensor_ == nullptr.
template < typename T >
platform::Size
RosettaTensorflowTensorContainer<T>::n_dimensions() const {
	if( tensor_ == nullptr ) return 0;
	return TF_NumDims( tensor_ );
}

/// @brief Get the measure of this tensor along the Nth dimension.
/// @note Dimension index is 1-based.
template < typename T >
platform::Size
RosettaTensorflowTensorContainer<T>::dimension(
	platform::Size const dimension_index
) const {
	debug_assert( tensor_ != nullptr );
	debug_assert( TF_NumDims(tensor_) >= static_cast<int>(dimension_index) );
	return TF_Dim( tensor_, static_cast<int>( dimension_index ) - 1 );
}

/// @brief Get the number of elements in a Tensor.
/// @details This is the product of the Tensor's dimensions.
/// @note The TF_TensorElementCount function is not available in all
/// versions of the Tensorflow library.
template < typename T >
platform::Size
RosettaTensorflowTensorContainer<T>::num_tensor_elements() const {
	platform::Size product(1);
	for( platform::Size i(0), imax( TF_NumDims(tensor_) ); i<imax; ++i ) {
		product *= TF_Dim( tensor_, i );
	}
	return product;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////


/// @brief Given a TF_Tensor, make a copy and return a pointer to the copy.
/// @details Intended ONLY for use in the context of the copy constructor of the RosettaTensorFlowContainer.  This function
/// creates a new object and returns a RAW pointer to it, meaning that there's no lifetime management!
/// @note Returns nullptr if src is nullptr.
template < typename T >
TF_Tensor *
RosettaTensorflowTensorContainer<T>::clone_tensor(
	TF_Tensor * src
) const {
	if( src == nullptr ) return nullptr;

	platform::Size const ndim( TF_NumDims(src) );

	utility::vector0< int64_t > dims(ndim);
	platform::Size nentries(1); //Number of entries in the tensor, used below.
	for( platform::Size i(0); i<ndim; ++i ){
		dims[i] = TF_Dim( src, i );
		nentries *= dims[i];
	}

	//Allocate a new tensor equal in size to the old one:
	TF_Tensor * new_tensor( TF_AllocateTensor(
		TF_TensorType( src ),
		dims.data(),
		static_cast< int >( ndim ),
		TF_TensorByteSize(src)
	) );

	//Copy the tensor data to the new tensor:
	T const * srcdata( static_cast< T const * >( TF_TensorData(src) ) );
	T* newdata( static_cast< T* >( TF_TensorData(new_tensor) ) );
	for( platform::Size i(0); i<nentries; ++i ) {
		newdata[i] = srcdata[i];
	}

	return new_tensor;
}

/// @brief Update the pointer to the tensor data.
/// @details Sets this to nullptr if tensor_ == nullptr; otherwise,
/// calls TF_TensorData.
template < typename T >
void
RosettaTensorflowTensorContainer<T>::update_tensor_data_pointer() {
	if( tensor_ == nullptr ) {
		tensor_data_ = nullptr;
	}
	else {
		tensor_data_ = static_cast< T* >( TF_TensorData( tensor_ ) );
	}
}

} //tensorflow_manager
} //basic

#endif //USE_TENSORFLOW
#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_tmpl_hh
