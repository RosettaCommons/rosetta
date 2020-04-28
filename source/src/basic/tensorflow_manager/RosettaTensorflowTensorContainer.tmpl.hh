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

#include <cstring> //memcpy

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

/// @brief combine utility::vector1< RosettaTensorflowTensorContainer<T1> > into one RosettaTensorflowTensorContainer<T1>
/// @details requires all tensors in the vector have dimensions { 1, x, y, ... , z } where x, y, z are THE SAME for all tensors
template< typename T >
RosettaTensorflowTensorContainer< T >
RosettaTensorflowTensorContainer< T >::combine_tensors(
	utility::vector1< RosettaTensorflowTensorContainer< T > > const & tensor_vector
){
	runtime_assert( ! tensor_vector.empty() );

#ifndef NDEBUG
	for( platform::Size t = 1; t <= tensor_vector.size(); ++t ){
		debug_assert( tensor_vector[ t ].dimension( 1 ) == 1 );
	}
#endif

	utility::vector1< int64_t > new_dimensions;
	//We already checked the vector is not empty, so [1] is safe
	new_dimensions.resize( tensor_vector[ 1 ].n_dimensions() );

	new_dimensions[ 1 ] = tensor_vector.size();

	for( platform::Size dim = 2; dim <= new_dimensions.size(); ++dim ){
		new_dimensions[ dim ] = tensor_vector[ 1 ].dimension( dim );
#ifndef NDEBUG
		for( platform::Size t = 2; t <= tensor_vector.size(); ++t ){
			debug_assert( tensor_vector[ t ].dimension( dim ) == platform::Size( new_dimensions[ dim ] ) );
		}
#endif
	}

	platform::Size const n_elements_per_tensor = tensor_vector[ 1 ].num_tensor_elements();
	for( platform::Size t = 2; t <= tensor_vector.size(); ++t ){
		debug_assert( tensor_vector[ t ].num_tensor_elements() == n_elements_per_tensor );
	}

	RosettaTensorflowTensorContainer< T > combined_tensor( new_dimensions );
	//tensor_ and tensor_data_ are updated at this point

	//Okay time for some power user optimizations
	//Be careful here
	//Copy data into new container
	platform::Size const size_of_tensor = sizeof( T ) * n_elements_per_tensor;
	for( platform::Size t = 1; t <= tensor_vector.size(); ++t ){
		platform::Size const zero_indexed_offset_in_dest = (t-1) * n_elements_per_tensor;
		T * const destination = & ( combined_tensor( zero_indexed_offset_in_dest + 1 ) );

		T const * const source = tensor_vector[ t ].raw_tensor_data_ptr();
		std::memcpy( destination, source, size_of_tensor );
	}

#ifndef NDEBUG
	for( platform::Size t = 1; t <= tensor_vector.size(); ++t ){
		platform::Size const offset = (t-1) * n_elements_per_tensor;

		//We want to check that all of the values are what they shoudl be
		//This is really messy when the dimension count is unknown,
		// so we're just going to treat each tensor as if it were 1D
		//This should work as of March 2020
		for( platform::Size element = 1; element <= n_elements_per_tensor; ++element ){
			debug_assert( tensor_vector[ t ]( element ) == combined_tensor( offset + element ) );
		}
	}
#endif

	return combined_tensor;
}

template< typename T >
void
RosettaTensorflowTensorContainer< T >::split_combined_tensors(
	RosettaTensorflowTensorContainer< T > combined_tensors,
	utility::vector1< RosettaTensorflowTensorContainer<T> > & tensor_vector
){
	//This makes all of the same assumptions as combine_tensors()
	debug_assert( ! tensor_vector.empty() );
	platform::Size const n_elements_per_tensor =
		tensor_vector[ 1 ].num_tensor_elements();
	platform::Size const size_of_tensor = sizeof( T ) * n_elements_per_tensor;
	for( platform::Size t = 1; t <= tensor_vector.size(); ++t ){
		T * const destination = tensor_vector[ t ].raw_tensor_data_ptr();

		platform::Size const zero_indexed_offset_in_source = (t-1) * n_elements_per_tensor;
		T * const source = & ( combined_tensors( zero_indexed_offset_in_source + 1 ) );
		std::memcpy( destination, source, size_of_tensor );
	}

#ifndef NDEBUG
	for( platform::Size t = 1; t <= tensor_vector.size(); ++t ){
		platform::Size const offset = (t-1) * n_elements_per_tensor;

		//We want to check that all of the values are what they shoudl be
		//This is really messy when the dimension count is unknown,
		// so we're just going to treat each tensor as if it were 1D
		//This should work as of March 2020
		for( platform::Size element = 1; element <= n_elements_per_tensor; ++element ){
			debug_assert( tensor_vector[ t ]( element ) == combined_tensors( offset + element ) );
		}
	}
#endif

}

template < typename T >
void
RosettaTensorflowTensorContainer<T>::copy_data( T const * src ){
 	debug_assert( tensor_data_ != nullptr );
 	std::memcpy( tensor_data_, src, sizeof(T) * num_tensor_elements() );
}

} //tensorflow_manager
} //basic

#endif //USE_TENSORFLOW
#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_tmpl_hh
