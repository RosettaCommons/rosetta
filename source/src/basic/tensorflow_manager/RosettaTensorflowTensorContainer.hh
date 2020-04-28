// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh
/// @brief A container class for Tensorflow's TF_Tensor objects, which manages access, creation,
/// and destruction in a way that prevents memory leaks.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_hh

#ifdef USE_TENSORFLOW

#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>
#include <basic/tensorflow_manager/TFDataTypeDetector.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// Platform headers
#include <platform/types.hh>

// External headers
#include <tensorflow/c/c_api.h>

namespace basic {
namespace tensorflow_manager {

/// @brief A container class for Tensorflow's TF_Tensor objects, which manages access, creation,
/// and destruction in a way that prevents memory leaks.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template < typename T >
class RosettaTensorflowTensorContainer : public utility::VirtualBase {

	//Needed to allow the TF_Session to access the TF_Tensor directly:
	friend class RosettaTensorflowSessionContainer;
	//Note that this is a container with only one thing that it manages, and the only thing
	//that it is supposed to do is keep the tensor that it manages private from everything
	//EXCEPT the RosettaTensorflowSessionContainer.

public:

	/// @brief Default constructor.
	/// @details Creates a nullptr pointer to the tensor object.  Tensor must
	/// be initialized before using it!
	RosettaTensorflowTensorContainer() = default;

	/// @brief Initialization constructor.
	/// @details Initializes this to be a tensor with N dimensions, given by the entries in a vector.
	/// @note This does not fill the tensor with anything.
	RosettaTensorflowTensorContainer( utility::vector1< int64_t > const & dimensions );

	/// @brief Initialization constructor.
	/// @details Initializes this to be a tensor with N dimensions, given by the entries in a vector.
	/// @note This fills the tensor with the value passed as the third parameter.
	RosettaTensorflowTensorContainer( utility::vector1< int64_t > const & dimensions, T const & init_value );

	/// @brief Copy constructor.
	RosettaTensorflowTensorContainer(RosettaTensorflowTensorContainer const & src);

	/// @brief Destructor.
	~RosettaTensorflowTensorContainer() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	RosettaTensorflowTensorContainerOP< T > clone() const;

	/// @brief Delete the assignment operator.
	RosettaTensorflowTensorContainer & operator=( RosettaTensorflowTensorContainer const & ) = delete;

public:

	/// @brief Initialize this to be a tensor with N dimensions, given by the entries in a vector.
	/// @details If already initialized, the old tensor is first deleted.
	/// @note This does not fill the tensor with anything.
	void initialize( TF_DataType const datatype, utility::vector1< int64_t > const & dimensions );

	/// @brief Initialize this to be a tensor with N dimensions, given by the entries in a vector.
	/// @details If already initialized, the old tensor is first deleted.
	/// @note This fills the tensor with the value passed as the third parameter.
	void initialize( TF_DataType const datatype, utility::vector1< int64_t > const & dimensions, T const & init_value  );

	/// @brief Is the TF_Tensor pointer pointing to something?
	/// @details Returns false if tensor_ == nullptr, true otherwise.
	bool contains_data() const;

	/// @brief Get the number of dimensions of this tensor.
	/// @details Returns 0 if tensor_ == nullptr.
	platform::Size n_dimensions() const;

	/// @brief Get the measure of this tensor along the Nth dimension.
	/// @note Dimension index is 1-based.
	platform::Size dimension( platform::Size const dimension_index ) const;

	/// @brief Get the number of elements in this Tensor.
	/// @details This is the product of the Tensor's dimensions.
	/// @note The TF_TensorElementCount function is not available in all
	/// versions of the Tensorflow library.
	platform::Size num_tensor_elements() const;

public: //Getters

	/// @brief Access an entry.  (Const.)
	/// @details Works with any tensor dimension.
	/// @note Indices are 1-based.
	inline
	T const &
	operator()(
		platform::Size const coord
	) const {
		debug_assert( tensor_ != nullptr );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord > 0  && coord <= static_cast<platform::Size>( num_tensor_elements() ) );
		return tensor_data_[coord - 1];
	}

	/// @brief Access an entry.  (Const.)
	/// @details Specific for 2D tensors (matrices).
	/// @note Indices are 1-based.
	inline
	T const &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2
	) const {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 2 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 2) + (coord2 - 1) ];
	}

	/// @brief Access an entry.  (Const.)
	/// @details Specific for 3D tensors.
	/// @note Indices are 1-based.
	inline
	T const &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2,
		platform::Size const coord3
	) const {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 3 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		debug_assert( coord2 > 0 && coord3 <= static_cast<platform::Size>( TF_Dim(tensor_, 3) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3) + (coord2 - 1)*TF_Dim(tensor_, 3) + (coord3 - 1) ];
	}

	/// @brief Access an entry.  (Const.)
	/// @details Specific for 4D tensors.
	/// @note Indices are 1-based.
	inline
	T const &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2,
		platform::Size const coord3,
		platform::Size const coord4
	) const {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 4 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		debug_assert( coord2 > 0 && coord3 <= static_cast<platform::Size>( TF_Dim(tensor_, 3) ) );
		debug_assert( coord2 > 0 && coord4 <= static_cast<platform::Size>( TF_Dim(tensor_, 4) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4) + (coord2 - 1)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4) + (coord3 - 1)*TF_Dim(tensor_, 4) + (coord4 - 1) ];
	}

	/// @brief Access an entry.  (Const.)
	/// @details Specific for 5D tensors.
	/// @note Indices are 1-based.
	inline
	T const &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2,
		platform::Size const coord3,
		platform::Size const coord4,
		platform::Size const coord5
	) const {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 5 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		debug_assert( coord2 > 0 && coord3 <= static_cast<platform::Size>( TF_Dim(tensor_, 3) ) );
		debug_assert( coord2 > 0 && coord4 <= static_cast<platform::Size>( TF_Dim(tensor_, 4) ) );
		debug_assert( coord2 > 0 && coord5 <= static_cast<platform::Size>( TF_Dim(tensor_, 5) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4)*TF_Dim(tensor_, 5) + (coord2 - 1)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4)*TF_Dim(tensor_, 5) + (coord3 - 1)*TF_Dim(tensor_, 4)*TF_Dim(tensor_, 5) + (coord4 - 1)*TF_Dim(tensor_, 5) + (coord5-1) ];
	}

	//More verbose option
	template< typename... Args>
	T const &
	value( Args&&... args ) const {
		return operator()( std::forward<Args>(args)... );
	}


public: //Setters

	/// @brief Access an entry.  (Non-const.)
	/// @details Works with any tensor dimension.
	/// @note Indices are 1-based.
	inline
	T &
	operator()(
		platform::Size const coord
	) {
		debug_assert( tensor_ != nullptr );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord > 0  && coord <= static_cast<platform::Size>( num_tensor_elements() ) );
		return tensor_data_[coord - 1];
	}

	/// @brief Access an entry.  (Non-const.)
	/// @details Specific for 2D tensors (matrices).
	/// @note Indices are 1-based.
	inline
	T &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2
	) {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 2 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 0) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 1) + (coord2 - 1) ];
	}

	/// @brief Access an entry.  (Non-const.)
	/// @details Specific for 3D tensors.
	/// @note Indices are 1-based.
	inline
	T &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2,
		platform::Size const coord3
	) {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 3 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 0) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord3 > 0 && coord3 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 1)*TF_Dim(tensor_, 2) + (coord2 - 1)*TF_Dim(tensor_, 2) + (coord3 - 1) ];
	}

	/// @brief Access an entry.  (Non-const.)
	/// @details Specific for 4D tensors.
	/// @note Indices are 1-based.
	inline
	T &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2,
		platform::Size const coord3,
		platform::Size const coord4
	) {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 4 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 0) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord3 > 0 && coord3 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		debug_assert( coord4 > 0 && coord4 <= static_cast<platform::Size>( TF_Dim(tensor_, 3) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3) + (coord2 - 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3) + (coord3 - 1)*TF_Dim(tensor_, 3) + (coord4 - 1) ];
	}

	/// @brief Access an entry.  (Non-const.)
	/// @details Specific for 5D tensors.
	/// @note Indices are 1-based.
	inline
	T &
	operator()(
		platform::Size const coord1,
		platform::Size const coord2,
		platform::Size const coord3,
		platform::Size const coord4,
		platform::Size const coord5
	) {
		debug_assert( tensor_ != nullptr );
		debug_assert( TF_NumDims(tensor_) == 5 );
		debug_assert( tensor_data_ != nullptr );
		debug_assert( coord1 > 0 && coord1 <= static_cast<platform::Size>( TF_Dim(tensor_, 0) ) );
		debug_assert( coord2 > 0 && coord2 <= static_cast<platform::Size>( TF_Dim(tensor_, 1) ) );
		debug_assert( coord3 > 0 && coord3 <= static_cast<platform::Size>( TF_Dim(tensor_, 2) ) );
		debug_assert( coord4 > 0 && coord4 <= static_cast<platform::Size>( TF_Dim(tensor_, 3) ) );
		debug_assert( coord5 > 0 && coord5 <= static_cast<platform::Size>( TF_Dim(tensor_, 4) ) );
		return tensor_data_[ (coord1 - 1)*TF_Dim(tensor_, 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4) + (coord2 - 1)*TF_Dim(tensor_, 2)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4) + (coord3 - 1)*TF_Dim(tensor_, 3)*TF_Dim(tensor_, 4) + (coord4 - 1)*TF_Dim(tensor_, 4) + (coord5-1) ];
	}

	//More verbose option
	template< typename... Args>
	T &
	value( Args&&... args ) {
		return operator()( std::forward<Args>(args)... );
	}

protected: //Danger Zone!
	TF_Tensor const * raw_tensor_ptr() const {
		return tensor_;
	}

	TF_Tensor * raw_tensor_ptr() {
		return tensor_;
	}

	T const * raw_tensor_data_ptr() const {
		return tensor_data_;
	}

	T * raw_tensor_data_ptr() {
		return tensor_data_;
	}

	/// @brief Get the Tensorflow data type given the C++ data type.
	static TF_DataType
	get_datatype() {
		TFDataTypeDetector< T >::validate();
		TFDataTypeDetector< T > t;
		return t.value;
	}

	/// @brief Given a TF_Tensor, make a copy and return a pointer to the copy.
	/// @details Intended ONLY for use in the context of the copy constructor of the RosettaTensorFlowContainer.  This function
	/// creates a new object and returns a RAW pointer to it, meaning that there's no lifetime management!
	/// @note Returns nullptr if src is nullptr.
	TF_Tensor * clone_tensor( TF_Tensor * src ) const;

	/// @brief Update the pointer to the tensor data.
	/// @details Sets this to nullptr if tensor_ == nullptr; otherwise,
	/// calls TF_TensorData.
	void update_tensor_data_pointer();

public: //static utilities

	static
	RosettaTensorflowTensorContainer< T >
	combine_tensors(
		utility::vector1< RosettaTensorflowTensorContainer< T > > const & tensor_vector
	);

	static
	void
	split_combined_tensors(
		RosettaTensorflowTensorContainer< T > combined_tensors,
		utility::vector1< RosettaTensorflowTensorContainer<T> > & tensor_vector
	);

	template< typename Container >
	void
	copy_data_from_container( Container const & src ){
		debug_assert( src.size() == num_tensor_elements() );
		copy_data( src.data() );
	}

protected:

	///@brief memcpy data from src to tensor_data_
 	void
 	copy_data( T const * src );

private:

	/// @brief A pointer to the Tensorflow tensor whose life we're managing.
	/// @details Can be nullptr if no tensor has been created yet.
	TF_Tensor * tensor_ = nullptr;

	/// @brief A pointer to the data in the tensor whose life we're managing.
	/// @details Will be nullptr if tensor_ is nullptr.  Updated to point to
	/// the data by any function that creates a tensor, EXCEPT when a tensor
	/// is created by the Tensorflow code called by the
	/// RosettaTensorflowSessionContainer.  In this case, the
	/// RosettaTensorflowSessionContainer is responsible for calling
	/// update_tensor_data_pointer().
	T * tensor_data_ = nullptr;

};

} //tensorflow_manager
} //basic

#endif //USE_TENSORFLOW
#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowTensorContainer_hh
