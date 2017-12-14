// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/util.cc
/// @brief small bundle of utilities for dealing with numbers.
/// @author James Thompson

#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <numeric/MathNTensorBase.hh>
#include <numeric/MathNTensor.hh>

#include <algorithm>

namespace numeric {

numeric::Real median( utility::vector1< numeric::Real > const & values ) {
	assert( values.size() ); // An empty list doesn't have a median
	utility::vector1< numeric::Real > vals = values;
	std::sort( vals.begin(), vals.end() );

	numeric::Size const n_vals( vals.size() );
	numeric::Real retval( 0.0 );
	if ( n_vals % 2 == 0 ) { // Even number of items
		retval += 0.5 * vals[ n_vals / 2 ];
		retval += 0.5 * vals[ n_vals / 2 + 1];
	} else { // Odd number of items
		retval = vals[ (n_vals - 1) / 2 + 1 ];
	}
	return retval;
}

numeric::Real mean( utility::vector1< numeric::Real > const & values ) {

	numeric::Size const n_vals( values.size() );
	numeric::Real total( 0.0 );
	for ( double value : values ) {
		total += value;
	}

	return static_cast< numeric::Real > ( total / n_vals );
}

/// @brief Utility function to access an entry in a MathNTensor of arbitrary dimensionality unknown at compile time, given a MathNTensorBaseOP.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Real &
access_Real_MathNTensor( MathNTensorBaseOP< Real > tensorbase, utility::vector1< Size > const &position ) {
	Size const dimension( tensorbase->dimensionality() );
	runtime_assert_string_msg( dimension == position.size(), "Error in numeric::access_Real_MathNTensor(): The position vector does not match the dimensionality of the tensor." );

	switch( dimension ) {
	case 1 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 1 > fixpos;
		fixpos[1] = position[1];
		MathNTensorOP<Real, 1> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 1> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
		break;
	case 2 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 2 > fixpos;
		fixpos[1] = position[1];
		fixpos[2] = position[2];
		MathNTensorOP<Real, 2> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 2> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
		break;
	case 3 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 3 > fixpos;
		for ( Size i=1; i<=3; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 3> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 3> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 4 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 4 > fixpos;
		for ( Size i=1; i<=4; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 4> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 4> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 5 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 5 > fixpos;
		for ( Size i=1; i<=5; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 5> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 5> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 6 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 6 > fixpos;
		for ( Size i=1; i<=6; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 6> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 6> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 7 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 7 > fixpos;
		for ( Size i=1; i<=7; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 7> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 7> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 8 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 8 > fixpos;
		for ( Size i=1; i<=8; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 8> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 8> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 9 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 9 > fixpos;
		for ( Size i=1; i<=9; ++i ) fixpos[i] = position[i];
		MathNTensorOP<Real, 9> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 9> >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
		break;
	default :
		utility_exit_with_message( "Error in numeric::access_Real_MathNTensor(): This access function only works for tensors of dimensionality 1 to 9." );
	}

	//The following should never end up being called, and would produce a memory leak if it were:
	auto * dummy = new Real(0.0);
	return (*dummy);
}


/// @brief Utility function to access an entry in a MathNTensor of arbitrary dimensionality unknown at compile time, given a MathNTensorBaseCOP.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Real const &
const_access_Real_MathNTensor( MathNTensorBaseCOP< Real > tensorbase, utility::vector1< Size > const &position ) {
	Size const dimension( tensorbase->dimensionality() );
	runtime_assert_string_msg( dimension == position.size(), "Error in numeric::const_access_Real_MathNTensor(): The position vector does not match the dimensionality of the tensor." );

	switch( dimension ) {
	case 1 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 1 > fixpos;
		fixpos[1] = position[1];
		MathNTensorCOP<Real, 1> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 1> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
		break;
	case 2 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 2 > fixpos;
		fixpos[1] = position[1];
		fixpos[2] = position[2];
		MathNTensorCOP<Real, 2> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 2> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
		break;
	case 3 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 3 > fixpos;
		for ( Size i=1; i<=3; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 3> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 3> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 4 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 4 > fixpos;
		for ( Size i=1; i<=4; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 4> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 4> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 5 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 5 > fixpos;
		for ( Size i=1; i<=5; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 5> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 5> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 6 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 6 > fixpos;
		for ( Size i=1; i<=6; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 6> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 6> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 7 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 7 > fixpos;
		for ( Size i=1; i<=7; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 7> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 7> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 8 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 8 > fixpos;
		for ( Size i=1; i<=8; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 8> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 8> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
	case 9 :
		{ //Scope for variable declaration.
		utility::fixedsizearray1< Size, 9 > fixpos;
		for ( Size i=1; i<=9; ++i ) fixpos[i] = position[i];
		MathNTensorCOP<Real, 9> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 9> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return (*tensor_ptr)( fixpos );
	}
		break;
	default :
		utility_exit_with_message( "Error in numeric::const_access_Real_MathNTensor(): This access function only works for tensors of dimensionality 1 to 9." );
	}

	//The following should never end up being called, and would produce a memory leak if it were:
	auto * dummy = new Real(0.0);
	return (*dummy);
}

/// @brief Given a MathNTensorBaseCOP, get the size along one dimension.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Size
get_Real_MathNTensor_dimension_size(
	MathNTensorBaseCOP< Real > tensorbase,
	Size const dimension_index
) {
	Size const dimension( tensorbase->dimensionality() );
	runtime_assert_string_msg( dimension_index > 0, "Error in numeric::get_Real_MathNTensor_dimension_size(): The dimension index must be greater than zero." );
	runtime_assert_string_msg( dimension_index <= dimension, "Error in numeric::get_Real_MathNTensor_dimension_size(): The dimension index must be less than the dimensionality of the tensor." );
	switch( dimension ) {
	case 1 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 1> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 1> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 2 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 2> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 2> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 3 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 3> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 3> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 4 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 4> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 4> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 5 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 5> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 5> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 6 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 6> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 6> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 7 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 7> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 7> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 8 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 8> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 8> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	case 9 :
		{ //Scope for variable declaration.
		MathNTensorCOP<Real, 9> tensor_ptr( utility::pointer::dynamic_pointer_cast< MathNTensor< Real, 9> const >(tensorbase) );
		runtime_assert( tensor_ptr );
		return tensor_ptr->n_bins( dimension_index );
	}
		break;
	default :
		utility_exit_with_message( "Error in numeric::get_Real_MathNTensor_dimension_size(): This access function only works for tensors of dimensionality 1 to 9." );
	} //switch

	return 0; //Should never occur, but needed to satisfy compiler.
}


} // numeric
