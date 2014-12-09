// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//////////////////////////////////////////////////////////////////////
/// @begin MathVector
///
/// @brief
/// Vector0's that can perform mathmatical functions
///
/// @detailed
/// This is an implementation of an algorithm that was taken from BCL (Jens Meiler)
/// The MathVector is constructed just like utility::vector0, however this class does not
/// inherit from utility::vector0. It is implemented this way to avoid confusion. Most functions
/// from the std::vector / utility::vector0 ARE NOT included. This is a vector that performs mathematical
/// functions, not a "storage" vector. Actual mathematical functions found in numeric/MathVector_operations.
/// To access specific values you must use the operator (). For example: vector(5), will give you the value at
/// index 5. This is done to distinguish from utility::vector!
///
///
/// @references
/// Nils Woetzl
/// Jens Meiler
///
/// @authors Steven Combs, Nils Woetzl, Jens Meiler
///
/// @last_modified August 19 2010
/////////////////////////////////////////////////////////////////////////



#ifndef INCLUDED_numeric_MathVector_hh
#define INCLUDED_numeric_MathVector_hh


// Package headers
#include <numeric/types.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <math.h>
#include <numeric>
#include <algorithm> //needed for std::transform, std::find_if
#include <functional>

namespace numeric{

template<typename T>
class MathVector
{

public:

	//////////
	// data //
	//////////


	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief default constructor
	MathVector< T>() :
	size_( 0),
	data_( NULL)
	{
	}

	/// @brief construct from size and possible filler
	/// @param SIZE number fo elements in Vector
	/// @param FILL_VALUE assign every element to that value
	explicit MathVector< T>( const Size SIZE, const T &FILL_VALUE= T( 0)) :
		size_( SIZE),
		data_( new T[ SIZE])
	{

		// set all values to FILL_VALUE
		std::fill( begin(), end(), FILL_VALUE);
	}

	/// @brief construct from length and pointer to data
	MathVector< T>( const Size SIZE, const T *DATA) :
		size_( SIZE),
		data_( new T[ SIZE])
	{
		std::copy( DATA, DATA + SIZE, data_);
	}

	/// @brief copy constructor
	/// @param VECTOR copy the given Vector
	MathVector< T>( const MathVector< T> &VECTOR) :
		size_( VECTOR.size_),
		data_( new T[ size_])
	{
		std::copy( VECTOR.data_, VECTOR.data_ + size_, data_);
	}


	/// @brief Clone function
	/// @return pointer to new Vector< T>
	MathVector< T> *clone() const
	{
		return new MathVector< T>( *this);
	}

	/// @ brief destructor
	~MathVector< T>()
	{
		delete[] data_;
	}

	/////////////////
	// data access //
	/////////////////



	/// @brief size of vector
	/// @return size of Vector
	Size size() const
	{
		return size_;
	}

	/// @brief pointer to First Element
	/// @return const pointer to first element in range containing all elements of Vector
	const T *begin() const
	{
		return data_;
	}

	/// @brief pointer to First Element
	/// @return pointer to first element in range containing all elements of Vector
	T *begin()
	{
		return data_;
	}

	/// @brief pointer to end of range
	/// @return const pointer to address one after last element in Vector
	const T *end() const
	{
		return data_ + size_;
	}

	/// @brief pointer to end of range
	/// @return pointer to address one after last element in Vector
	T *end()
	{
		return data_ + size_;
	}

	/*     /// return number of elements
Size size() const
{
return size_;
}*/

	////////////////
	// operations //
	////////////////

	/// @brief norm = length of vector
	/// @return length of vector
	T norm() const
	{
		return T( sqrt( square_norm()));
	}

	/// @brief square norm = square length of vector
	/// @return square length of vector
	T square_norm() const
	{
		return std::inner_product( begin(), end(), begin(), T( 0));
	}

	/// @brief sum up all elements
	/// @return sum of all elements
	// starting value is always set to 0
	T sum() const
	{
		return std::accumulate( data_, data_ + size_, 0);
	}




	////////////////
	// operations //
	////////////////

	/// construct vector from one element
	MathVector< T> MakeVector( const T &X)
	{
		return MathVector< T>( 1, X);
	}

	/// construct vector from two elements
	inline MathVector< T> MakeVector( const T &X, const T &Y)
	{
		MathVector< T> newvector( 2);
		newvector( 0) = X;
		newvector( 1) = Y;

		return newvector;
	}

	/// construct vector from three elements
	inline	MathVector< T> MakeVector( const T &X, const T &Y, const T &Z)
	{
		MathVector< T> newvector( 3);
		newvector( 0) = X;
		newvector( 1) = Y;
		newvector( 2) = Z;

		return newvector;
	}
	///////////////
	// operators //
	///////////////




	/// return reference to changeable element ( POS)
	T &operator()( const Size POS)
	{
		assert_valid_position( POS);
		return data_[ POS];
	}

	/// return copy of element ( POS)
	const T &operator()( const Size POS) const
	{
		assert_valid_position( POS);
		return data_[ POS];
	}

	/// @brief equal operator
	/// @param VECTOR source vector
	/// @return reference to this assigned Vector
	MathVector< T> &operator =( const MathVector< T> &VECTOR)
	{
		// check that data is different
		if( data_ != VECTOR.data_)
		{
			// compare sizes
			if( size_ != VECTOR.size_)
			{
				// delete data
				delete[] data_;

				// assign size and allocate
				size_ = VECTOR.size_;
				data_ = new T[ size_];

			}

			// copy all elements
			std::copy( VECTOR.data_, VECTOR.data_ + size_, data_);
		}

		// return reference to this Vector
		return *this;
	}

	/// @brief equal operator
	/// @param VALUE all elements are set to that value
	/// @return reference to this assigned Vector
	MathVector< T> &operator =( const T &VALUE)
	{
		// set all element to given VALUE
		std::fill( data_, data_ + size_, VALUE);

		// return reference to this Vector
		return *this;
	}





	/// operator += VALUE
	MathVector<T> & operator += ( const T &VALUE)
	{
		// transform all elements
		std::transform
		(
				data_, data_ + size_, // input
				data_,                  // output
				std::binder2nd<  std::plus< T> >( std::plus< T>(), VALUE)
		);

		//end
		return *this;
	}

	/// operator -= VALUE
	MathVector< T> &operator -= ( const T &VALUE)
	{
		// transform all elements
		std::transform
		(
				data_, data_ + size_, // input
				data_,                  // output
				std::binder2nd< std::minus< T> >( std::minus< T>(), VALUE)
		);

		//end
		return *this;
	}

	/// operator *= VALUE
	MathVector< T> &operator *= ( const T &VALUE)
	{
		// transform all elements
		std::transform
		(
				data_, data_ + size_, // input
				data_,                  // output
				std::binder2nd<  std::multiplies< T> >( std::multiplies< T>(), VALUE)
		);

		//end
		return *this;
	}

	/// operator /= VALUE
	MathVector< T> &operator /= ( const T &VALUE)
	{
		// transform all elements
		std::transform
		(
				data_, data_ + size_, // input
				data_,                  // output
				std::binder2nd<  std::divides< T> >( std::divides< T>(), VALUE)
		);

		//end
		return *this;
	}


	/// operator -= BASE
	MathVector< T> &operator -= ( const MathVector< T> &BASE)
	{


		// transform
		std::transform( data_, data_ + size_, BASE.data_, data_, std::minus< T>());

		//end
		return *this;
	}

	/// operator += BASE
	MathVector< T> &operator += ( const MathVector< T> &BASE)
	{

		// transform
		std::transform( data_, data_ + size_, BASE.data_, data_, std::plus< T>());

		//end
		return *this;
	}


protected:

	/// check whether position is valid
	bool assert_valid_position( const Size POS) const
	{
		if(POS > size_){
			utility_exit_with_message("cannot access element outside of range!");
			return false;
		}
		else return true;

	}

private:

	//////////
	// data //
	//////////

	/// length of vector
	Size size_;

	/// range of dynamically allocated memory of size size_
	T *data_;



};



}




#endif /* MATHVECTOR_HH_ */
