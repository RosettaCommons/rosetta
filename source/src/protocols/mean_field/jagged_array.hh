// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    jagged_array.hh

/// @brief   Declarations and simple accessor/mutator definitions for jagged_array. Generally copied from utility::vector1.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_jagged_array_HH
#define INCLUDED_protocols_mean_field_jagged_array_HH

// Unit header
#include <protocols/mean_field/jagged_array.fwd.hh>

// Package headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C++ headers
#include <iostream>
#include <iomanip>

namespace protocols {
namespace mean_field {

/// @brief utility::vector1<utility::vector1>> vector of vectors to create the effect of a jagged array
template<typename T, typename A >

class jagged_array : public utility::vector1 < utility::vector1 < T, A > > {

private: // Types


	typedef utility::vector1 < utility::vector1 < T, A > >  super;
	typedef T ( *one_param_func ) (T const &);
	typedef T ( *two_param_func ) (T const &, T const &);
	typedef T ( *three_param_func) (T const &, T const &, T const &);

	template<typename T_2>
	struct func_T2
	{
		typedef T ( *three_param_func) (T const &, T const &, T_2 const &);
		typedef T ( *two_param_func) (T const &, T_2 const &);
	};

public: // Types


	// STL/boost style
	typedef  typename super::value_type  value_type;
	typedef  typename super::reference  reference;
	typedef  typename super::const_reference  const_reference;
	typedef  typename super::pointer  pointer;
	typedef  typename super::const_pointer  const_pointer;
	typedef  typename super::iterator  iterator;
	typedef  typename super::const_iterator  const_iterator;
	typedef  typename super::reverse_iterator  reverse_iterator;
	typedef  typename super::const_reverse_iterator  const_reverse_iterator;
	typedef  typename super::size_type  size_type;
	typedef  typename super::difference_type  difference_type;
	typedef  typename super::allocator_type  allocator_type;
	typedef  typename super::index_type  index_type;
	typedef  typename super::ssize_type  ssize_type;

	// Project style
	typedef  typename super::Value  Value;
	typedef  typename super::Reference  Reference;
	typedef  typename super::ConstReference  ConstReference;
	typedef  typename super::Pointer  Pointer;
	typedef  typename super::ConstPointer  ConstPointer;
	typedef  typename super::Iterator  Iterator;
	typedef  typename super::ConstIterator  ConstIterator;
	typedef  typename super::ReverseIterator  ReverseIterator;
	typedef  typename super::ConstReverseIterator  ConstReverseIterator;
	typedef  typename super::Size  Size;
	typedef  typename super::Difference  Difference;
	typedef  typename super::Allocator  Allocator;
	typedef  typename super::Index  Index;
	typedef  typename super::SSize  SSize;


public: // Methods: imports


	using super::assign;
	using super::at;
	using super::back;
	using super::begin;
	using super::capacity;
	using super::clear;
	using super::empty;
	using super::end;
	using super::erase;
	using super::front;
	using super::get_allocator;
	using super::insert;
	using super::max_size;
	using super::operator [];
	using super::operator =;
	using super::pop_back;
	using super::push_back;
	using super::rbegin;
	using super::rend;
	using super::reserve;
	using super::resize;
	using super::size;
	using super::swap;


public: // Creation


	/// @brief Default constructor
	inline
	explicit
	jagged_array( allocator_type const & alloc = allocator_type() ) :
		super( alloc )
	{}


	/// @brief Copy constructor
	inline
	jagged_array( jagged_array < T > const & v ) :
		super( v.begin(), v.end() )
	{}


	/// @brief Assignable copy constructor
	//template< ssize_type L_, typename T_, typename A_ >
	inline
	jagged_array( utility::vector1 < utility::vector1 < T, A > > const & v ) :
		super( v.begin(), v.end() )
	{}


	/// @brief Size constructor constructs an empty jagged_array of size num
	/// @param [in] num - determines size of jagged_array
	inline
	explicit
	jagged_array( size_type const num ) :
		super( num )
	{}

	/// @brief Size and vals constructor constructs a jagged_array of size num filled with copies of vals
	/// @param [in] num - determines size of jagged_array
	/// @param [in] vals - used as initial value of vector1's in jagged_array
	inline
	explicit
	jagged_array( size_type const num, utility::vector1 < T > const & vals  ) :
		super( num, vals )
	{}

	/// @brief Uniform value constructor
	/// @details constructs a vector1 of size dim that holds vector1s that each contain dim2 copies of value
	/// @param [in] dim - first dimension of jagged_array (determines size of jagged_array itself)
	/// @param [in] dim2 - second dimension of jagged_array (determines size of each vector1 within jagged_array)
	/// @remarks currently not operational due to PyRosetta build issues with this constructor
	//inline
	//jagged_array(
	// size_type const dim,
	// size_type const dim2,
	// T const & value,
	// allocator_type const & alloc = allocator_type()
	//) :
	// super( dim, super( dim2, value, alloc ), alloc )
	//{}

	/// @brief Constructor that reserves size based on vector of size types
	/// @param [in] dims - dimensions used to determine size of empty vectors within jagged_array
	inline
	jagged_array(
		utility::vector1 < size_type > const & dims
	) : super()
	{

		reserve( dims.size() );

		for ( Size vec_ind = 1 ; vec_ind <= dims.size() ; ++vec_ind ) {
			push_back( utility::vector1 < T > ( dims[ vec_ind ] ) );
		}
	}

	/// @brief Constructor that reserves size based on vector of size types, fills vectors with vals
	/// @details Constructs a vector of size dims.size().  Each vector within this vector has the size of the element
	/// @details with its index within the vector dims and is filled with the element with its index within the vector vals
	inline
	jagged_array(
		utility::vector1 < size_type > const & dims,
		utility::vector1 < T > const & vals
	) : super()
	{

		reserve( dims.size() );

		for ( Size vec_ind = 1 ; vec_ind <= dims.size() ; ++vec_ind ) {
			push_back( utility::vector1 < T > ( dims[ vec_ind ], vals[ vec_ind ] ) );
		}
	}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	jagged_array(
		InputIterator const beg,
		InputIterator const end,
		allocator_type const & alloc = allocator_type()
	) :
		super( beg, end, alloc )
	{}


	/// @brief Destructor
	inline
	virtual
	~jagged_array()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	jagged_array &
	operator =( jagged_array < T > const & v )
	{
		if ( this != &v ) {
			super::operator =( v );
		}
		return *this;
	}


	/// @brief Assignable copy assignment
	template< ssize_type L_, typename T_, typename A_ >
	inline
	jagged_array &
	operator =( utility::vector1 < utility::vector1 < T, A > > const & v )
	{
		super::assign( v.begin(), v.end() );
		return *this;
	}

	/// @brief clears jagged_array and reserves size based on vector of size types
	inline
	void assign(
		utility::vector1 < size_type > const & dims
	)
	{
		clear();

		reserve( dims.size() );

		for ( Size vec_ind = 1 ; vec_ind <= dims.size() ; ++vec_ind ) {
			push_back( utility::vector1 < T > ( dims[vec_ind] ) );
		}
	}

	/// @brief clears jagged_array and reserves size based on vector of size types, fills vectors with vals
	inline
	void assign(
		utility::vector1 < size_type > const & dims,
		utility::vector1 < T > const & vals
	)
	{
		clear();

		reserve( dims.size() );

		for ( Size vec_ind = 1 ; vec_ind <= dims.size() ; ++vec_ind ) {
			push_back( utility::vector1 < T > ( dims[ vec_ind ], vals[ vec_ind ] ) );
		}
	}

	//untested
#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		for ( Iterator iter = begin() ; iter != end() ; ++iter )
			ar & *iter;
	}
#endif

	/// @brief Find the index of an element. If not found then return 0;
	inline
	int
	index( T const & t ) const {

		Size idx = ( *this )[ 1 ].size() + 1;
		Size vec_ind = 0;

		do {
			++vec_ind;
			idx = ( ( *this )[ vec_ind ] ).index( t );
		} while ( idx == 0 && vec_ind < size() );

		return idx;
	}

	/// @brief Checks if element is present and returns bool depending whether present
	inline
	bool
	has_value( T const & t ) const {

		bool has = false;
		Size vec_ind = 0;

		while ( !has && vec_ind <= size() ) {

			has = ( *this )[ ++vec_ind ].has_value( t );

		}

		return has;
	}

	/// @brief shows a string representation of the jagged_array
	/// @details each column represents one of the utility::vector1's within the vector that comprises the jagged_array
	/// @details therefore, the "jagged" edge of the array is the bottom edge of the matrix
	/// @remarks does nothing if jagged_array is empty
	inline
	void
	show(std::ostream & output=std::cout) const {

		if ( ! empty() ) {
			for ( Size ii = 1 ; ii <= max_size_col() ; ++ii ) {
				for ( Size jj = 1 ; jj <= size() ; ++jj ) {
					if ( ii <= ( *this )[ jj ].size() ) {
						output << ( *this )[ jj ][ ii ] << "\t";
					} else {
						output << "\t\t";
					}
				}
				output << std::endl;
			}
		}
	}

	/// @brief applies param *func to all values in the jagged_array (1-param func)
	//inline
	//void
	//apply_func_to_all( one_param_func func ) {

	// for ( Size ii = 1 ; ii <= size() ; ++ii ) {
	//  for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
	//   ( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ] );
	//  }
	// }

	//}

	/// @brief applies param *func to all values in the jagged_array (2-param func)
	/// @remarks currently not operational due to cppcheck errors
	//inline
	//void
	//apply_func_to_all( two_param_func func, T operand ) {
	// for ( Size ii = 1 ; ii <= size() ; ++ii ) {
	//  for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
	//   ( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ] , operand );
	//  }
	// }
	//}

	/// @brief applies param *func to all values in the jagged_array (2-param func) using different template variable as other operand
	template < typename T_2 >
	inline
	void
	apply_func_to_all( typename func_T2<T_2>::two_param_func func, T_2 operand ) {
		for ( Size ii = 1 ; ii <= size() ; ++ii ) {
			for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
				( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ] , operand );
			}
		}
	}

	/// @brief applies param *func to all values in the jagged_array (2-param func) using column_dependent operand
	/// @remarks currently not operational due to cppcheck errors
	//inline
	//void
	//apply_func_to_each_col ( two_param_func func, utility::vector1 < T > operands ) {
	//
	// assert( size() == operands.size() );

	// for ( Size ii = 1 ; ii <= size() ; ++ii ) {
	//  for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
	//   ( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ], operands[ ii ] );
	//  }
	// }

	//}


	/// @brief applies param *func to all values in the jagged_array (2-param func) using different template variable as other, column_dependent operand
	template < typename T_2 >
	inline
	void
	apply_func_to_each_col ( typename func_T2<T_2>::two_param_func func, utility::vector1 < T_2 > operands ) {

		assert( size() == operands.size() );

		for ( Size ii = 1 ; ii <= size() ; ++ii ) {
			for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
				( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ], operands[ ii ] );
			}
		}

	}

	/// @brief applies param *func to all values in the jagged_array (2-param func) using 2D-index dependent operand (i.e. matrix arithmetic)
	/// @remarks currently not operational due to cppcheck errors
	//inline
	//void
	//apply_func_to_each_elem ( two_param_func func, jagged_array < T > operands ) {

	// assert( equal_size( operands ) );

	// for ( Size ii = 1 ; ii <= size() ; ++ii ) {
	//  for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
	//   ( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ], operands[ ii ][ jj ] );
	//  }
	// }

	//}

	/// @brief applies param *func to all values in the jagged_array (2-param func) using different template variable as other, 2D-index dependent operand
	template < typename T_2 >
	inline
	void
	apply_func_to_each_elem ( typename func_T2<T_2>::two_param_func func, jagged_array < T_2 > operands ) {

		assert( equal_size( operands ) );

		for ( Size ii = 1 ; ii <= size() ; ++ii ) {
			for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
				( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ], operands[ ii ][ jj ] );
			}
		}

	}

	/// @brief applies param *func to all values in the jagged_array (3-param func) using 2D-index dependent operand and a third, general operand
	/// @remarks currently not operational due to cppcheck errors
	//inline
	//void
	//apply_func_to_each_elem ( three_param_func func, jagged_array < T > operands, T third_oper ) {

	// assert( equal_size( operands ) );

	// for ( Size ii = 1 ; ii <= size() ; ++ii ) {
	//  for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
	//   ( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ], operands[ ii ][ jj ], third_oper );
	//  }
	// }

	//}

	/// @brief applies param *func to all values in the jagged_array (2-param func) using 2D-index dependent operand and a different template variable operand
	template < typename T_2 >
	inline
	void
	apply_func_to_each_elem ( typename func_T2<T_2>::three_param_func func, jagged_array < T > operands, T_2 third_oper ) {

		assert( equal_size( operands ) );

		for ( Size ii = 1 ; ii <= size() ; ++ii ) {
			for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
				( *this ) [ ii ][ jj ] = ( *func )( ( *this ) [ ii ][ jj ], operands[ ii ][ jj ], third_oper );
			}
		}

	}

	/// @brief returns vector1 of totals of columns of jagged_array
	/// @details class T must provide operator +=
	/// @details Error if size is 0
	/// @remarks very useful in mean-field algorithm
	inline
	utility::vector1 < T >
	get_totals_columns () const {

		assert( size() > 0 );

		utility::vector1 < T > totals( size() );

		for ( Size ii = 1 ; ii <= size() ; ++ii ) {
			T total = ( *this )[ ii ][1];

			for ( Size jj = 2 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
				total += ( *this )[ ii ][ jj ];
			}

			totals[ ii ] = total;
		}
		return totals;
	}

	/// @brief returns total of all elements in jagged array
	/// @details class T must provide operator += and -=
	/// @details Error if size is 0
	inline
	T
	get_total () const {

		assert( size() > 0 );

		T total;

		if ( size() > 0 ) {
			total = ( *this )[1][1];
			for ( Size ii = 1 ; ii <= size() ; ++ii ) {
				for ( Size jj = 1 ; jj <= ( *this )[ ii ].size() ; ++jj ) {
					total+=( *this )[ ii ][ jj ];
				}
			}
			total -= ( *this )[1][1];
		}
		return total;
	}

	/// @brief returns size of specified column
	/// @param [in] col - column for which to determine size
	/// @details Error if col > size
	inline
	Size
	size_col ( Size col ) const {

		assert( size() >= col );

		return ( *this )[ col ].size();
	}

	/// @brief returns the maximum size of any columns
	/// @details Error if jagged_array is empty
	inline
	Size
	max_size_col () const {

		assert( ! empty() );

		Size max_size = size_col( 1 );

		for ( Size i = 2; i <= size(); ++i ) {
			if ( max_size < size_col( i ) ) max_size = size_col( i );
		}

		return max_size;
	}

	/// @brief checks if a second jagged_array is of equal size for all columns
	/// @param [in] other - second jagged_array with which to compare
	/// @remarks useful for matrix arithmetic
	inline
	bool
	equal_size ( jagged_array <T> const & other ) const
	{
		if ( size() != other.size() ) return false;

		for ( Size ii = 1; ii <= size(); ++ii ) {
			if ( ( *this )[ ii ].size() != other[ ii ].size() ) return false;
		}

		return true;
	}

	/// @brief checks if a second jagged_array (of different type) is of equal size for all columns
	/// @param [in] other - second jagged_array with which to compare
	/// @remarks useful for matrix arithmetic
	template < typename T_2 >
	inline
	bool
	equal_size ( jagged_array <T_2> const & other ) const
	{
		if ( size() != other.size() ) return false;

		for ( Size ii = 1; ii <= size(); ++ii ) {
			if ( ( *this )[ ii ].size() != other[ ii ].size() ) return false;
		}

		return true;
	}

}; // jagged_array

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_jagged_array_HH
