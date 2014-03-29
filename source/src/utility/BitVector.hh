// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/BitVector.hh
/// @brief  Simple bit vector
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li   Wraps std::vector<bool> with a more convenient interface for bit fields
///  @li   Bits not in vector are considered false
///  @li   Bit type must be convertible to vector::size_type: making this a template arg allows the
///        use of classes with private conversion to an integral type that make BitVector a friend
///  @li   Stores a vector<bool> of the bits
///  @li   Faster but less space efficient than BitSet for sparse sets (mostly false bits)
///  @li   Bits are the index of the bits in the vector
///  @li   There must be a way to generate a BitVector from 2 Bits to construct a BitVector with more than
///        5 Bits, such as:
///             inline
///             utility::BitVector< Bit >
///             operator |( Bit const & i, Bit const & j )
///             {
///             	return utility::BitVector< Bit >( i, j );
///             }
///  @li   Construction with more than 5 Bits can be done most efficiently as
///             BitVector< Bit >( i | j |= k |= l )
///        assuming an operator| is defined as above: the use of |= instead of | after the first |
///        avoids generating additional BitVector temporaries


#ifndef INCLUDED_utility_BitVector_hh
#define INCLUDED_utility_BitVector_hh


// Unit headers
#include <utility/BitVector.fwd.hh>

// Package headers
#include <utility/utility.functions.hh>

// C++ headers
#include <vector>


namespace utility {


/// @brief Simple bit vector
template< typename B >
class BitVector
{


public: // Types


	typedef  B  Bit;
	typedef  std::vector< bool >  Bits;

	// STL/boost style
	typedef  bool  value_type;
	typedef  bool &  reference;
	typedef  bool const &  const_reference;
	typedef  bool *  pointer;
	typedef  bool const *  const_pointer;
	typedef  typename Bits::iterator  iterator;
	typedef  typename Bits::const_iterator  const_iterator;
	typedef  typename Bits::size_type  std::size_type;

	// Project style
	typedef  bool  Value;
	typedef  bool &  Reference;
	typedef  bool const &  ConstReference;
	typedef  bool *  Pointer;
	typedef  bool const *  ConstPointer;
	typedef  typename Bits::iterator  Iterator;
	typedef  typename Bits::const_iterator  ConstIterator;
	typedef  typename Bits::size_type  Size;


public: // Creation


	/// @brief Default constructor
	inline
	BitVector()
	{}


	/// @brief Bit constructor (implicit)
	inline
	BitVector( Bit const & i ) :
		bits_( i + 1 )
	{
		bits_[ i ] = true;
	}


	/// @brief 2 Bit constructor
	inline
	BitVector( Bit const & i, Bit const & j ) :
		bits_( utility::max( i, j ) + 1 )
	{
		bits_[ i ] = true;
		bits_[ j ] = true;
	}


	/// @brief 3 Bit constructor
	inline
	BitVector( Bit const & i, Bit const & j, Bit const & k ) :
		bits_( utility::max( i, j, k ) + 1 )
	{
		bits_[ i ] = true;
		bits_[ j ] = true;
		bits_[ k ] = true;
	}


	/// @brief 4 Bit constructor
	inline
	BitVector( Bit const & i, Bit const & j, Bit const & k, Bit const & l ) :
		bits_( utility::max( i, j, k, l ) + 1 )
	{
		bits_[ i ] = true;
		bits_[ j ] = true;
		bits_[ k ] = true;
		bits_[ l ] = true;
	}


	/// @brief 5 Bit constructor
	inline
	BitVector( Bit const & i, Bit const & j, Bit const & k, Bit const & l, Bit const & m ) :
		bits_( utility::max( i, j, k, l, m ) + 1 )
	{
		bits_[ i ] = true;
		bits_[ j ] = true;
		bits_[ k ] = true;
		bits_[ l ] = true;
		bits_[ m ] = true;
	}


	/// @brief Destructor
	inline
	~BitVector()
	{}


public: // Assignment


	/// @brief += BitVector: Union
	inline
	BitVector &
	operator +=( BitVector const & s )
	{
		expand_tight( s.bits_.size() );
		for ( Size i = 0, e = s.bits_.size(); i < e; ++i ) {
			bits_[ i ] = ( bits_[ i ] || s.bits_[ i ] );
		}
		return *this;
	}


	/// @brief |= BitVector: Union
	inline
	BitVector &
	operator |=( BitVector const & s )
	{
		expand_tight( s.bits_.size() );
		for ( Size i = 0, e = s.bits_.size(); i < e; ++i ) {
			bits_[ i ] = ( bits_[ i ] || s.bits_[ i ] );
		}
		return *this;
	}


	/// @brief -= BitVector: Difference
	inline
	BitVector &
	operator -=( BitVector const & s )
	{
		expand_tight( s.bits_.size() );
		for ( Size i = 0, e = s.bits_.size(); i < e; ++i ) {
			bits_[ i ] = ( bits_[ i ] && ! s.bits_[ i ] );
		}
		return *this;
	}


	/// @brief += Bit
	inline
	BitVector &
	operator +=( Bit const & i )
	{
		expand( i + 1 );
		bits_[ i ] = true;
		return *this;
	}


	/// @brief |= Bit
	inline
	BitVector &
	operator |=( Bit const & i )
	{
		expand( i + 1 );
		bits_[ i ] = true;
		return *this;
	}


	/// @brief -= Bit
	inline
	BitVector &
	operator -=( Bit const & i )
	{
		expand( i + 1 );
		bits_[ i ] = false;
		return *this;
	}


public: // Methods


	/// @brief BitVector + BitVector: Union
	friend
	inline
	BitVector
	operator +( BitVector const & a, BitVector const & b )
	{
		BitVector s( a );
		s += b;
		return s;
	}


	/// @brief BitVector | BitVector: Union
	friend
	inline
	BitVector
	operator |( BitVector const & a, BitVector const & b )
	{
		BitVector s( a );
		s |= b;
		return s;
	}


	/// @brief BitVector - BitVector: Difference
	friend
	inline
	BitVector
	operator -( BitVector const & a, BitVector const & b )
	{
		BitVector s( a );
		s -= b;
		return s;
	}


	/// @brief Shrink the bit vector to remove unused capacity
	inline
	void
	shrink()
	{
		if ( bits_.size() < bits_.capacity() ) Bits( bits_ ).swap( bits_ );
	}


	/// @brief Expand the bit vector if necessary to the specified size
	inline
	void
	expand( Size const & n )
	{
		if ( n > bits_.size() ) bits_.resize( n );
	}


	/// @brief Expand the bit vector if necessary to the specified size and remove excess capacity
	inline
	void
	expand_tight( Size const & n )
	{
		if ( n > bits_.size() ) {
			bits_.resize( n );
			shrink();
		}
	}


	/// @brief swap( BitVector )
	inline
	void
	swap( BitVector & s )
	{
		bits_.swap( s.bits_ );
	}


	/// @brief swap( BitVector, BitVector )
	friend
	inline
	void
	swap( BitVector & a, BitVector & b )
	{
		a.bits_.swap( b.bits_ );
	}


public: // Properties


	/// @brief Size
	inline
	Size
	size() const
	{
		return bits_.size();
	}


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return bits_.empty();
	}


public: // Indexers


	/// @brief BitVector[ i ] const
	inline
	bool
	operator []( Bit const & i ) const
	{
		return ( ( Size( i ) < bits_.size() ) && ( bits_[ i ] ) );
	}


public: // Iterators


	/// @brief Begin iterator
	inline
	ConstIterator
	begin() const
	{
		return bits_.begin();
	}


	/// @brief Begin iterator
	inline
	Iterator
	begin()
	{
		return bits_.begin();
	}


	/// @brief End iterator
	inline
	ConstIterator
	end() const
	{
		return bits_.end();
	}


	/// @brief End iterator
	inline
	Iterator
	end()
	{
		return bits_.end();
	}


public: // Comparison


	/// @brief BitVector == BitVector
	friend
	inline
	bool
	operator ==( BitVector const & a, BitVector const & b )
	{
		return ( a.bits_ == b.bits_ );
	}


	/// @brief BitVector != BitVector
	friend
	inline
	bool
	operator !=( BitVector const & a, BitVector const & b )
	{
		return ( a.bits_ != b.bits_ );
	}


private: // Fields


	/// @brief Bit vector
	Bits bits_;


}; // BitVector


// Friend function namespace declarations


/// @brief BitVector + BitVector: Union
template< typename B >
BitVector< B >
operator +( BitVector< B > const & a, BitVector< B > const & b );


/// @brief BitVector | BitVector: Union
template< typename B >
BitVector< B >
operator |( BitVector< B > const & a, BitVector< B > const & b );


/// @brief BitVector - BitVector: Difference
template< typename B >
BitVector< B >
operator -( BitVector< B > const & a, BitVector< B > const & b );


/// @brief swap( BitVector, BitVector )
template< typename B >
void
swap( BitVector< B > & a, BitVector< B > & b );


/// @brief BitVector == BitVector
template< typename B >
bool
operator ==( BitVector< B > const & a, BitVector< B > const & b );


/// @brief BitVector != BitVector
template< typename B >
bool
operator !=( BitVector< B > const & a, BitVector< B > const & b );


} // namespace utility


#endif // INCLUDED_utility_BitVector_HH
