// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/BitSet.hh
/// @brief  Simple bit set
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_BitSet_hh
#define INCLUDED_utility_BitSet_hh


// Unit headers
#include <utility/BitSet.fwd.hh>

// C++ headers
#include <set>


namespace utility {


/// @brief Simple bit set
/// @note
///  @li   Bit set with a convenient interface for bit fields
///  @li   Stores a set of the bit numbers that are on (present) so space efficient for sparse sets
///        but could have slower lookup for large sets than std::bitset, vector<bool>, or BitVector
///  @li   Bits not in set are considered false
///  @li   Bits are the index of the stored (true) bits
///  @li   There must be a way to generate a BitSet from 2 Bits to construct a BitSet with more than
///        5 Bits, such as:
///             inline
///             utility::BitSet< Bit >
///             operator |( Bit const & i, Bit const & j )
///             {
///             	return utility::BitSet< Bit >( i, j );
///             }
///  @li   Construction with more than 5 Bits can be done most efficiently as
///             BitSet< Bit >( i | j |= k |= l )
///        assuming an operator| is defined as above: the use of |= instead of | after the first |
///        avoids generating additional BitSet temporaries
template< typename B >
class BitSet
{


public: // Types


	typedef  std::set< B >  Bits;

	// STL/boost style
	typedef  B  value_type;
	typedef  B &  reference;
	typedef  B const &  const_reference;
	typedef  B *  pointer;
	typedef  B const *  const_pointer;
	typedef  typename Bits::iterator  iterator;
	typedef  typename Bits::const_iterator  const_iterator;
	typedef  typename Bits::size_type  std::size_type;

	// Project style
	typedef  B  Bit;
	typedef  B &  Reference;
	typedef  B const &  ConstReference;
	typedef  B *  Pointer;
	typedef  B const *  ConstPointer;
	typedef  typename Bits::iterator  Iterator;
	typedef  typename Bits::const_iterator  ConstIterator;
	typedef  typename Bits::size_type  Size;


public: // Creation


	/// @brief Default constructor
	inline
	BitSet()
	{}


	/// @brief Bit constructor (implicit)
	inline
	BitSet( Bit const & i )
	{
		bits_.insert( i );
	}


	/// @brief 2 Bit constructor
	inline
	BitSet( Bit const & i, Bit const & j )
	{
		bits_.insert( i );
		bits_.insert( j );
	}


	/// @brief 3 Bit constructor
	inline
	BitSet( Bit const & i, Bit const & j, Bit const & k )
	{
		bits_.insert( i );
		bits_.insert( j );
		bits_.insert( k );
	}


	/// @brief 4 Bit constructor
	inline
	BitSet( Bit const & i, Bit const & j, Bit const & k, Bit const & l )
	{
		bits_.insert( i );
		bits_.insert( j );
		bits_.insert( k );
		bits_.insert( l );
	}


	/// @brief 5 Bit constructor
	inline
	BitSet( Bit const & i, Bit const & j, Bit const & k, Bit const & l, Bit const & m )
	{
		bits_.insert( i );
		bits_.insert( j );
		bits_.insert( k );
		bits_.insert( l );
		bits_.insert( m );
	}


	/// @brief Destructor
	inline
	~BitSet()
	{}


public: // Assignment


	/// @brief += BitSet: Union
	inline
	BitSet &
	operator +=( BitSet const & s )
	{
		bits_.insert( s.begin(), s.end() );
		return *this;
	}


	/// @brief |= BitSet: Union
	inline
	BitSet &
	operator |=( BitSet const & s )
	{
		bits_.insert( s.begin(), s.end() );
		return *this;
	}


	/// @brief -= BitSet: Difference
	inline
	BitSet &
	operator -=( BitSet const & s )
	{
		for ( ConstIterator i = s.begin(), e = s.end(); i != e; ++i ) {
			bits_.erase( *i );
		}
		return *this;
	}


	/// @brief += Bit
	inline
	BitSet &
	operator +=( Bit const & i )
	{
		bits_.insert( i );
		return *this;
	}


	/// @brief |= Bit
	inline
	BitSet &
	operator |=( Bit const & i )
	{
		bits_.insert( i );
		return *this;
	}


	/// @brief -= Bit
	inline
	BitSet &
	operator -=( Bit const & i )
	{
		bits_.erase( i );
		return *this;
	}


public: // Methods


	/// @brief BitSet + BitSet: Union
	friend
	inline
	BitSet
	operator +( BitSet const & a, BitSet const & b )
	{
		BitSet s( a );
		s.bits_.insert( b.begin(), b.end() );
		return s;
	}


	/// @brief BitSet | BitSet: Union
	friend
	inline
	BitSet
	operator |( BitSet const & a, BitSet const & b )
	{
		BitSet s( a );
		s.bits_.insert( b.begin(), b.end() );
		return s;
	}


	/// @brief BitSet - BitSet: Difference
	friend
	inline
	BitSet
	operator -( BitSet const & a, BitSet const & b )
	{
		BitSet s( a );
		for ( ConstIterator i = b.begin(), e = b.end(); i != e; ++i ) {
			s.bits_.erase( *i );
		}
		return s;
	}


	/// @brief swap( BitSet )
	inline
	void
	swap( BitSet & s )
	{
		bits_.swap( s.bits_ );
	}


	/// @brief swap( BitSet, BitSet )
	friend
	inline
	void
	swap( BitSet & a, BitSet & b )
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


	/// @brief BitSet[ i ] const
	inline
	bool
	operator []( Bit const & i ) const
	{
		return ( bits_.find( i ) != bits_.end() );
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


	/// @brief BitSet == BitSet
	friend
	inline
	bool
	operator ==( BitSet const & a, BitSet const & b )
	{
		return ( a.bits_ == b.bits_ );
	}


	/// @brief BitSet != BitSet
	friend
	inline
	bool
	operator !=( BitSet const & a, BitSet const & b )
	{
		return ( a.bits_ != b.bits_ );
	}


private: // Fields


	/// @brief Bit set
	Bits bits_;


}; // BitSet


// Friend function namespace declarations


/// @brief BitSet + BitSet: Union
template< typename B >
BitSet< B >
operator +( BitSet< B > const & a, BitSet< B > const & b );


/// @brief BitSet | BitSet: Union
template< typename B >
BitSet< B >
operator |( BitSet< B > const & a, BitSet< B > const & b );


/// @brief BitSet - BitSet: Difference
template< typename B >
BitSet< B >
operator -( BitSet< B > const & a, BitSet< B > const & b );


/// @brief swap( BitSet, BitSet )
template< typename B >
void
swap( BitSet< B > & a, BitSet< B > & b );


/// @brief BitSet == BitSet
template< typename B >
bool
operator ==( BitSet< B > const & a, BitSet< B > const & b );


/// @brief BitSet != BitSet
template< typename B >
bool
operator !=( BitSet< B > const & a, BitSet< B > const & b );


} // namespace utility


#endif // INCLUDED_utility_BitSet_HH
