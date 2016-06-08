#ifndef INCLUDED_ObjexxFCL_IndexRange_hh
#define INCLUDED_ObjexxFCL_IndexRange_hh


// IndexRange: Index Range Abstract Base Class
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/Dimension.fwd.hh>

// C++ Headers
#include <cassert>
#include <cstddef>
#include <iosfwd>
// only available in cxx11
#ifdef CXX11
#include <utility> // for swap
#else
#include <algorithm>
#endif


namespace ObjexxFCL {

/// @brief IndexRange: Index Range Abstract Base Class
///
/// @remarks
///  @li Zero-size range is indicated by ( l - 1 == u ) and ( size == 0 )
///  @li Upper-unbounded range is indicated by ( l - 2 == u ) and ( size == npos )
///  @li Legal ranges have ( l - 2 <= u ) with l and u in their allowed ranges
class IndexRange
{


private: // Friend


	friend class StaticIndexRange;
	friend class DynamicIndexRange;


public: // Types


	// STL style
	typedef  std::size_t  size_type;

	// C++ style
	typedef  std::size_t  Size;


protected: // Creation


	/// @brief Default Constructor
	inline
	IndexRange() :
		l_( 1 ),
		u_( 0 ),
		size_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	IndexRange( IndexRange const & I ) :
		l_( I.l_ ),
		u_( I.u_ ),
		size_( I.size_ )
	{}


	/// @brief Upper Index Constructor
	inline
	IndexRange( int const u_a ) :
		l_( 1 ),
		u_( u_a ),
		size_( u_ )
	{}


	/// @brief Unbounded Upper Index Constructor
	inline
	IndexRange( Star const & ) :
		l_( 1 ),
		u_( -1 ),
		size_( npos )
	{}


	/// @brief Index Range Constructor
	inline
	IndexRange( int const l_a, int const u_a ) :
		l_( l_a ),
		u_( u_a ),
		size_( u_ - l_ + 1 )
	{}


	/// @brief Index and Unbounded Upper Index Constructor
	inline
	IndexRange( int const l_a, Star const & ) :
		l_( l_a ),
		u_( l_ - 2 ),
		size_( npos )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~IndexRange()
	{}


protected: // Assignment


	/// @brief Copy Assignment
	inline
	IndexRange &
	operator =( IndexRange const & I )
	{
		l_ = I.l_;
		u_ = I.u_;
		size_ = I.size_;
		return *this;
	}


public: // Assignment


	/// @brief Upper Index Assignment
	inline
	virtual
	IndexRange &
	operator =( int const u_a )
	{
		l_ = 1;
		u_ = u_a;
		size_ = u_;
		return *this;
	}


	/// @brief Unbounded Upper Index Assignment
	inline
	virtual
	IndexRange &
	operator =( Star const & )
	{
		l_ = 1;
		u_ = -1;
		size_ = npos;
		return *this;
	}


	/// @brief Upper Index Assignment
	inline
	virtual
	IndexRange &
	assign( int const u_a )
	{
		l_ = 1;
		u_ = u_a;
		size_ = u_;
		return *this;
	}


	/// @brief Unbounded Upper Index Assignment
	inline
	virtual
	IndexRange &
	assign( Star const & )
	{
		l_ = 1;
		u_ = -1;
		size_ = npos;
		return *this;
	}


	/// @brief Index Range Assignment
	inline
	virtual
	IndexRange &
	assign( int const l_a, int const u_a )
	{
		l_ = l_a;
		u_ = u_a;
		size_ = u_ - l_ + 1;
		return *this;
	}


	/// @brief Index and Unbounded Upper Index Assignment
	inline
	virtual
	IndexRange &
	assign( int const l_a, Star const & )
	{
		l_ = l_a;
		u_ = l_ - 2;
		size_ = npos;
		return *this;
	}


public: // Predicate


	/// @brief Initialized?
	inline
	virtual
	bool
	initialized() const
	{
		return true;
	}


	/// @brief Lower Initialized?
	inline
	virtual
	bool
	l_initialized() const
	{
		return true;
	}


	/// @brief Upper Initialized?
	inline
	virtual
	bool
	u_initialized() const
	{
		return true;
	}


	/// @brief Legal?
	inline
	virtual
	bool
	legal() const
	{
		return ( ( l_ >= l_min ) && ( u_ <= u_max ) && ( l_ - 2 <= u_ ) );
	}


	/// @brief Bounded?
	inline
	virtual
	bool
	bounded() const
	{
		return ( l_ - 1 <= u_ );
	}


	/// @brief Bounded?
	inline
	bool
	bounded_value() const
	{
		return ( l_ - 1 <= u_ );
	}


	/// @brief Unbounded?
	inline
	virtual
	bool
	unbounded() const
	{
		return ( l_ - 2 == u_ );
	}


	/// @brief Unbounded?
	inline
	bool
	unbounded_value() const
	{
		return ( l_ - 2 == u_ );
	}


	/// @brief Not Unbounded?
	inline
	virtual
	bool
	not_unbounded() const
	{
		return ( l_ - 1 <= u_ );
	}


	/// @brief Bounded with Positive Size?
	inline
	virtual
	bool
	positive() const
	{
		return ( l_ <= u_ );
	}


	/// @brief Bounded with Positive Size?
	inline
	bool
	positive_value() const
	{
		return ( l_ <= u_ );
	}


	/// @brief Contains an Index?
	inline
	virtual
	bool
	contains( int const i ) const
	{
		return ( ( l_ <= i ) && ( ( i <= u_ ) || ( size_ == npos ) ) );
	}


	/// @brief Contains Another IndexRange?
	virtual
	bool
	contains( IndexRange const & I ) const;


	/// @brief Intersects Another IndexRange?
	virtual
	bool
	intersects( IndexRange const & I ) const;


public: // Inspector


	/// @brief Lower Index
	inline
	int
	l() const
	{
		assert( l_initialized() );
		return l_;
	}


	/// @brief Lower Index (Zero if Uninitialized)
	inline
	int
	lz() const
	{
		return l_;
	}


	/// @brief Upper Index
	inline
	int
	u() const
	{
		assert( u_initialized() );
		return u_;
	}


	/// @brief Upper Index (Zero if Uninitialized)
	inline
	int
	uz() const
	{
		return u_;
	}


	/// @brief Size
	inline
	size_type
	size() const
	{
		return size_; // Unbounded => npos
	}


	/// @brief Offset of an Index
	inline
	int
	offset( int const i ) const
	{
		assert( l_initialized() );
		return ( i - l_ ); // Doesn't check/require that IndexRange includes i
	}


public: // Modifier


	/// @brief Lower Index Set
	inline
	virtual
	IndexRange &
	l( int const l_a )
	{
		if ( l_ - 2 == u_ ) { // Unbounded range
			l_ = l_a;
			u_ = l_ - 2; // Reset u_ to maintain unbounded state
		} else { // Bounded
			l_ = l_a;
			size_ = u_ - l_ + 1;
		}
		return *this;
	}


	/// @brief Upper Index Set
	inline
	virtual
	IndexRange &
	u( int const u_a )
	{
		u_ = u_a;
		size_ = u_ - l_ + 1;
		return *this;
	}


	/// @brief Unbounded Upper Index Set
	inline
	virtual
	IndexRange &
	u( Star const & )
	{
		u_ = l_ - 2;
		size_ = npos;
		return *this;
	}


	/// @brief Expand to Contain an Index
	inline
	virtual
	IndexRange &
	contain( int const i )
	{
		if ( l_ - 1 <= u_ ) { // Bounded
			if ( l_ > i ) l_ = i;
			if ( u_ < i ) u_ = i;
			size_ = u_ - l_ + 1;
		} else { // Unbounded
			if ( l_ > i ) {
				l_ = i;
				u_ = l_ - 2; // Reset u_ to maintain unbounded state
			}
		}
		return *this;
	}


	/// @brief Expand to Contain Another IndexRange
	virtual
	IndexRange &
	contain( IndexRange const & I );


	/// @brief Intersect With Another IndexRange
	virtual
	IndexRange &
	intersect( IndexRange const & I );


	/// @brief Clear
	inline
	virtual
	IndexRange &
	clear()
	{
		l_ = 1;
		u_ = 0;
		size_ = 0;
		return *this;
	}


public: // Comparison


	/// @brief IndexRange == IndexRange
	friend
	inline
	bool
	operator ==( IndexRange const & I, IndexRange const & J )
	{
		return ( ( I.initialized() ) && ( J.initialized() ) && ( I.l_ == J.l_ ) && ( I.u_ == J.u_ ) );
	}


	/// @brief IndexRange != IndexRange
	friend
	inline
	bool
	operator !=( IndexRange const & I, IndexRange const & J )
	{
		return !( I == J );
	}


	/// @brief IndexRange < IndexRange
	friend
	inline
	bool
	operator <( IndexRange const & I, IndexRange const & J )
	{
		return ( ( I.positive() ) && ( J.positive() ) && ( I.u_ < J.l_ ) );
	}


	/// @brief IndexRange <= IndexRange
	friend
	inline
	bool
	operator <=( IndexRange const & I, IndexRange const & J )
	{
		return ( ( I.positive() ) && ( J.positive() ) && ( I.u_ <= J.l_ ) );
	}


	/// @brief IndexRange > IndexRange
	friend
	inline
	bool
	operator >( IndexRange const & I, IndexRange const & J )
	{
		return ( ( I.positive() ) && ( J.positive() ) && ( I.l_ > J.u_ ) );
	}


	/// @brief IndexRange >= IndexRange
	friend
	inline
	bool
	operator >=( IndexRange const & I, IndexRange const & J )
	{
		return ( ( I.positive() ) && ( J.positive() ) && ( I.l_ >= J.u_ ) );
	}


public: // I/O


	/// @brief Stream Input
	friend
	std::istream &
	operator >>( std::istream & stream, IndexRange & I );


	/// @brief Stream Output
	friend
	std::ostream &
	operator <<( std::ostream & stream, IndexRange const & I );


protected: // Inspector


	/// @brief Legal Static Range?
	inline
	bool
	legal_static() const
	{
		return ( ( l_ >= l_min ) && ( u_ <= u_max ) && ( l_ - 2 <= u_ ) );
	}


	/// @brief Lower Dimension Clone
	inline
	virtual
	Dimension *
	l_dim_clone() const
	{
		return 0;
	}


	/// @brief Upper Dimension Clone
	inline
	virtual
	Dimension *
	u_dim_clone() const
	{
		return 0;
	}


protected: // Modifier


	/// @brief Swap
	inline
	void
	swap( IndexRange & I )
	{
		if ( this != &I ) {
			std::swap( l_, I.l_ );
			std::swap( u_, I.u_ );
			std::swap( size_, I.size_ );
		}
	}


public: // Data


	static size_type const npos = static_cast< size_type >( -1 ); // Unbounded "size"

	static int const l_min = -( static_cast< int >( ( static_cast< unsigned int >( -1 ) / 2u ) ) - 1 ); // Min lower index

	static int const u_max = static_cast< int >( ( static_cast< unsigned int >( -1 ) / 2u ) ); // Max upper index


private: // Data


	/// @brief Lower index
	int l_;

	/// @brief Upper index
	int u_;

	/// @brief Size (npos iff unbounded)
	size_type size_;


}; // IndexRange


/// @brief IndexRange == IndexRange
#ifndef __clang__
bool
operator ==( IndexRange const & I, IndexRange const & J );
#endif


/// @brief IndexRange != IndexRange
bool
operator !=( IndexRange const & I, IndexRange const & J );


/// @brief IndexRange < IndexRange
bool
operator <( IndexRange const & I, IndexRange const & J );


/// @brief IndexRange <= IndexRange
bool
operator <=( IndexRange const & I, IndexRange const & J );


/// @brief IndexRange > IndexRange
bool
operator >( IndexRange const & I, IndexRange const & J );


/// @brief IndexRange >= IndexRange
bool
operator >=( IndexRange const & I, IndexRange const & J );


/// @brief Stream Input
std::istream &
operator >>( std::istream & stream, IndexRange & I );


/// @brief Stream Output
std::ostream &
operator <<( std::ostream & stream, IndexRange const & I );


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_IndexRange_HH
