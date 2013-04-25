#ifndef INCLUDED_ObjexxFCL_StaticIndexRange_hh
#define INCLUDED_ObjexxFCL_StaticIndexRange_hh


// StaticIndexRange: Static Index Range
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
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Dimension.hh>


namespace ObjexxFCL {


/// @brief StaticIndexRange: Static Index Range
///
/// @remarks
///  @li Zero-size range is indicated by ( l - 1 == u ) and ( size == 0 )
///  @li Upper-unbounded range is indicated by ( l - 2 == u ) and ( size == npos )
///  @li Legal ranges have ( l - 2 <= u ) with l and u in their allowed ranges
class StaticIndexRange :
	public IndexRange
{


private: // Types


	typedef  IndexRange  Super;


public: // Types


	using Super::l;
	using Super::u;

	typedef  DimensionExpression  Expression;


public: // Creation


	/// @brief Default Constructor
	inline
	StaticIndexRange()
	{}


	/// @brief Copy Constructor
	inline
	StaticIndexRange( StaticIndexRange const & I ) :
		Super( I )
	{}


	/// @brief IndexRange Constructor
	inline
	StaticIndexRange( IndexRange const & I ) :
		Super( I )
	{
		assert( I.initialized() );
	}


	/// @brief Upper Index Constructor
	inline
	StaticIndexRange( int const u_a ) :
		Super( u_a )
	{
		assert( legal_static() );
	}


	/// @brief Unbounded Upper Index Constructor
	inline
	StaticIndexRange( Star const & star ) :
		Super( star )
	{}


	/// @brief Upper Dimension Constructor
	inline
	StaticIndexRange( Dimension const & u_dim_a ) :
		Super( u_dim_a.value() )
	{
		assert( legal_static() );
	}


	/// @brief Upper Expression Constructor
	inline
	StaticIndexRange( Expression const & u_exp_a ) :
		Super( u_exp_a.ivalue() )
	{
		assert( legal_static() );
	}


	/// @brief Index Range Constructor
	inline
	StaticIndexRange( int const l_a, int const u_a ) :
		Super( l_a, u_a )
	{
		assert( legal_static() );
	}


	/// @brief Dimension Range Constructor
	inline
	StaticIndexRange( Dimension const & l_dim_a, Dimension const & u_dim_a ) :
		Super( l_dim_a.value(), u_dim_a.value() )
	{
		assert( legal_static() );
	}


	/// @brief Expression Range Constructor
	inline
	StaticIndexRange( Expression const & l_exp_a, Expression const & u_exp_a ) :
		Super( l_exp_a.ivalue(), u_exp_a.ivalue() )
	{
		assert( legal_static() );
	}


	/// @brief Index and Dimension Constructor
	inline
	StaticIndexRange( int const l_a, Dimension const & u_dim_a ) :
		Super( l_a, u_dim_a.value() )
	{
		assert( legal_static() );
	}


	/// @brief Dimension and Index Constructor
	inline
	StaticIndexRange( Dimension const & l_dim_a, int const u_a ) :
		Super( l_dim_a.value(), u_a )
	{
		assert( legal_static() );
	}


	/// @brief Index and Expression Constructor
	inline
	StaticIndexRange( int const l_a, Expression const & u_exp_a ) :
		Super( l_a, u_exp_a.ivalue() )
	{
		assert( legal_static() );
	}


	/// @brief Expression and Index Constructor
	inline
	StaticIndexRange( Expression const & l_exp_a, int const u_a ) :
		Super( l_exp_a.ivalue(), u_a )
	{
		assert( legal_static() );
	}


	/// @brief Dimension and Expression Constructor
	inline
	StaticIndexRange( Dimension const & l_dim_a, Expression const & u_exp_a ) :
		Super( l_dim_a.value(), u_exp_a.ivalue() )
	{
		assert( legal_static() );
	}


	/// @brief Expression and Dimension Constructor
	inline
	StaticIndexRange( Expression const & l_exp_a, Dimension const & u_dim_a ) :
		Super( l_exp_a.ivalue(), u_dim_a.value() )
	{
		assert( legal_static() );
	}


	/// @brief Index and Unbounded Upper Index Constructor
	inline
	StaticIndexRange( int const l_a, Star const & star ) :
		Super( l_a, star )
	{
		assert( legal_static() );
	}


	/// @brief Dimension and Unbounded Upper Index Constructor
	inline
	StaticIndexRange( Dimension const & l_dim_a, Star const & star ) :
		Super( l_dim_a.value(), star )
	{
		assert( legal_static() );
	}


	/// @brief Expression and Unbounded Upper Index Constructor
	inline
	StaticIndexRange( Expression const & l_exp_a, Star const & star ) :
		Super( l_exp_a.ivalue(), star )
	{
		assert( legal_static() );
	}


	/// @brief Destructor
	inline
	virtual
	~StaticIndexRange()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	StaticIndexRange &
	operator =( StaticIndexRange const & I )
	{
		if ( this != &I ) {
			Super::operator =( I );
			assert( legal_static() );
		}
		return *this;
	}


	/// @brief IndexRange Assignment
	inline
	StaticIndexRange &
	operator =( IndexRange const & I )
	{
		if ( this != &I ) {
			assert( I.initialized() );
			Super::operator =( I );
			assert( legal_static() );
		}
		return *this;
	}


	/// @brief Upper Index Assignment
	inline
	StaticIndexRange &
	operator =( int const u_a )
	{
		Super::operator =( u_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Unbounded Upper Index Assignment
	inline
	StaticIndexRange &
	operator =( Star const & star )
	{
		Super::operator =( star );
		return *this;
	}


	/// @brief Upper Dimension Assignment
	inline
	StaticIndexRange &
	operator =( Dimension const & u_dim_a )
	{
		Super::operator =( u_dim_a.value() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Upper Expression Assignment
	inline
	StaticIndexRange &
	operator =( Expression const & u_exp_a )
	{
		Super::operator =( u_exp_a.ivalue() );
		assert( legal_static() );
		return *this;
	}


	/// @brief StaticIndexRange Assignment
	inline
	StaticIndexRange &
	assign( StaticIndexRange const & I )
	{
		if ( this != &I ) {
			Super::operator =( I );
			assert( legal_static() );
		}
		return *this;
	}


	/// @brief IndexRange Assignment
	inline
	StaticIndexRange &
	assign( IndexRange const & I )
	{
		if ( this != &I ) {
			assert( I.initialized() );
			Super::operator =( I );
			assert( legal_static() );
		}
		return *this;
	}


	/// @brief Upper Index Assignment
	inline
	StaticIndexRange &
	assign( int const u_a )
	{
		Super::operator =( u_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Unbounded Upper Index Assignment
	inline
	StaticIndexRange &
	assign( Star const & star )
	{
		Super::operator =( star );
		return *this;
	}


	/// @brief Upper Dimension Assignment
	inline
	StaticIndexRange &
	assign( Dimension const & u_dim_a )
	{
		Super::operator =( u_dim_a.value() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Upper Expression Assignment
	inline
	StaticIndexRange &
	assign( Expression const & u_exp_a )
	{
		Super::operator =( u_exp_a.ivalue() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Index Range Assignment
	inline
	StaticIndexRange &
	assign( int const l_a, int const u_a )
	{
		Super::assign( l_a, u_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Dimension Range Assignment
	inline
	StaticIndexRange &
	assign( Dimension const & l_dim_a, Dimension const & u_dim_a )
	{
		Super::assign( l_dim_a.value(), u_dim_a.value() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Expression Range Assignment
	inline
	StaticIndexRange &
	assign( Expression const & l_exp_a, Expression const & u_exp_a )
	{
		Super::assign( l_exp_a.ivalue(), u_exp_a.ivalue() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Index and Dimension Assignment
	inline
	StaticIndexRange &
	assign( int const l_a, Dimension const & u_dim_a )
	{
		Super::assign( l_a, u_dim_a.value() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Dimension and Index Assignment
	inline
	StaticIndexRange &
	assign( Dimension const & l_dim_a, int const u_a )
	{
		Super::assign( l_dim_a.value(), u_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Index and Expression Assignment
	inline
	StaticIndexRange &
	assign( int const l_a, Expression const & u_exp_a )
	{
		Super::assign( l_a, u_exp_a.ivalue() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Expression and Index Assignment
	inline
	StaticIndexRange &
	assign( Expression const & l_exp_a, int const u_a )
	{
		Super::assign( l_exp_a.ivalue(), u_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Dimension and Expression Assignment
	inline
	StaticIndexRange &
	assign( Dimension const & l_dim_a, Expression const & u_exp_a )
	{
		Super::assign( l_dim_a.value(), u_exp_a.ivalue() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Expression and Dimension Assignment
	inline
	StaticIndexRange &
	assign( Expression const & l_exp_a, Dimension const & u_dim_a )
	{
		Super::assign( l_exp_a.ivalue(), u_dim_a.value() );
		assert( legal_static() );
		return *this;
	}


	/// @brief Index and Unbounded Upper Index Assignment
	inline
	StaticIndexRange &
	assign( int const l_a, Star const & star )
	{
		Super::assign( l_a, star );
		assert( legal_static() );
		return *this;
	}


	/// @brief Dimension and Unbounded Upper Index Assignment
	inline
	StaticIndexRange &
	assign( Dimension const & l_dim_a, Star const & star )
	{
		Super::assign( l_dim_a.value(), star );
		assert( legal_static() );
		return *this;
	}


	/// @brief Expression and Unbounded Upper Index Assignment
	inline
	StaticIndexRange &
	assign( Expression const & l_exp_a, Star const & star )
	{
		Super::assign( l_exp_a.ivalue(), star );
		assert( legal_static() );
		return *this;
	}


	/// @brief Assign Static Value of Another IndexRange: Faster Than operator =( I )
	inline
	void
	assign_value_of( IndexRange const & I )
	{ // Skips self-assignment check for speed
		l_ = I.l_;
		u_ = I.u_;
		size_ = I.size_;
		assert( legal_static() );
	}


public: // Inspector


	/// @brief Lower Index
	inline
	int
	l() const
	{
		return l_;
	}


	/// @brief Upper Index
	inline
	int
	u() const
	{
		return u_;
	}


	/// @brief Offset of an Index
	inline
	int
	offset( int const i ) const
	{
		return ( i - l_ ); // Doesn't check/require that IndexRange includes i
	}


public: // Modifier


	/// @brief Lower Index Set
	inline
	StaticIndexRange &
	l( int const l_a )
	{
		Super::l( l_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Upper Index Set
	inline
	StaticIndexRange &
	u( int const u_a )
	{
		Super::u( u_a );
		assert( legal_static() );
		return *this;
	}


	/// @brief Unbounded Upper Index Set
	inline
	StaticIndexRange &
	u( Star const & star )
	{
		Super::u( star );
		assert( legal_static() );
		return *this;
	}


	/// @brief Expand to Contain an Index
	inline
	StaticIndexRange &
	contain( int const i )
	{
		Super::contain( i );
		assert( legal_static() );
		return *this;
	}


	/// @brief Expand to Contain Another IndexRange
	inline
	StaticIndexRange &
	contain( IndexRange const & I )
	{
		Super::contain( I );
		assert( legal_static() );
		return *this;
	}


	/// @brief Intersect With Another IndexRange
	inline
	StaticIndexRange &
	intersect( IndexRange const & I )
	{
		Super::intersect( I );
		assert( legal_static() );
		return *this;
	}


	/// @brief Clear
	inline
	StaticIndexRange &
	clear()
	{
		Super::clear();
		return *this;
	}


	/// @brief Swap
	inline
	StaticIndexRange &
	swap( StaticIndexRange & I )
	{
		if ( this != &I ) {
			Super::swap( I );
			assert( legal_static() );
		}
		return *this;
	}


public: // Friend


	/// @brief Swap
	friend
	inline
	void
	swap( StaticIndexRange & a, StaticIndexRange & b )
	{
		a.swap( b );
	}


}; // StaticIndexRange


/// @brief Swap
void
swap( StaticIndexRange & a, StaticIndexRange & b );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.


namespace std {


/// @brief std::swap( StaticIndexRange, StaticIndexRange )
inline
void
swap( ObjexxFCL::StaticIndexRange & a, ObjexxFCL::StaticIndexRange & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_StaticIndexRange_HH
