#ifndef INCLUDED_ObjexxFCL_FArray1_hh
#define INCLUDED_ObjexxFCL_FArray1_hh


// FArray1: Fortran-Compatible 1D Array Abstract Base Class
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
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray.hh>

// C++ Headers
#include <cmath>


namespace ObjexxFCL {


// Forward Declarations
template< typename > class FArray1D;
template< typename > class FArray1P;
template< typename > class FArray1A;
template< typename > class KeyFArray1D;


/// @brief FArray1: Fortran-Compatible 1D Array Abstract Base Class
template< typename T >
class FArray1 :
	public FArray< T >
{


private: // Types


	typedef  FArray< T >  Super;
	typedef  FArray1D< T >  real_FArray;
	typedef  FArray1P< T >  proxy_FArray;
	typedef  FArray1A< T >  arg_FArray;


private: // Friend


	template< typename > friend class FArray1;
	template< typename > friend class FArray1D;
	template< typename > friend class FArray1P;
	template< typename > friend class FArray1A;
	template< typename > friend class KeyFArray1D;


protected: // Types


	typedef  internal::InitializerSentinel  InitializerSentinel;
	typedef  internal::ProxySentinel  ProxySentinel;


public: // Types


	typedef  typename Super::Base  Base;
	typedef  typename Base::Section  Section;
	typedef  typename Base::IR  IR;

	// STL Style
	typedef  typename Base::value_type  value_type;
	typedef  typename Base::reference  reference;
	typedef  typename Base::const_reference  const_reference;
	typedef  typename Base::pointer  pointer;
	typedef  typename Base::const_pointer  const_pointer;
	typedef  typename Base::size_type  size_type;
	typedef  typename Base::difference_type  difference_type;

	// C++ Style
	typedef  typename Base::Value  Value;
	typedef  typename Base::Reference  Reference;
	typedef  typename Base::ConstReference  ConstReference;
	typedef  typename Base::Pointer  Pointer;
	typedef  typename Base::ConstPointer  ConstPointer;
	typedef  typename Base::Size  Size;
	typedef  typename Base::Difference  Difference;

	using Super::array_;
	using Super::array_size_;
	using Super::dimensions_initialized;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
	using Super::not_const_proxy;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	using Super::npos;
	using Super::sarray_;
	using Super::shift_;

	// Types to prevent compile failure when std::distance is in scope
	typedef  void  iterator_category;


protected: // Creation


	/// @brief Default Constructor
	inline
	FArray1()
	{}


	/// @brief Copy Constructor
	inline
	FArray1( FArray1 const & a ) :
		Super( a )
	{}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	FArray1( FArray1< U > const & a ) :
		Super( a )
	{}


	/// @brief Size Constructor
	inline
	explicit
	FArray1( size_type const size_a ) :
		Super( size_a )
	{}


	/// @brief Size + InitializerSentinel Constructor
	inline
	FArray1( size_type const size_a, InitializerSentinel const & initialized ) :
		Super( size_a, initialized )
	{}


	/// @brief Default Proxy Constructor
	inline
	FArray1( ProxySentinel const & proxy ) :
		Super( proxy )
	{}


	/// @brief Copy Proxy Constructor
	inline
	FArray1( FArray1 const & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Non-Const Copy Proxy Constructor
	inline
	FArray1( FArray1 & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Base Proxy Constructor
	inline
	FArray1( Base const & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Non-Const Base Proxy Constructor
	inline
	FArray1( Base & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Section Proxy Constructor
	inline
	FArray1( Section const & s, ProxySentinel const & proxy ) :
		Super( s, proxy )
	{}


	/// @brief Non-Const Section Proxy Constructor
	inline
	FArray1( Section & s, ProxySentinel const & proxy ) :
		Super( s, proxy )
	{}


	/// @brief Value Proxy Constructor
	inline
	FArray1( T const & t, ProxySentinel const & proxy ) :
		Super( t, proxy )
	{}


	/// @brief Non-Const Value Proxy Constructor
	inline
	FArray1( T & t, ProxySentinel const & proxy ) :
		Super( t, proxy )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~FArray1()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray1 &
	operator =( FArray1 const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension_assign( a.I() );
			Super::operator =( a );
		}
		return *this;
	}


	/// @brief Copy Assignment Template
	template< typename U >
	inline
	FArray1 &
	operator =( FArray1< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension_assign( a.I() );
		Super::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray1 &
	operator +=( FArray1< U > const & a )
	{
		assert( equal_dimension( a ) );
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray1 &
	operator -=( FArray1< U > const & a )
	{
		assert( equal_dimension( a ) );
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray1 &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray1 &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray1 &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray1 &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray1 &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i ) const
	inline
	T const &
	operator ()( int const i ) const
	{
		assert( contains( i ) );
		return sarray_[ i ];
	}


	/// @brief array( i )
	inline
	T &
	operator ()( int const i )
	{
		proxy_const_assert( not_const_proxy() );
		assert( contains( i ) );
		return sarray_[ i ];
	}


	/// @brief Const Section Starting at array( i )
	inline
	Section const
	a( int const i ) const
	{
		assert( contains( i ) );
		return Section( static_cast< T const * >( sarray_ + i ), ( array_size_ != npos ) ? array_size_ - ( i - shift_ ) : npos );
	}


	/// @brief Section Starting at array( i )
	inline
	Section
	a( int const i )
	{
		proxy_const_assert( not_const_proxy() );
		assert( contains( i ) );
		return Section( sarray_ + i, ( array_size_ != npos ) ? array_size_ - ( i - shift_ ) : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i ) const
	{
		assert( dimensions_initialized() );
		return ( i - shift_ );
	}


public: // Predicate


	/// @brief Contains Indexed Element?
	virtual
	inline
	bool
	contains( int const i ) const
	{
		return ( I().contains( i ) );
	}


	/// @brief Equal Dimension?
	template< typename U >
	inline
	bool
	equal_dimension( FArray1< U > const & a ) const
	{
		return ( I() == a.I() );
	}


public: // Inspector


	/// @brief IndexRange
	virtual
	IR const &
	I1() const = 0;


	/// @brief Lower Index
	virtual
	int
	l1() const = 0;


	/// @brief Upper Index
	virtual
	int
	u1() const = 0;


	/// @brief Size
	virtual
	size_type
	size1() const = 0;


	/// @brief IndexRange
	virtual
	IR const &
	I() const = 0;


	/// @brief Lower Index
	virtual
	int
	l() const = 0;


	/// @brief Upper Index
	virtual
	int
	u() const = 0;


	/// @brief Length
	inline
	T
	length() const
	{
		T length_sq( T( 0 ) );
		for ( int i = l(), e = u(); i <= e; ++i ) {
			T const length_i( sarray_[ i ] );
			length_sq += length_i * length_i;
		}
		return std::sqrt( length_sq );
	}


	/// @brief Length Squared
	inline
	T
	length_squared() const
	{
		T length_sq( T( 0 ) );
		for ( int i = l(), e = u(); i <= e; ++i ) {
			T const length_i( sarray_[ i ] );
			length_sq += length_i * length_i;
		}
		return length_sq;
	}


public: // Modifier


	/// @brief Clear
	inline
	FArray1 &
	clear()
	{
		Super::clear();
		return *this;
	}


	/// @brief Assign Default Value to all Elements
	inline
	FArray1 &
	to_default()
	{
		Super::to_default();
		return *this;
	}


	/// @brief Normalize to Unit Length
	inline
	FArray1 &
	normalize()
	{
		proxy_const_assert( not_const_proxy() );
		T const length_( length() );
		assert( length_ > T( 0 ) );
		operator /=( length_ );
		return *this;
	}


public: // Comparison


	/// @brief FArray1 == FArray1
	friend
	inline
	bool
	operator ==( FArray1 const & a, FArray1 const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( equal_dimensions( a, b ) ) { // Index ranges match
			return ( static_cast< Super const & >( a ) == static_cast< Super const & >( b ) );
		} else { // Index ranges differ
			return false;
		}
	}


	/// @brief FArray1 != FArray1
	friend
	inline
	bool
	operator !=( FArray1 const & a, FArray1 const & b )
	{
		return !( a == b );
	}


public: // Friend


	/// @brief Dot Product
	template< typename U >
	friend
	U
	dot_product( FArray1<U> const & a, FArray1<U> const & b );


	/// @brief Dot Product
	template< typename U >
	friend
	U
	dot( FArray1<U> const & a, FArray1<U> const & b );


	/// @brief Distance
	template< typename U >
	friend
	U
	distance( FArray1<U> const & a, FArray1<U> const & b );


	/// @brief Distance Squared
	template< typename U >
	friend
	U
	distance_squared( FArray1<U> const & a, FArray1<U> const & b );

protected: // Functions


	/// @brief Dimension by IndexRanges
	virtual
	void
	dimension_assign( IR const & I_a ) = 0;


	/// @brief Swap
	inline
	void
	swap1DB( FArray1 & v )
	{
		this->swapB( v );
	}


}; // FArray1


/// @brief FArray1 == FArray1
template< typename T >
bool
operator ==( FArray1< T > const & a, FArray1< T > const & b );


/// @brief FArray1 != FArray1
template< typename T >
bool
operator !=( FArray1< T > const & a, FArray1< T > const & b );


/// @brief Dot Product
template< typename T >
T
dot_product( FArray1< T > const & a, FArray1< T > const & b );


/// @brief Dot Product
template< typename T >
T
dot( FArray1< T > const & a, FArray1< T > const & b );


/// @brief Distance
template< typename T >
T
distance( FArray1< T > const & a, FArray1< T > const & b );


/// @brief Distance Squared
template< typename T >
T
distance_squared( FArray1< T > const & a, FArray1< T > const & b );


/// @brief Equal Dimensions?
template< typename U, typename V >
inline
bool
equal_dimensions( FArray1< U > const & a, FArray1< V > const & b )
{
	return ( a.I() == b.I() );
}


/// @brief Dot Product
template< typename U >
U
dot_product( FArray1<U> const & a, FArray1<U> const & b )
{
	assert( equal_dimensions( a, b ) );
	U sum( U( 0 ) );
	for ( int i = a.l(), e = a.u(); i <= e; ++i ) {
		sum += a( i ) * b( i );
	}
	return sum;
}


/// @brief Dot Product
template< typename U >
U
dot( FArray1<U> const & a, FArray1<U> const & b )
{
	assert( equal_dimensions( a, b ) );
	U sum( U( 0 ) );
	for ( int i = a.l(), e = a.u(); i <= e; ++i ) {
		sum += a( i ) * b( i );
	}
	return sum;
}


/// @brief Distance
template< typename U >
U
distance( FArray1<U> const & a, FArray1<U> const & b )
{
	assert( equal_dimensions( a, b ) );
	U distance_sq( U( 0 ) );
	for ( int i = a.l(), e = a.u(); i <= e; ++i ) {
		U const distance_i( a( i ) - b( i ) );
		distance_sq += distance_i * distance_i;
	}
	return std::sqrt( distance_sq );
}


/// @brief Distance Squared
template< typename U >
U
distance_squared( FArray1<U> const & a, FArray1<U> const & b )
{
	assert( equal_dimensions( a, b ) );
	U distance_sq( U( 0 ) );
	for ( int i = a.l(), e = a.u(); i <= e; ++i ) {
		U const distance_i( a( i ) - b( i ) );
		distance_sq += distance_i * distance_i;
	}
	return distance_sq;
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1_HH
