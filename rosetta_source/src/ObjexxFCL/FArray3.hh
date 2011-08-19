#ifndef INCLUDED_ObjexxFCL_FArray3_hh
#define INCLUDED_ObjexxFCL_FArray3_hh


// FArray3: Fortran-Compatible 3D Array Abstract Base Class
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
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray.hh>


namespace ObjexxFCL {


// Forward Declarations
template< typename > class FArray3D;
template< typename > class FArray3P;
template< typename > class FArray3A;
template< typename > class KeyFArray3D;


/// @brief FArray3: Fortran-Compatible 3D Array Abstract Base Class
template< typename T >
class FArray3 :
	public FArray< T >
{


private: // Types


	typedef  FArray< T >  Super;
	typedef  FArray3D< T >  real_FArray;
	typedef  FArray3P< T >  proxy_FArray;
	typedef  FArray3A< T >  arg_FArray;


private: // Friend


	template< typename > friend class FArray3;
	template< typename > friend class FArray3D;
	template< typename > friend class FArray3P;
	template< typename > friend class FArray3A;
	template< typename > friend class KeyFArray3D;


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


protected: // Creation


	/// @brief Default Constructor
	inline
	FArray3() :
		s1_( 0 ),
		s2_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	FArray3( FArray3 const & a ) :
		Super( a ),
		s1_( a.s1_ ),
		s2_( a.s2_ )
	{}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	FArray3( FArray3< U > const & a ) :
		Super( a ),
		s1_( a.s1_ ),
		s2_( a.s2_ )
	{}


	/// @brief Size Constructor
	inline
	explicit
	FArray3( size_type const size_a ) :
		Super( size_a )
	{}


	/// @brief Size + InitializerSentinel Constructor
	inline
	FArray3( size_type const size_a, InitializerSentinel const & initialized ) :
		Super( size_a, initialized )
	{}


	/// @brief Default Proxy Constructor
	inline
	FArray3( ProxySentinel const & proxy ) :
		Super( proxy ),
		s1_( 0 ),
		s2_( 0 )
	{}


	/// @brief Copy Proxy Constructor
	inline
	FArray3( FArray3 const & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Non-Const Copy Proxy Constructor
	inline
	FArray3( FArray3 & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Base Proxy Constructor
	inline
	FArray3( Base const & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Non-Const Base Proxy Constructor
	inline
	FArray3( Base & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Section Proxy Constructor
	inline
	FArray3( Section const & s, ProxySentinel const & proxy ) :
		Super( s, proxy )
	{}


	/// @brief Non-Const Section Proxy Constructor
	inline
	FArray3( Section & s, ProxySentinel const & proxy ) :
		Super( s, proxy )
	{}


	/// @brief Value Proxy Constructor
	inline
	FArray3( T const & t, ProxySentinel const & proxy ) :
		Super( t, proxy )
	{}


	/// @brief Non-Const Value Proxy Constructor
	inline
	FArray3( T & t, ProxySentinel const & proxy ) :
		Super( t, proxy )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~FArray3()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray3 &
	operator =( FArray3 const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension_assign( a.I1(), a.I2(), a.I3() );
			Super::operator =( a );
		}
		return *this;
	}


	/// @brief Copy Assignment Template
	template< typename U >
	inline
	FArray3 &
	operator =( FArray3< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension_assign( a.I1(), a.I2(), a.I3() );
		Super::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray3 &
	operator +=( FArray3< U > const & a )
	{
		assert( equal_dimension( a ) );
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray3 &
	operator -=( FArray3< U > const & a )
	{
		assert( equal_dimension( a ) );
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray3 &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray3 &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray3 &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray3 &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray3 &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2, i3 ) const
	inline
	T const &
	operator ()( int const i1, int const i2, int const i3 ) const
	{
		assert( contains( i1, i2, i3 ) );
		return sarray_[ ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2, i3 )
	inline
	T &
	operator ()( int const i1, int const i2, int const i3 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( contains( i1, i2, i3 ) );
		return sarray_[ ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief Const Section Starting at array( i1, i2, i3 )
	inline
	Section const
	a( int const i1, int const i2, int const i3 ) const
	{
		assert( contains( i1, i2, i3 ) );
		size_type const offset( ( ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( static_cast< T const * >( array_ + offset ), ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Section Starting at array( i1, i2, i3 )
	inline
	Section
	a( int const i1, int const i2, int const i3 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( contains( i1, i2, i3 ) );
		size_type const offset( ( ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( array_ + offset, ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i1, int const i2, int const i3 ) const
	{
		assert( dimensions_initialized() );
		return ( ( ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
	}


public: // Predicate


	/// @brief Contains Indexed Element?
	virtual
	inline
	bool
	contains( int const i1, int const i2, int const i3 ) const
	{
		return ( ( I1().contains( i1 ) ) && ( I2().contains( i2 ) ) && ( I3().contains( i3 ) ) );
	}


	/// @brief Equal Dimension?
	template< typename U >
	inline
	bool
	equal_dimension( FArray3< U > const & a ) const
	{
		return ( ( I1() == a.I1() ) && ( I2() == a.I2() ) && ( I3() == a.I3() ) );
	}


public: // Inspector


	/// @brief IndexRange of Dimension 1
	virtual
	IR const &
	I1() const = 0;


	/// @brief Lower Index of Dimension 1
	virtual
	int
	l1() const = 0;


	/// @brief Upper Index of Dimension 1
	virtual
	int
	u1() const = 0;


	/// @brief Size of Dimension 1
	inline
	size_type
	size1() const
	{
		return s1_;
	}


	/// @brief IndexRange of Dimension 2
	virtual
	IR const &
	I2() const = 0;


	/// @brief Lower Index of Dimension 2
	virtual
	int
	l2() const = 0;


	/// @brief Upper Index of Dimension 2
	virtual
	int
	u2() const = 0;


	/// @brief Size of Dimension 2
	inline
	size_type
	size2() const
	{
		return s2_;
	}


	/// @brief IndexRange of Dimension 3
	virtual
	IR const &
	I3() const = 0;


	/// @brief Lower Index of Dimension 3
	virtual
	int
	l3() const = 0;


	/// @brief Upper Index of Dimension 3
	virtual
	int
	u3() const = 0;


	/// @brief Size of Dimension 3
	virtual
	size_type
	size3() const = 0;


public: // Modifier


	/// @brief Clear
	inline
	FArray3 &
	clear()
	{
		Super::clear();
		s1_ = s2_ = 0;
		return *this;
	}


	/// @brief Assign Default Value to all Elements
	inline
	FArray3 &
	to_default()
	{
		Super::to_default();
		return *this;
	}


public: // Comparison


	/// @brief FArray3 == FArray3
	friend
	inline
	bool
	operator ==( FArray3 const & a, FArray3 const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( equal_dimensions( a, b ) ) { // Index ranges match
			return ( static_cast< Super const & >( a ) == static_cast< Super const & >( b ) );
		} else { // Index ranges differ
			return false;
		}
	}


	/// @brief FArray3 != FArray3
	friend
	inline
	bool
	operator !=( FArray3 const & a, FArray3 const & b )
	{
		return !( a == b );
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	virtual
	void
	dimension_assign( IR const & I1_a, IR const & I2_a, IR const & I3_a ) = 0;


	/// @brief Swap
	inline
	void
	swap3DB( FArray3 & v )
	{
		this->swapB( v );
		std::swap( s1_, v.s1_ );
		std::swap( s2_, v.s2_ );
	}


protected: // Data


	/// @brief Dimension 1 size
	size_type s1_;

	/// @brief Dimension 2 size
	size_type s2_;


}; // FArray3


/// @brief FArray3 == FArray3
template< typename T >
bool
operator ==( FArray3< T > const & a, FArray3< T > const & b );


/// @brief FArray3 != FArray3
template< typename T >
bool
operator !=( FArray3< T > const & a, FArray3< T > const & b );


/// @brief Equal Dimensions?
template< typename U, typename V >
inline
bool
equal_dimensions( FArray3< U > const & a, FArray3< V > const & b )
{
	return ( ( a.I1() == b.I1() ) && ( a.I2() == b.I2() ) && ( a.I3() == b.I3() ) );
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray3_HH
