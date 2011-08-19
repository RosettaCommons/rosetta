#ifndef INCLUDED_ObjexxFCL_FArray2_hh
#define INCLUDED_ObjexxFCL_FArray2_hh


// FArray2: Fortran-Compatible 2D Array Abstract Base Class
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
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray.hh>


namespace ObjexxFCL {


// Forward Declarations
template< typename > class FArray2D;
template< typename > class FArray2P;
template< typename > class FArray2A;
template< typename > class KeyFArray2D;


/// @brief FArray2: Fortran-Compatible 2D Array Abstract Base Class
template< typename T >
class FArray2 :
	public FArray< T >
{


private: // Types


	typedef  FArray< T >  Super;
	typedef  FArray2D< T >  real_FArray;
	typedef  FArray2P< T >  proxy_FArray;
	typedef  FArray2A< T >  arg_FArray;


private: // Friend


	template< typename > friend class FArray2;
	template< typename > friend class FArray2D;
	template< typename > friend class FArray2P;
	template< typename > friend class FArray2A;
	template< typename > friend class KeyFArray2D;


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
	FArray2() :
		s1_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	FArray2( FArray2 const & a ) :
		Super( a ),
		s1_( a.s1_ )
	{}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	FArray2( FArray2< U > const & a ) :
		Super( a ),
		s1_( a.s1_ )
	{}


	/// @brief Size Constructor
	inline
	explicit
	FArray2( size_type const size_a ) :
		Super( size_a )
	{}


	/// @brief Size + InitializerSentinel Constructor
	inline
	FArray2( size_type const size_a, InitializerSentinel const & initialized ) :
		Super( size_a, initialized )
	{}


	/// @brief Default Proxy Constructor
	inline
	FArray2( ProxySentinel const & proxy ) :
		Super( proxy ),
		s1_( 0 )
	{}


	/// @brief Copy Proxy Constructor
	inline
	FArray2( FArray2 const & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Non-Const Copy Proxy Constructor
	inline
	FArray2( FArray2 & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Base Proxy Constructor
	inline
	FArray2( Base const & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Non-Const Base Proxy Constructor
	inline
	FArray2( Base & a, ProxySentinel const & proxy ) :
		Super( a, proxy )
	{}


	/// @brief Section Proxy Constructor
	inline
	FArray2( Section const & s, ProxySentinel const & proxy ) :
		Super( s, proxy )
	{}


	/// @brief Non-Const Section Proxy Constructor
	inline
	FArray2( Section & s, ProxySentinel const & proxy ) :
		Super( s, proxy )
	{}


	/// @brief Value Proxy Constructor
	inline
	FArray2( T const & t, ProxySentinel const & proxy ) :
		Super( t, proxy )
	{}


	/// @brief Non-Const Value Proxy Constructor
	inline
	FArray2( T & t, ProxySentinel const & proxy ) :
		Super( t, proxy )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~FArray2()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray2 &
	operator =( FArray2 const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension_assign( a.I1(), a.I2() );
			Super::operator =( a );
		}
		return *this;
	}


	/// @brief Copy Assignment Template
	template< typename U >
	inline
	FArray2 &
	operator =( FArray2< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension_assign( a.I1(), a.I2() );
		Super::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray2 &
	operator +=( FArray2< U > const & a )
	{
		assert( equal_dimension( a ) );
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray2 &
	operator -=( FArray2< U > const & a )
	{
		assert( equal_dimension( a ) );
		Super::operator -=( a );
		return *this;
	}


	/// @brief *= Array Template
	template< typename U >
	inline
	FArray2 &
	operator *=( FArray2< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( a.square() );
		assert( a.size2() == size1() );

		T * const row( new T[ size2() ] ); // Temporary row
		int const off( l2() );

		FArray2 & A( *this ); // Shorthand name
		for ( int i = l1(), ie = u1(); i <= ie; ++i ) { // Row i
			for ( int j = l2(), je = u2(); j <= je; ++j ) { // Col j
				row[ j - off ] = A( i, j ); // Save the row
			}
			for ( int j = l2(), je = u2(), jj = a.l2(); j <= je; ++j, ++jj ) { // Col j
				T sum( T( 0 ) );
				for ( int k = 0, ke = size2(), ii = a.l1(); k < ke; ++k, ++ii ) {
					sum += row[ k ] * a( ii, jj );
				}
				A( i, j ) = sum;
			}
		}

		delete[] row;
		return *this;
	}


	/// @brief = Value
	inline
	FArray2 &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray2 &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray2 &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray2 &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray2 &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2 ) const
	inline
	T const &
	operator ()( int const i1, int const i2 ) const
	{
		assert( contains( i1, i2 ) );
		return sarray_[ ( i2 * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2 )
	inline
	T &
	operator ()( int const i1, int const i2 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( contains( i1, i2 ) );
		return sarray_[ ( i2 * s1_ ) + i1 ];
	}


	/// @brief Const Section Starting at array( i1, i2 )
	inline
	Section const
	a( int const i1, int const i2 ) const
	{
		assert( contains( i1, i2 ) );
		size_type const offset( ( ( i2 * s1_ ) + i1 ) - shift_ );
		return Section( static_cast< T const * >( array_ + offset ), ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Section Starting at array( i1, i2 )
	inline
	Section
	a( int const i1, int const i2 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( contains( i1, i2 ) );
		size_type const offset( ( ( i2 * s1_ ) + i1 ) - shift_ );
		return Section( array_ + offset, ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i1, int const i2 ) const
	{
		assert( dimensions_initialized() );
		return ( ( ( i2 * s1_ ) + i1 ) - shift_ );
	}


public: // Predicate


	/// @brief Contains Indexed Element?
	virtual
	inline
	bool
	contains( int const i1, int const i2 ) const
	{
		return ( ( I1().contains( i1 ) ) && ( I2().contains( i2 ) ) );
	}


	/// @brief Equal Dimension?
	template< typename U >
	inline
	bool
	equal_dimension( FArray2< U > const & a ) const
	{
		return ( ( I1() == a.I1() ) && ( I2() == a.I2() ) );
	}


	/// @brief Is Identity?
	inline
	bool
	is_identity() const
	{
		static T const ZERO( 0 );
		static T const ONE( 1 );
		FArray2 const & A( *this ); // Shorthand name
		if ( ! square() ) { // Non-square
			return false;
		} else { // Square
			for ( int j = l2(), je = u2(); j <= je; ++j ) {
				int const id( l1() + ( j - l2() ) ); // Row index of diagonal
				for ( int i = l1(), ie = u1(); i <= ie; ++i ) {
					if ( A( i, j ) != ( i == id ? ONE : ZERO ) ) return false;
				}
			}
			return true;
		}
	}


	/// @brief Square?
	inline
	bool
	square() const
	{
		return ( ( dimensions_initialized() ) && ( s1_ == size2() ) );
	}


	/// @brief Symmetric?
	inline
	bool
	symmetric() const
	{
		FArray2 const & A( *this ); // Shorthand name
		if ( I1() != I2() ) { // Unequal index ranges
			return false;
		} else { // Equal index ranges
			for ( int i = l1(), ie = u1(); i <= ie; ++i ) {
				for ( int j = l2(), je = i - 1; j <= je; ++j ) {
					if ( A( i, j ) != A( j, i ) ) return false;
				}
			}
			return true;
		}
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
	virtual
	size_type
	size2() const = 0;


public: // Modifier


	/// @brief Clear
	inline
	FArray2 &
	clear()
	{
		Super::clear();
		s1_ = 0;
		return *this;
	}


	/// @brief Assign Default Value to all Elements
	inline
	FArray2 &
	to_default()
	{
		Super::to_default();
		return *this;
	}


	/// @brief Set to the Identity Matrix
	inline
	FArray2 &
	to_identity()
	{
		proxy_const_assert( not_const_proxy() );
		assert( square() );
		FArray2 & A( *this ); // Shorthand name
		A = T( 0 ); // Zero the array
		T const One( T( 1 ) );
		for ( int i = l1(), j = l2(), e = u1(); i <= e; ++i, ++j ) { // Set diagonal to unity
			A( i, j ) = One;
		}
		return *this;
	}


	/// @brief Set to Diagonal Matrix with Uniform Value
	inline
	FArray2 &
	to_diag( T const & d )
	{
		proxy_const_assert( not_const_proxy() );
		assert( square() );
		FArray2 & A( *this ); // Shorthand name
		A = T( 0 ); // Zero the array
		for ( int i = l1(), j = l2(), e = u1(); i <= e; ++i, ++j ) { // Set diagonal value
			A( i, j ) = d;
		}
		return *this;
	}


	/// @brief Set Diagonal of Matrix to a Uniform Value
	inline
	FArray2 &
	set_diagonal( T const & d )
	{
		proxy_const_assert( not_const_proxy() );
		assert( square() );
		FArray2 & A( *this ); // Shorthand name
		for ( int i = l1(), j = l2(), e = u1(); i <= e; ++i, ++j ) { // Set diagonal value
			A( i, j ) = d;
		}
		return *this;
	}


	/// @brief Transpose
	inline
	FArray2 &
	transpose()
	{
		proxy_const_assert( not_const_proxy() );
		assert( square() );
		FArray2 & A( *this ); // Shorthand name
		int const l_off( l2() - l1() );
		for ( int i = l1(), ie = u1(); i <= ie; ++i ) {
			for ( int j = l2(), je = i + l_off; j < je; ++j ) {
				T const A_ij( A( i, j ) );
				A( i, j ) = A( j, i );
				A( j, i ) = A_ij;
			}
		}
		return *this;
	}


	/// @brief Right Multiply By Array
	template< typename U >
	inline
	FArray2 &
	right_multiply_by( FArray2< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		return ( operator *=( a ) );
	}


	/// @brief Right Multiply By Transpose of Array
	template< typename U >
	inline
	FArray2 &
	right_multiply_by_transpose( FArray2< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( a.square() );
		assert( a.size2() == size2() );

		T * const row( new T[ size2() ] ); // Temporary row
		int const off( l2() );

		FArray2 & A( *this ); // Shorthand name
		for ( int i = l1(), ie = u1(); i <= ie; ++i ) { // Row i
			for ( int j = l2(), je = u2(); j <= je; ++j ) { // Col j
				row[ j - off ] = A( i, j ); // Save the row
			}
			for ( int j = l2(), je = u2(), ii = a.l1(); j <= je; ++j, ++ii ) { // Col j
				T sum( T( 0 ) );
				for ( int k = 0, ke = size2(), jj = a.l2(); k < ke; ++k, ++jj ) {
					sum += row[ k ] * a( ii, jj );
				}
				A( i, j ) = sum;
			}
		}

		delete[] row;
		return *this;
	}


public: // Comparison


	/// @brief FArray2 == FArray2
	friend
	inline
	bool
	operator ==( FArray2 const & a, FArray2 const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( equal_dimensions( a, b ) ) { // Index ranges match
			return ( static_cast< Super const & >( a ) == static_cast< Super const & >( b ) );
		} else { // Index ranges differ
			return false;
		}
	}


	/// @brief FArray2 != FArray2
	friend
	inline
	bool
	operator !=( FArray2 const & a, FArray2 const & b )
	{
		return !( a == b );
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	virtual
	void
	dimension_assign( IR const & I1_a, IR const & I2_a ) = 0;


	/// @brief Swap
	inline
	void
	swap2DB( FArray2 & v )
	{
		this->swapB( v );
		std::swap( s1_, v.s1_ );
	}


protected: // Data


	/// @brief Dimension 1 size
	size_type s1_;


}; // FArray2


/// @brief FArray2 == FArray2
template< typename T >
bool
operator ==( FArray2< T > const & a, FArray2< T > const & b );


/// @brief FArray2 != FArray2
template< typename T >
bool
operator !=( FArray2< T > const & a, FArray2< T > const & b );


/// @brief Equal Dimensions?
template< typename U, typename V >
inline
bool
equal_dimensions( FArray2< U > const & a, FArray2< V > const & b )
{
	return ( ( a.I1() == b.I1() ) && ( a.I2() == b.I2() ) );
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray2_HH
