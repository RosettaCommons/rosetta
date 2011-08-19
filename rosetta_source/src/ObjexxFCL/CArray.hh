#ifndef INCLUDED_ObjexxFCL_CArray_hh
#define INCLUDED_ObjexxFCL_CArray_hh


// CArray: Memory-Managed C Array Wrapper
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
#include <ObjexxFCL/CArray.fwd.hh>

// C++ Headers
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>


namespace ObjexxFCL {


/// @brief CArray: Memory-Managed C Array Wrapper
template< typename T >
class CArray
{


private: // Friend


	template< typename > friend class CArray; // Friendship across value types


public: // Types


	// STL style
	typedef  T  value_type;
	typedef  T &  reference;
	typedef  T const &  const_reference;
	typedef  T *  pointer;
	typedef  T const *  const_pointer;
	typedef  T *  iterator;
	typedef  T const *  const_iterator;
	typedef  std::size_t  size_type;
	typedef  std::ptrdiff_t  difference_type;

	// C++ style
	typedef  T  Value;
	typedef  T &  Reference;
	typedef  T const &  ConstReference;
	typedef  T *  Pointer;
	typedef  T const *  ConstPointer;
	typedef  T *  Iterator;
	typedef  T const *  ConstIterator;
	typedef  std::size_t  Size;
	typedef  std::ptrdiff_t  Difference;

	// Types to prevent compile failure when std::distance is in scope
	typedef  void  iterator_category;


public: // Creation


	/// @brief Default constructor
	inline
	CArray() :
		size_( 0 ),
		array_( 0 )
	{}


	/// @brief Copy constructor
	inline
	CArray( CArray const & a ) :
		size_( a.size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = a.array_[ i ];
		}
	}


	/// @brief Copy constructor template
	template< typename U >
	inline
	CArray( CArray< U > const & a ) :
		size_( a.size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( a.array_[ i ] );
		}
	}


	/// @brief Pointer + size constructor
	inline
	CArray(
		T const * const p,
		size_type const size_a
	) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = p[ i ];
		}
	}


	/// @brief Pointer + size constructor template
	template< typename U >
	inline
	CArray(
		U const * const p,
		size_type const size_a
	) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( p[ i ] );
		}
	}


	/// @brief Iterator range constructor template
	template< typename InputIterator >
	inline
	CArray(
		InputIterator const beg,
		InputIterator const end
	) :
		size_( end - beg ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		if ( size_ > 0 ) {
			InputIterator k( beg );
			for ( size_type i = 0; i < size_; ++i, ++k ) {
				array_[ i ] = T( *k );
			}
		}
	}


	/// @brief Size constructor
	/// @note Built-in value types are not initialized
	inline
	explicit
	CArray( size_type const size_a ) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{}


	/// @brief Size + uniform value constructor
	inline
	CArray(
		size_type const size_a,
		T const & t
	) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = t;
		}
	}


	/// @brief Destructor
	inline
	~CArray()
	{
		delete[] array_;
	}


public: // Conversion


	/// @brief Active?
	inline
	operator bool() const
	{
		return ( array_ != 0 );
	}


public: // Assignment


	/// @brief Copy assignment
	inline
	CArray &
	operator =( CArray const & a )
	{
		if ( this != &a ) {
			if ( size_ != a.size_ ) {
				size_ = a.size_;
				delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
			}
			for ( size_type i = 0; i < size_; ++i ) {
				array_[ i ] = a.array_[ i ];
			}
		}
		return *this;
	}


	/// @brief Copy assignment template
	template< typename U >
	inline
	CArray &
	operator =( CArray< U > const & a )
	{
		if ( size_ != a.size_ ) {
			size_ = a.size_;
			delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
		}
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( a.array_[ i ] );
		}
		return *this;
	}


	/// @brief Uniform value assignment
	inline
	CArray &
	operator =( T const & t )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = t;
		}
		return *this;
	}


	/// @brief Pointer + size assignment
	inline
	CArray &
	assign(
		T const * const p,
		size_type const size_a
	)
	{
		if ( size_ != size_a ) {
			size_ = size_a;
			delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
		}
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = p[ i ];
		}
		return *this;
	}


	/// @brief Pointer + size assignment template
	template< typename U >
	inline
	CArray &
	assign(
		U const * const p,
		size_type const size_a
	)
	{
		if ( size_ != size_a ) {
			size_ = size_a;
			delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
		}
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( p[ i ] );
		}
		return *this;
	}


	/// @brief Iterator range assignment template
	template< typename InputIterator >
	inline
	CArray &
	assign(
		InputIterator const beg,
		InputIterator const end
	)
	{
		size_type const size_a( end - beg );
		if ( size_ != size_a ) {
			size_ = size_a;
			delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
		}
		if ( size_ > 0 ) {
			InputIterator k( beg );
			for ( size_type i = 0; i < size_; ++i, ++k ) {
				array_[ i ] = T( *k );
			}
		}
		return *this;
	}


	/// @brief Size + value assignment
	inline
	CArray &
	assign(
		size_type const size_a,
		T const & value
	)
	{
		if ( size_ != size_a ) { // Set to new array with uniform values
			CArray( size_a, value ).swap( *this );
		} else { // Set to uniform value
			(*this) = value;
		}
		return *this;
	}


	/// @brief += CArray
	template< typename U >
	inline
	CArray &
	operator +=( CArray< U > const & a )
	{
		assert( size_ == a.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += T( a.array_[ i ] );
		}
		return *this;
	}


	/// @brief -= CArray
	template< typename U >
	inline
	CArray &
	operator -=( CArray< U > const & a )
	{
		assert( size_ == a.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= T( a.array_[ i ] );
		}
		return *this;
	}


	/// @brief += value
	inline
	CArray &
	operator +=( T const & t )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += t;
		}
		return *this;
	}


	/// @brief -= value
	inline
	CArray &
	operator -=( T const & t )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= t;
		}
		return *this;
	}


	/// @brief *= value
	inline
	CArray &
	operator *=( T const & t )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] *= t;
		}
		return *this;
	}


	/// @brief /= value
	inline
	CArray &
	operator /=( T const & t )
	{
		assert( t != T( 0 ) );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] /= t;
		}
		return *this;
	}


public: // Predicate


	/// @brief Active?
	inline
	bool
	active() const
	{
		return ( array_ != 0 );
	}


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return ( size_ == 0 );
	}


public: // Inspector


	/// @brief Size
	inline
	size_type
	size() const
	{
		return size_;
	}


	/// @brief First element
	inline
	T const &
	front() const
	{
		assert( size_ > 0 );
		return array_[ 0 ];
	}


	/// @brief Last element
	inline
	T const &
	back() const
	{
		assert( size_ > 0 );
		return array_[ size_ - 1 ];
	}


	/// @brief Length
	inline
	T
	length() const
	{
		T length_sq( T( 0 ) );
		for ( size_type i = 0; i < size_; ++i ) {
			length_sq += square( array_[ i ] );
		}
		return std::sqrt( length_sq );
	}


	/// @brief Length squared
	inline
	T
	length_squared() const
	{
		T length_sq( T( 0 ) );
		for ( size_type i = 0; i < size_; ++i ) {
			length_sq += square( array_[ i ] );
		}
		return length_sq;
	}


public: // Modifier


	/// @brief First element
	inline
	T &
	front()
	{
		assert( size_ > 0 );
		return array_[ 0 ];
	}


	/// @brief Last element
	inline
	T &
	back()
	{
		assert( size_ > 0 );
		return array_[ size_ - 1 ];
	}


	/// @brief Resize: Values not preserved
	/// @note Built-in values are uninitialized if size changes
	inline
	CArray &
	size( size_type const size_a )
	{
		if ( size_ != size_a ) { // Set to new array
			CArray( size_a ).swap( *this );
		}
		return *this;
	}


	/// @brief Resize to size with fill value: Values preserved
	inline
	CArray &
	resize(
		size_type const size_a,
		T const & fill = T()
	)
	{
		if ( size_ < size_a ) {
			CArray a( size_a, fill ); // New array: Elements set to fill fill
			for ( size_type i = 0; i < size_; ++i ) { // Copy current values
				a.array_[ i ] = array_[ i ];
			}
			swap( a ); // Swap in new array
		} else if ( size_ > size_a ) {
			CArray a( size_a ); // New array
			for ( size_type i = 0; i < size_a; ++i ) { // Copy current values within new range
				a.array_[ i ] = array_[ i ];
			}
			swap( a ); // Swap in new array
		}
		return *this;
	}


	/// @brief Swap
	inline
	void
	swap( CArray & a )
	{
		std::swap( size_, a.size_ );
		std::swap( array_, a.array_ );
	}


	/// @brief Clear
	inline
	CArray &
	clear()
	{
		size_ = 0;
		delete[] array_; array_ = 0;
		return *this;
	}


	/// @brief Normalize to unit length
	inline
	CArray &
	normalize()
	{
		T const length_( length() );
		assert( length_ > T( 0 ) );
		operator /=( length_ );
		return *this;
	}


public: // Subscript


	/// @brief CArray[ i ] const: 0-based indexing
	inline
	T const &
	operator []( size_type const i ) const
	{
		assert( i < size_ );
		return array_[ i ];
	}


	/// @brief CArray[ i ]: 0-based indexing
	inline
	T &
	operator []( size_type const i )
	{
		assert( i < size_ );
		return array_[ i ];
	}


#ifdef OBJEXXFCL_CARRAY_1_BASED_LOOKUP


	/// @brief CArray( i ) const: 1-based indexing
	inline
	T const &
	operator ()( size_type const i ) const
	{
		assert( ( i > 0 ) && ( i <= size_ ) );
		return array_[ i - 1 ];
	}


	/// @brief CArray( i ): 1-based indexing
	inline
	T &
	operator ()( size_type const i )
	{
		assert( ( i > 0 ) && ( i <= size_ ) );
		return array_[ i - 1 ];
	}


#endif // OBJEXXFCL_CARRAY_1_BASED_LOOKUP


public: // Iterator


	/// @brief const_iterator to beginning of array
	inline
	const_iterator
	begin() const
	{
		return array_;
	}


	/// @brief iterator to beginning of array
	inline
	iterator
	begin()
	{
		return array_;
	}


	/// @brief const_iterator to element past end of array
	inline
	const_iterator
	end() const
	{
		return array_ + size_;
	}


	/// @brief iterator to element past end of array
	inline
	iterator
	end()
	{
		return array_ + size_;
	}


public: // Array Accessor


	/// @brief C array const accessor
	inline
	T const &
	operator ()() const
	{
		return array_;
	}


	/// @brief C array non-const accessor
	inline
	T &
	operator ()()
	{
		return array_;
	}


public: // Comparison


	/// @brief Are two CArrays comparable?
	friend
	inline
	bool
	comparable( CArray const & a, CArray const & b )
	{
		return ( a.size_ == b.size_ );
	}


	/// @brief CArray == CArray
	friend
	inline
	bool
	operator ==( CArray const & a, CArray const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a.array_[ i ] == b.array_[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief CArray != CArray
	friend
	inline
	bool
	operator !=( CArray const & a, CArray const & b )
	{
		return !( a == b );
	}


	/// @brief CArray < CArray
	friend
	inline
	bool
	operator <( CArray const & a, CArray const & b )
	{
		if ( &a == &b ) { // Same objects
			return false;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a.array_[ i ] < b.array_[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief CArray <= CArray
	friend
	inline
	bool
	operator <=( CArray const & a, CArray const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a.array_[ i ] <= b.array_[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief CArray >= CArray
	friend
	inline
	bool
	operator >=( CArray const & a, CArray const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a.array_[ i ] >= b.array_[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief CArray > CArray
	friend
	inline
	bool
	operator >( CArray const & a, CArray const & b )
	{
		if ( &a == &b ) { // Same objects
			return false;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a.array_[ i ] > b.array_[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief CArray == T
	friend
	inline
	bool
	operator ==( CArray const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( a.array_[ i ] != t ) return false;
		}
		return true;
	}


	/// @brief CArray != T
	friend
	inline
	bool
	operator !=( CArray const & a, T const & t )
	{
		return !( a == t );
	}


	/// @brief CArray < T
	friend
	inline
	bool
	operator <( CArray const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] < t ) ) return false;
		}
		return true;
	}


	/// @brief CArray <= T
	friend
	inline
	bool
	operator <=( CArray const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] <= t ) ) return false;
		}
		return true;
	}


	/// @brief CArray >= T
	friend
	inline
	bool
	operator >=( CArray const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] >= t ) ) return false;
		}
		return true;
	}


	/// @brief CArray > T
	friend
	inline
	bool
	operator >( CArray const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] > t ) ) return false;
		}
		return true;
	}


	/// @brief T == CArray
	friend
	inline
	bool
	operator ==( T const & t, CArray const & a )
	{
		return ( a == t );
	}


	/// @brief T != CArray
	friend
	inline
	bool
	operator !=( T const & t, CArray const & a )
	{
		return !( t == a );
	}


	/// @brief T < CArray
	friend
	inline
	bool
	operator <( T const & t, CArray const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t < a.array_[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T <= CArray
	friend
	inline
	bool
	operator <=( T const & t, CArray const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t <= a.array_[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T >= CArray
	friend
	inline
	bool
	operator >=( T const & t, CArray const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t >= a.array_[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T > CArray
	friend
	inline
	bool
	operator >( T const & t, CArray const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t > a.array_[ i ] ) ) return false;
		}
		return true;
	}


public: // Generator


	/// @brief -CArray
	friend
	inline
	CArray
	operator -( CArray const & a )
	{
		CArray r( a );
		r *= T( -1 );
		return r;
	}


	/// @brief CArray + CArray
	friend
	inline
	CArray
	operator +( CArray const & a, CArray const & b )
	{
		CArray r( a );
		r += b;
		return r;
	}


	/// @brief CArray - CArray
	friend
	inline
	CArray
	operator -( CArray const & a, CArray const & b )
	{
		CArray r( a );
		r -= b;
		return r;
	}


	/// @brief CArray + Value
	friend
	inline
	CArray
	operator +( CArray const & a, T const & t )
	{
		CArray r( a );
		r += t;
		return r;
	}


	/// @brief Value + CArray
	friend
	inline
	CArray
	operator +( T const & t, CArray const & a )
	{
		CArray r( a );
		r += t;
		return r;
	}


	/// @brief CArray - Value
	friend
	inline
	CArray
	operator -( CArray const & a, T const & t )
	{
		CArray r( a );
		r -= t;
		return r;
	}


	/// @brief Value - CArray
	friend
	inline
	CArray
	operator -( T const & t, CArray const & a )
	{
		CArray r( -a );
		r += t;
		return r;
	}


	/// @brief CArray * Value
	friend
	inline
	CArray
	operator *( CArray const & a, T const & t )
	{
		CArray r( a );
		r *= t;
		return r;
	}


	/// @brief Value * CArray
	friend
	inline
	CArray
	operator *( T const & t, CArray const & a )
	{
		CArray r( a );
		r *= t;
		return r;
	}


	/// @brief CArray / Value
	friend
	inline
	CArray
	operator /( CArray const & a, T const & t )
	{
		CArray r( a );
		r /= t;
		return r;
	}


public: // Friend


	/// @brief Dot product
	friend
	inline
	T
	dot_product( CArray const & a, CArray const & b )
	{
		assert( a.size() == b.size() );
		T sum( T( 0 ) );
		for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
			sum += a[ i ] * b[ i ];
		}
		return sum;
	}


	/// @brief Dot product
	friend
	inline
	T
	dot( CArray const & a, CArray const & b )
	{
		assert( a.size() == b.size() );
		T sum( T( 0 ) );
		for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
			sum += a[ i ] * b[ i ];
		}
		return sum;
	}


	/// @brief Distance
	friend
	inline
	T
	distance( CArray const & a, CArray const & b )
	{
		assert( a.size() == b.size() );
		T distance_sq( T( 0 ) );
		for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
			distance_sq += square( a[ i ] - b[ i ] );
		}
		return std::sqrt( distance_sq );
	}


	/// @brief Distance squared
	friend
	inline
	T
	distance_squared( CArray const & a, CArray const & b )
	{
		assert( a.size() == b.size() );
		T distance_sq( T( 0 ) );
		for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
			distance_sq += square( a[ i ] - b[ i ] );
		}
		return distance_sq;
	}


	/// @brief Swap
	friend
	inline
	void
	swap( CArray & a, CArray & b )
	{
		a.swap( b );
	}


private: // Static Functions


	/// @brief square( x ) == x^2
	inline
	static
	T
	square( T const & x )
	{
		return x * x;
	}


private: // Data


	/// @brief Number of array elements
	size_type size_;

	/// @brief C array
	T * array_;


}; // CArray


/// @brief Are two CArrays comparable?
template< typename T >
bool
comparable( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > == CArray< T >
template< typename T >
bool
operator ==( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > != CArray< T >
template< typename T >
bool
operator !=( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > < CArray< T >
template< typename T >
bool
operator <( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > <= CArray< T >
template< typename T >
bool
operator <=( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > >= CArray< T >
template< typename T >
bool
operator >=( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > > CArray< T >
template< typename T >
bool
operator >( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > == T
template< typename T >
bool
operator ==( CArray< T > const & a, T const & t );


/// @brief CArray< T > != T
template< typename T >
bool
operator !=( CArray< T > const & a, T const & t );


/// @brief CArray< T > < T
template< typename T >
bool
operator <( CArray< T > const & a, T const & t );


/// @brief CArray< T > <= T
template< typename T >
bool
operator <=( CArray< T > const & a, T const & t );


/// @brief CArray< T > >= T
template< typename T >
bool
operator >=( CArray< T > const & a, T const & t );


/// @brief CArray< T > > T
template< typename T >
bool
operator >( CArray< T > const & a, T const & t );


/// @brief T == CArray< T >
template< typename T >
bool
operator ==( T const & t, CArray< T > const & a );


/// @brief T != CArray< T >
template< typename T >
bool
operator !=( T const & t, CArray< T > const & a );


/// @brief T < CArray< T >
template< typename T >
bool
operator <( T const & t, CArray< T > const & a );


/// @brief T <= CArray< T >
template< typename T >
bool
operator <=( T const & t, CArray< T > const & a );


/// @brief T >= CArray< T >
template< typename T >
bool
operator >=( T const & t, CArray< T > const & a );


/// @brief T > CArray< T >
template< typename T >
bool
operator >( T const & t, CArray< T > const & a );


/// @brief -CArray< T >
template< typename T >
CArray< T >
operator -( CArray< T > const & a );


/// @brief CArray< T > + CArray< T >
template< typename T >
CArray< T >
operator +( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > - CArray< T >
template< typename T >
CArray< T >
operator -( CArray< T > const & a, CArray< T > const & b );


/// @brief CArray< T > + Value
template< typename T >
CArray< T >
operator +( CArray< T > const & a, T const & t );


/// @brief Value + CArray< T >
template< typename T >
CArray< T >
operator +( T const & t, CArray< T > const & a );


/// @brief CArray< T > - Value
template< typename T >
CArray< T >
operator -( CArray< T > const & a, T const & t );


/// @brief Value - CArray< T >
template< typename T >
CArray< T >
operator -( T const & t, CArray< T > const & a );


/// @brief CArray< T > * Value
template< typename T >
CArray< T >
operator *( CArray< T > const & a, T const & t );


/// @brief Value * CArray< T >
template< typename T >
CArray< T >
operator *( T const & t, CArray< T > const & a );


/// @brief CArray< T > / Value
template< typename T >
CArray< T >
operator /( CArray< T > const & a, T const & t );


/// @brief Dot product
template< typename T >
T
dot_product( CArray< T > const & a, CArray< T > const & b );


/// @brief Dot product
template< typename T >
T
dot( CArray< T > const & a, CArray< T > const & b );


/// @brief Distance
template< typename T >
T
distance( CArray< T > const & a, CArray< T > const & b );


/// @brief Distance squared
template< typename T >
T
distance_squared( CArray< T > const & a, CArray< T > const & b );


/// @brief Swap
template< typename T >
void
swap( CArray< T > & a, CArray< T > & b );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.  The legal alternative would be
// to add specializations of swap for each anticipated instantiation.


namespace std {


/// @brief std::swap( CArray, CArray )
template< typename T >
inline
void
swap( ObjexxFCL::CArray< T > & a, ObjexxFCL::CArray< T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_CArray_HH
