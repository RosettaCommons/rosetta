#ifndef INCLUDED_ObjexxFCL_FArray_hh
#define INCLUDED_ObjexxFCL_FArray_hh


// FArray: Fortran-Compatible Array Abstract Base Class
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
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/proxy_const_assert.hh>

// C++ Headers
#include <algorithm>
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
#include <iostream>
#include <typeinfo>
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT


namespace ObjexxFCL {


/// @brief FArray: Fortran-Compatible Array Hierarchy Abstract Base Class
///
/// @remarks
///  @li Can hold numeric or non-numeric values: Numeric operations on non-numeric values won't compile
///  @li Any meaningful array index ranges can be specified as in Fortran
///  @li Column-major storage order used as in Fortran
///  @li Zero-sized arrays are supported but have no valid indices
///  @li Argument/proxy arrays can have unbounded/unknown size
///  @li For efficiency constructors without initializer function or value do not initialize the array
///  @li For efficiency bounds checking is only active via asserts in debug builds
template< typename T >
class FArray
{


private: // Friend


	template< typename > friend class FArray;


protected: // Types


	typedef  internal::InitializerSentinel  InitializerSentinel;
	typedef  internal::ProxySentinel  ProxySentinel;


public: // Types


	typedef  FArray< T >  Base;
	typedef  FArrayTraits< T >  Traits;
	typedef  FArraySection< T >  Section;
	typedef  IndexRange  IR;

	// STL style
	typedef  T  value_type;
	typedef  T &  reference;
	typedef  T const &  const_reference;
	typedef  T *  pointer;
	typedef  T const *  const_pointer;
	typedef  std::size_t  size_type;
	typedef  std::ptrdiff_t  difference_type;

	// C++ style
	typedef  T  Value;
	typedef  T &  Reference;
	typedef  T const &  ConstReference;
	typedef  T *  Pointer;
	typedef  T const *  ConstPointer;
	typedef  std::size_t  Size;
	typedef  std::ptrdiff_t  Difference;


protected: // Creation


	/// @brief Default Constructor
	inline
	FArray() :
		array_size_( 0 ),
		array_( 0 ),
		size_( 0 ),
		owner_( true ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
		size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
	}


	/// @brief Copy Constructor
	inline
	FArray( FArray const & a ) :
		array_size_( size_of( a.size_ ) ),
		array_( array_size_ > 0 ? new T[ array_size_ ] : 0 ),
		size_( array_size_ ),
		owner_( true ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( a.shift_ ),
		sarray_( array_ - shift_ )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = a[ i ];
		}
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
		size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
	}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	FArray( FArray< U > const & a ) :
		array_size_( size_of( a.size() ) ),
		array_( array_size_ > 0 ? new T[ array_size_ ] : 0 ),
		size_( array_size_ ),
		owner_( true ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( a.shift_ ),
		sarray_( array_ - shift_ )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( a[ i ] );
		}
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
		size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
	}


	/// @brief Size Constructor
	inline
	explicit
	FArray( size_type const size_a ) :
		array_size_( size_of( size_a ) ),
		array_( array_size_ > 0 ? new T[ array_size_ ] : 0 ),
		size_( array_size_ ),
		owner_( true ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{
#ifdef OBJEXXFCL_FARRAY_INIT
		std::fill_n( array_, size_, Traits::initial_value() );
#endif // OBJEXXFCL_FARRAY_INIT
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
		size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
	}


	/// @brief Size + InitializerSentinel Constructor
	inline
	FArray( size_type const size_a, InitializerSentinel const & ) :
		array_size_( size_of( size_a ) ),
		array_( array_size_ > 0 ? new T[ array_size_ ] : 0 ),
		size_( array_size_ ),
		owner_( true ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
		size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
	}


	/// @brief Default Proxy Constructor
	inline
	FArray( ProxySentinel const & ) :
		array_size_( 0 ),
		array_( 0 ),
		size_( 0 ),
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


	/// @brief Array Proxy Constructor
	inline
	FArray( FArray const & a, ProxySentinel const & ) :
		array_size_( a.array_size_ ),
		array_( a.array_ ),
		size_( a.size_ ),
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( true ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


	/// @brief Non-Const Array Proxy Constructor
	inline
	FArray( FArray & a, ProxySentinel const & ) :
		array_size_( a.array_size_ ),
		array_( a.array_ ),
		size_( a.size_ ),
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( a.const_proxy_ ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


	/// @brief Section Proxy Constructor
	inline
	FArray( Section const & s, ProxySentinel const & ) :
		array_size_( s.size() ),
		array_( s.array_ ),
		size_( array_size_ ),
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( true ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


	/// @brief Non-Const Section Proxy Constructor
	inline
	FArray( Section & s, ProxySentinel const & ) :
		array_size_( s.size() ),
		array_( s.array_ ),
		size_( array_size_ ),
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( s.const_proxy_ ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


	/// @brief Value Proxy Constructor
	inline
	FArray( T const & t, ProxySentinel const & ) :
		array_size_( npos ), // Unknown
		array_( const_cast< T * >( &t ) ),
		size_( npos ), // Unbounded
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( true ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


	/// @brief Non-Const Value Proxy Constructor
	inline
	FArray( T & t, ProxySentinel const & ) :
		array_size_( npos ), // Unknown
		array_( &t ),
		size_( npos ), // Unbounded
		owner_( false ),
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_( false ),
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_( 0 ),
		sarray_( 0 )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~FArray()
	{
		if ( owner_ ) delete[] array_;
	}


protected: // Assignment


	/// @brief Copy Assignment
	inline
	FArray &
	operator =( FArray const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		assert( size_ == a.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = a[ i ];
		}
		return *this;
	}


	/// @brief Copy Assignment Template
	template< typename U >
	inline
	void
	operator =( FArray< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		assert( size_ == a.size() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( a[ i ] );
		}
	}


	/// @brief += Array Template
	template< typename U >
	inline
	void
	operator +=( FArray< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		assert( size_ == a.size() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += T( a[ i ] );
		}
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	void
	operator -=( FArray< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		assert( size_ == a.size() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= T( a[ i ] );
		}
	}


public: // Assignment


	/// @brief = Value
	inline
	FArray &
	operator =( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		std::fill_n( array_, size_, t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray &
	operator +=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += t;
		}
		return *this;
	}


	/// @brief -= Value
	inline
	FArray &
	operator -=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= t;
		}
		return *this;
	}


	/// @brief *= Value
	inline
	FArray &
	operator *=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] *= t;
		}
		return *this;
	}


	/// @brief /= Value
	inline
	FArray &
	operator /=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_bounded() );
		assert( t != T( 0 ) );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] /= t;
		}
		return *this;
	}


public: // Subscript


	/// @brief array[ i ] const: Linear Subscript
	inline
	T const &
	operator []( size_type const i ) const
	{
		assert( ( i < size_ ) || ( size_ == npos ) );
		return array_[ i ];
	}


	/// @brief array[ i ]: Linear Subscript
	inline
	T &
	operator []( size_type const i )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( i < size_ ) || ( size_ == npos ) );
		return array_[ i ];
	}


public: // Predicate


	/// @brief Dimensions Initialized?
	virtual
	bool
	dimensions_initialized() const = 0;


	/// @brief Initializer Active?
	virtual
	bool
	initializer_active() const = 0;


	/// @brief Active?
	inline
	bool
	active() const
	{
		return ( array_ != 0 );
	}


	/// @brief Array Size Bounded?
	inline
	bool
	array_size_bounded() const
	{
		return ( array_size_ != npos );
	}


	/// @brief Array Size Unbounded?
	inline
	bool
	array_size_unbounded() const
	{
		return ( array_size_ == npos );
	}


	/// @brief Active Array Size Bounded?
	inline
	bool
	size_bounded() const
	{
		return ( size_ != npos );
	}


	/// @brief Active Array Size Unbounded?
	inline
	bool
	size_unbounded() const
	{
		return ( size_ == npos );
	}


	/// @brief Owner?
	inline
	bool
	owner() const
	{
		return owner_;
	}


	/// @brief Proxy?
	inline
	bool
	proxy() const
	{
		return ( ! owner_ );
	}


	/// @brief All Elements Default Valued?
	inline
	bool
	is_default() const
	{
		for ( size_type i = 0; i < size_; ++i ) {
			if ( array_[ i ] != Traits::initial_value() ) return false;
		}
		return true;
	}


	/// @brief All Elements Zero?
	inline
	bool
	is_zero() const
	{
		for ( size_type i = 0; i < size_; ++i ) {
			if ( array_[ i ] != T( 0 ) ) return false;
		}
		return true;
	}


	/// @brief Uniform Valued?
	inline
	bool
	is_uniform() const
	{
		if ( size_ <= 1 ) return true;
		T const t( array_[ 0 ] );
		for ( size_type i = 1; i < size_; ++i ) {
			if ( array_[ i ] != t ) return false;
		}
		return true;
	}


	/// @brief Uniform Valued with Specified Value?
	inline
	bool
	is_uniform( T const & t ) const
	{
		for ( size_type i = 0; i < size_; ++i ) {
			if ( array_[ i ] != t ) return false;
		}
		return true;
	}


public: // Inspector


	/// @brief Array Size
	inline
	size_type
	array_size() const
	{
		return array_size_;
	}


	/// @brief Active Array Size
	inline
	size_type
	size() const
	{
		return size_;
	}


public: // Modifier


	/// @brief Clear
	inline
	virtual
	FArray &
	clear()
	{
		array_size_ = 0;
		if ( owner_ ) delete[] array_;
    array_ = 0;
		size_ = 0;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = false;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = 0;
		sarray_ = array_;
		return *this;
	}


	/// @brief Assign Default Value to all Elements
	inline
	virtual
	FArray &
	to_default()
	{
		proxy_const_assert( not_const_proxy() );
		std::fill_n( array_, size_, Traits::initial_value() );
		return *this;
	}


	/// @brief Assign Zero to all Elements
	/// @note  Can't be virtual (for covariant return) or will try to instantiate for all value types
	inline
	void
	zero()
	{
		proxy_const_assert( not_const_proxy() );
		std::fill_n( array_, size_, T( 0 ) );
	}


	/// @brief Assign Zero to all Elements
	/// @note  Can't be virtual (for covariant return) or will try to instantiate for all value types
	inline
	void
	to_zero()
	{
		proxy_const_assert( not_const_proxy() );
		std::fill_n( array_, size_, T( 0 ) );
	}


protected: // Comparison


	/// @brief FArray == FArray
	friend
	inline
	bool
	operator ==( FArray const & a, FArray const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ == b.size_ ) { // Sizes match
			for ( typename FArray::size_type i = 0, e = a.size_; i < e; ++i ) {
				if ( a[ i ] != b[ i ] ) return false;
			}
			return true;
		} else { // Sizes differ
			return false;
		}
	}


	/// @brief FArray != FArray
	friend
	inline
	bool
	operator !=( FArray const & a, FArray const & b )
	{
		return !( a == b );
	}


protected: // Functions


	/// @brief Array Size Product of Specified Bounded Dimensional Sizes
	inline
	static
	size_type
	size_of( size_type const s1 )
	{
		assert( s1 != npos );
		return s1;
	}


	/// @brief Array Size Product of Specified Bounded Dimensional Sizes
	inline
	static
	size_type
	size_of( size_type const s1, size_type const s2 )
	{
		assert( s1 != npos );
		assert( s2 != npos );
		assert( ( s2 == 0 ) || ( s1 <= max_size / s2 ) );
		return s1 * s2;
	}


	/// @brief Array Size Product of Specified Bounded Dimensional Sizes
	inline
	static
	size_type
	size_of( size_type const s1, size_type const s2, size_type const s3 )
	{
		return size_of( size_of( s1, s2 ), s3 );
	}


	/// @brief Array Size Product of Specified Bounded Dimensional Sizes
	inline
	static
	size_type
	size_of( size_type const s1, size_type const s2, size_type const s3, size_type const s4 )
	{
		return size_of( size_of( s1, s2 ), size_of( s3, s4 ) );
	}


	/// @brief Array Size Product of Specified Bounded Dimensional Sizes
	inline
	static
	size_type
	size_of( size_type const s1, size_type const s2, size_type const s3, size_type const s4, size_type const s5 )
	{
		return size_of( size_of( size_of( s1, s2 ), size_of( s3, s4 ) ), s5 );
	}


	/// @brief Array Size Product of Specified Bounded Dimensional Sizes
	inline
	static
	size_type
	size_of( size_type const s1, size_type const s2, size_type const s3, size_type const s4, size_type const s5, size_type const s6 )
	{
		return size_of( size_of( size_of( s1, s2 ), size_of( s3, s4 ) ), size_of( s5, s6 ) );
	}


	/// @brief Array Size Product of Specified Bounded IndexRanges
	inline
	static
	size_type
	size_of( IR const & I1 )
	{
		return size_of( I1.size() );
	}


	/// @brief Array Size Product of Specified Bounded IndexRanges
	inline
	static
	size_type
	size_of( IR const & I1, IR const & I2 )
	{
		return size_of( I1.size(), I2.size() );
	}


	/// @brief Array Size Product of Specified Bounded IndexRanges
	inline
	static
	size_type
	size_of( IR const & I1, IR const & I2, IR const & I3 )
	{
		return size_of( I1.size(), I2.size(), I3.size() );
	}


	/// @brief Array Size Product of Specified Bounded IndexRanges
	inline
	static
	size_type
	size_of( IR const & I1, IR const & I2, IR const & I3, IR const & I4 )
	{
		return size_of( I1.size(), I2.size(), I3.size(), I4.size() );
	}


	/// @brief Array Size Product of Specified Bounded IndexRanges
	inline
	static
	size_type
	size_of( IR const & I1, IR const & I2, IR const & I3, IR const & I4, IR const & I5 )
	{
		return size_of( I1.size(), I2.size(), I3.size(), I4.size(), I5.size() );
	}


	/// @brief Array Size Product of Specified Bounded IndexRanges
	inline
	static
	size_type
	size_of( IR const & I1, IR const & I2, IR const & I3, IR const & I4, IR const & I5, IR const & I6 )
	{
		return size_of( I1.size(), I2.size(), I3.size(), I4.size(), I5.size(), I6.size() );
	}


	/// @brief Shift Setup
	inline
	void
	shift_set( int const shift_a )
	{
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Active Array Size Setup
	inline
	void
	size_set( size_type const size_a )
	{
		assert( size_a <= array_size_ );
		size_ = size_a;
	}


	/// @brief Resize a Real Array
	inline
	FArray &
	resize( size_type const size_a )
	{
		assert( owner_ );
		assert( size_a != npos );
		if ( array_size_ != size_a ) {
			array_size_ = size_a;
			delete[] array_; array_ = ( array_size_ > 0 ? new T[ array_size_ ] : 0 );
			size_ = size_a;
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
			size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
		}
#ifdef OBJEXXFCL_FARRAY_INIT
		if ( ! initializer_active() ) {
			std::fill_n( array_, size_, Traits::initial_value() );
		}
#endif // OBJEXXFCL_FARRAY_INIT
		return *this;
	}


	/// @brief Attach Proxy/Argument Array to Const Array of Same Rank
	inline
	void
	attach( FArray const & a )
	{
		assert( ! owner_ );
		array_size_ = a.array_size_;
		array_ = a.array_;
		size_ = a.size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = a.shift_;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Array of Same Rank
	inline
	void
	attach( FArray & a )
	{
		assert( ! owner_ );
		array_size_ = a.array_size_;
		array_ = a.array_;
		size_ = a.size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = a.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = a.shift_;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Const Array
	inline
	void
	attach( FArray const & a, int const shift_a )
	{
		assert( ! owner_ );
		array_size_ = a.array_size_;
		array_ = a.array_;
		size_ = a.size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Array
	inline
	void
	attach( FArray & a, int const shift_a )
	{
		assert( ! owner_ );
		array_size_ = a.array_size_;
		array_ = a.array_;
		size_ = a.size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = a.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Const Section
	inline
	void
	attach( Section const & s, int const shift_a )
	{
		assert( ! owner_ );
		array_size_ = s.size();
		array_ = s.array_;
		size_ = array_size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Section
	inline
	void
	attach( Section & s, int const shift_a )
	{
		assert( ! owner_ );
		array_size_ = s.size();
		array_ = s.array_;
		size_ = array_size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = s.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Const Value
	inline
	void
	attach( T const & t, int const shift_a )
	{
		assert( ! owner_ );
		array_size_ = npos; // Unknown
		array_ = const_cast< T * >( &t );
		size_ = npos; // Unbounded
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Attach Proxy/Argument Array to Value
	inline
	void
	attach( T & t, int const shift_a )
	{
		assert( ! owner_ );
		array_size_ = npos; // Unknown
		array_ = &t;
		size_ = npos; // Unbounded
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = false;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = shift_a;
		sarray_ = array_ - shift_;
	}


	/// @brief Detach Proxy/Argument Array
	inline
	void
	detach()
	{
		assert( ! owner_ );
		array_size_ = 0;
		array_ = 0;
		size_ = 0;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = false;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		shift_ = 0;
		sarray_ = 0;
	}


	/// @brief Update Proxy Array Attachment to Const Array
	inline
	void
	update_to( FArray const & a )
	{
		assert( ! owner_ );
		array_size_ = a.array_size_;
		array_ = a.array_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	}


	/// @brief Update Proxy Array Attachment to Array
	inline
	void
	update_to( FArray & a )
	{
		assert( ! owner_ );
		array_size_ = a.array_size_;
		array_ = a.array_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = a.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	}


	/// @brief Swap
	inline
	void
	swapB( FArray & v )
	{
		assert( owner_ );
		assert( v.owner_ );
		std::swap( array_size_, v.array_size_ );
		std::swap( array_, v.array_ );
		std::swap( size_, v.size_ );
		std::swap( shift_, v.shift_ );
		std::swap( sarray_, v.sarray_ );
	}


#ifdef OBJEXXFCL_PROXY_CONST_CHECKS

	/// @brief Const Proxy?
	inline
	bool
	const_proxy() const
	{
		return const_proxy_;
	}

#endif // OBJEXXFCL_PROXY_CONST_CHECKS


#ifdef OBJEXXFCL_PROXY_CONST_CHECKS

	/// @brief Not a Const Proxy Under Strict Const-Correctness?
	inline
	bool
	not_const_proxy() const
	{
		return ( ! const_proxy_ );
	}

#endif // OBJEXXFCL_PROXY_CONST_CHECKS


private: // Properties


#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT

	/// @brief Report size if at least value defined for OBJEXXFCL_FARRAY_SIZE_REPORT
	/// @note  Size is based on sizeof( T ) so T-controlled heap memory is not counted
	inline
	void
	size_report() const
	{
		if ( size_ * sizeof( T ) >= OBJEXXFCL_FARRAY_SIZE_REPORT ) {
			std::cout << "FArray< " << typeid( T ).name() << " >"
			 << "  Size: " << size_ * sizeof( T ) << "  Elements: " << size_;
		}
	}

#endif // OBJEXXFCL_FARRAY_SIZE_REPORT


public: // Data


	/// @brief Unbounded "size"
	static size_type const npos = static_cast< size_type >( -1 );

	/// @brief Max array size
	static size_type const max_size = npos - static_cast< size_type >( 1 );


protected: // Data


	/// @brief Size of data array
	size_type array_size_;

	/// @brief Pointer to data array
	T * array_;

	/// @brief Size of active array
	size_type size_;

	/// @brief Owner of data array?
	bool const owner_;

#ifdef OBJEXXFCL_PROXY_CONST_CHECKS

	/// @brief Proxy for const data array?
	bool const_proxy_;

#endif // OBJEXXFCL_PROXY_CONST_CHECKS

	/// @brief Array shift
	int shift_;

	/// @brief Shifted pointer to data array
	T * sarray_;


}; // FArray


// Static Data Member Template Definitions

template< typename T > typename FArray< T >::size_type const FArray< T >::npos;

template< typename T > typename FArray< T >::size_type const FArray< T >::max_size;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray_HH
