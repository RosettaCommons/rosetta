#ifndef INCLUDED_ObjexxFCL_CArrayP_hh
#define INCLUDED_ObjexxFCL_CArrayP_hh


// CArrayP: Memory-Managed C Array Wrapper Supporting Proxies
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
#include <ObjexxFCL/CArrayP.fwd.hh>
#include <ObjexxFCL/proxy_const_assert.hh>

// C++ Headers
#include <algorithm>
#include <cmath>


namespace ObjexxFCL {


/// @brief CArrayP: Memory-Managed C Array Wrapper Supporting Proxies
///
/// @note  Proxy CArrayPs are invalidated if the underlying (owning) array data is deleted
/// @note  Proxy CPArras can be created at construction with the Proxy named constructors
/// @note  CArrayPs can become proxies with the attach() member function
/// @note  CArrayPs can stop being proxies with the detach() member function
template< typename T >
class CArrayP
{


private: // Friend


	template< typename > friend class CArrayP; // Friendship across value types


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
	CArrayP() :
		size_( 0 ),
		array_( 0 ),
		owner_( true )
	{}


	/// @brief Copy constructor
	inline
	CArrayP( CArrayP const & a ) :
		size_( a.size_ ),
		array_( a.owner_ ? ( size_ > 0 ? new T[ size_ ] : 0 ) : a.array_ ),
		owner_( a.owner_ )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( ! a.owner_ )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{
		if ( owner_ ) {
			for ( size_type i = 0; i < size_; ++i ) {
				array_[ i ] = a.array_[ i ];
			}
		}
	}


	/// @brief Non-Const Copy constructor
	inline
	CArrayP( CArrayP & a ) :
		size_( a.size_ ),
		array_( a.owner_ ? ( size_ > 0 ? new T[ size_ ] : 0 ) : a.array_ ),
		owner_( a.owner_ )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( a.const_proxy_ )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{
		if ( owner_ ) {
			for ( size_type i = 0; i < size_; ++i ) {
				array_[ i ] = a.array_[ i ];
			}
		}
	}


	/// @brief Copy constructor template
	template< typename U >
	inline
	CArrayP( CArrayP< U > const & a ) :
		size_( a.size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 ),
		owner_( true )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( a.array_[ i ] );
		}
	}


	/// @brief Pointer + size constructor
	inline
	CArrayP(
		T const * const p,
		size_type const size_a
	) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 ),
		owner_( true )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = p[ i ];
		}
	}


	/// @brief Pointer + size constructor template
	template< typename U >
	inline
	CArrayP(
		U const * const p,
		size_type const size_a
	) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 ),
		owner_( true )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( p[ i ] );
		}
	}


	/// @brief Iterator range constructor template
	template< typename InputIterator >
	inline
	CArrayP(
		InputIterator const beg,
		InputIterator const end
	) :
		size_( end - beg ),
		array_( size_ > 0 ? new T[ size_ ] : 0 ),
		owner_( true )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
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
	CArrayP( size_type const size_a ) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 ),
		owner_( true )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Size + uniform value constructor
	inline
	CArrayP(
		size_type const size_a,
		T const & t
	) :
		size_( size_a ),
		array_( size_ > 0 ? new T[ size_ ] : 0 ),
		owner_( true )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = t;
		}
	}


	/// @brief Proxy const copy named constructor
	static
	inline
	CArrayP
	Proxy( CArrayP const & a )
	{
		CArrayP p;
		p.size_ = a.size_;
		p.array_ = a.array_;
		p.owner_ = false;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		p.const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return p;
	}


	/// @brief Proxy copy named constructor
	static
	inline
	CArrayP
	Proxy( CArrayP & a )
	{
		CArrayP p;
		p.size_ = a.size_;
		p.array_ = a.array_;
		p.owner_ = false;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		p.const_proxy_ = a.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return p;
	}


	/// @brief Proxy const copy + size named constructor
	static
	inline
	CArrayP
	Proxy(
		CArrayP const & a,
		size_type const size_a
	)
	{
		assert( size_a <= a.size_ );
		CArrayP p;
		p.size_ = size_a;
		p.array_ = a.array_;
		p.owner_ = false;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		p.const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return p;
	}


	/// @brief Proxy copy + size named constructor
	static
	inline
	CArrayP
	Proxy(
		CArrayP & a,
		size_type const size_a
	)
	{
		assert( size_a <= a.size_ );
		CArrayP p;
		p.size_ = size_a;
		p.array_ = a.array_;
		p.owner_ = false;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		p.const_proxy_ = a.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return p;
	}


	/// @brief Destructor
	inline
	~CArrayP()
	{
		if ( owner_ ) delete[] array_;
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
	CArrayP &
	operator =( CArrayP const & a )
	{
		proxy_const_assert( not_const_proxy() );
		if ( this != &a ) {
			if ( size_ != a.size_ ) {
				assert( owner_ );
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
	CArrayP &
	operator =( CArrayP< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		if ( size_ != a.size_ ) {
			assert( owner_ );
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
	CArrayP &
	operator =( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = t;
		}
		return *this;
	}


	/// @brief Pointer + size assignment
	inline
	CArrayP &
	assign(
		T const * const p,
		size_type const size_a
	)
	{
		proxy_const_assert( not_const_proxy() );
		if ( size_ != size_a ) {
			assert( owner_ );
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
	CArrayP &
	assign(
		U const * const p,
		size_type const size_a
	)
	{
		proxy_const_assert( not_const_proxy() );
		if ( size_ != size_a ) {
			assert( owner_ );
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
	CArrayP &
	assign(
		InputIterator const beg,
		InputIterator const end
	)
	{
		proxy_const_assert( not_const_proxy() );
		size_type const size_a( end - beg );
		if ( size_ != size_a ) {
			assert( owner_ );
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
	CArrayP &
	assign(
		size_type const size_a,
		T const & value
	)
	{
		proxy_const_assert( not_const_proxy() );
		if ( size_ != size_a ) { // Set to new array with uniform values
			assert( owner_ );
			CArrayP( size_a, value ).swap( *this );
		} else { // Set to uniform value
			(*this) = value;
		}
		return *this;
	}


	/// @brief += CArrayP
	template< typename U >
	inline
	CArrayP &
	operator +=( CArrayP< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_ == a.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += T( a.array_[ i ] );
		}
		return *this;
	}


	/// @brief -= CArrayP
	template< typename U >
	inline
	CArrayP &
	operator -=( CArrayP< U > const & a )
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_ == a.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= T( a.array_[ i ] );
		}
		return *this;
	}


	/// @brief += value
	inline
	CArrayP &
	operator +=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += t;
		}
		return *this;
	}


	/// @brief -= value
	inline
	CArrayP &
	operator -=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= t;
		}
		return *this;
	}


	/// @brief *= value
	inline
	CArrayP &
	operator *=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] *= t;
		}
		return *this;
	}


	/// @brief /= value
	inline
	CArrayP &
	operator /=( T const & t )
	{
		proxy_const_assert( not_const_proxy() );
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


public: // Inspector


	/// @brief Size
	inline
	size_type
	size() const
	{
		return size_;
	}


	/// @brief Lower index
	inline
	size_type
	l() const
	{
		return 0u;
	}


	/// @brief Upper index
	inline
	size_type
	u() const
	{
		return size_ - 1; // npos if size_ == 0
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
		proxy_const_assert( not_const_proxy() );
		assert( size_ > 0 );
		return array_[ 0 ];
	}


	/// @brief Last element
	inline
	T &
	back()
	{
		proxy_const_assert( not_const_proxy() );
		assert( size_ > 0 );
		return array_[ size_ - 1 ];
	}


	/// @brief Resize: Values not preserved
	/// @note Built-in values are uninitialized if size changes
	inline
	CArrayP &
	size( size_type const size_a )
	{
		assert( owner_ );
		if ( size_ != size_a ) { // Set to new array
			CArrayP( size_a ).swap( *this );
		}
		return *this;
	}


	/// @brief Resize to size with fill value: Values preserved
	inline
	CArrayP &
	resize(
		size_type const size_a,
		T const & fill = T()
	)
	{
		assert( owner_ );
		if ( size_ < size_a ) {
			CArrayP a( size_a, fill ); // New array: Elements set to fill fill
			for ( size_type i = 0; i < size_; ++i ) { // Copy current values
				a.array_[ i ] = array_[ i ];
			}
			swap( a ); // Swap in new array
		} else if ( size_ > size_a ) {
			CArrayP a( size_a ); // New array
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
	swap( CArrayP & a )
	{
		std::swap( size_, a.size_ );
		std::swap( array_, a.array_ );
		std::swap( owner_, a.owner_ );
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		std::swap( const_proxy_, a.const_proxy_ );
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	}


	/// @brief Clear
	inline
	CArrayP &
	clear()
	{
		size_ = 0;
		if ( owner_ ) delete[] array_; array_ = 0;
		owner_ = true;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = false;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return *this;
	}


	/// @brief Normalize to unit length
	inline
	CArrayP &
	normalize()
	{
		T const length_( length() );
		assert( length_ > T( 0 ) );
		operator /=( length_ );
		return *this;
	}


	/// @brief Attach as proxy to a const CArrayP
	inline
	CArrayP &
	attach( CArrayP const & a )
	{
		size_ = a.size_;
		if ( owner_ ) delete[] array_; array_ = a.array_;
		owner_ = false;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = true;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return *this;
	}


	/// @brief Attach as proxy to a CArrayP
	inline
	CArrayP &
	attach( CArrayP & a )
	{
		size_ = a.size_;
		if ( owner_ ) delete[] array_; array_ = a.array_;
		owner_ = false;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		const_proxy_ = a.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		return *this;
	}


	/// @brief Detach as proxy to a CArrayP
	inline
	CArrayP &
	detach()
	{
		if ( ! owner_ ) clear(); // Proxy: Clear fields
		return *this;
	}


public: // Subscript


	/// @brief CArrayP[ i ] const: 0-based indexing
	inline
	T const &
	operator []( size_type const i ) const
	{
		assert( i < size_ );
		return array_[ i ];
	}


	/// @brief CArrayP[ i ]: 0-based indexing
	inline
	T &
	operator []( size_type const i )
	{
		proxy_const_assert( not_const_proxy() );
		assert( i < size_ );
		return array_[ i ];
	}


#ifdef OBJEXXFCL_CARRAY_1_BASED_LOOKUP


	/// @brief CArrayP( i ) const: 1-based indexing
	inline
	T const &
	operator ()( size_type const i ) const
	{
		assert( ( i > 0 ) && ( i <= size_ ) );
		return array_[ i - 1 ];
	}


	/// @brief CArrayP( i ): 1-based indexing
	inline
	T &
	operator ()( size_type const i )
	{
		proxy_const_assert( not_const_proxy() );
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
		proxy_const_assert( not_const_proxy() );
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
		proxy_const_assert( not_const_proxy() );
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
		proxy_const_assert( not_const_proxy() );
		return array_;
	}


public: // Comparison


	/// @brief Are two CArrayPs comparable?
	friend
	inline
	bool
	comparable( CArrayP const & a, CArrayP const & b )
	{
		return ( a.size_ == b.size_ );
	}


	/// @brief CArrayP == CArrayP
	friend
	inline
	bool
	operator ==( CArrayP const & a, CArrayP const & b )
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


	/// @brief CArrayP != CArrayP
	friend
	inline
	bool
	operator !=( CArrayP const & a, CArrayP const & b )
	{
		return !( a == b );
	}


	/// @brief CArrayP < CArrayP
	friend
	inline
	bool
	operator <( CArrayP const & a, CArrayP const & b )
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


	/// @brief CArrayP <= CArrayP
	friend
	inline
	bool
	operator <=( CArrayP const & a, CArrayP const & b )
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


	/// @brief CArrayP >= CArrayP
	friend
	inline
	bool
	operator >=( CArrayP const & a, CArrayP const & b )
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


	/// @brief CArrayP > CArrayP
	friend
	inline
	bool
	operator >( CArrayP const & a, CArrayP const & b )
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


	/// @brief CArrayP == T
	friend
	inline
	bool
	operator ==( CArrayP const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( a.array_[ i ] != t ) return false;
		}
		return true;
	}


	/// @brief CArrayP != T
	friend
	inline
	bool
	operator !=( CArrayP const & a, T const & t )
	{
		return !( a == t );
	}


	/// @brief CArrayP < T
	friend
	inline
	bool
	operator <( CArrayP const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] < t ) ) return false;
		}
		return true;
	}


	/// @brief CArrayP <= T
	friend
	inline
	bool
	operator <=( CArrayP const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] <= t ) ) return false;
		}
		return true;
	}


	/// @brief CArrayP >= T
	friend
	inline
	bool
	operator >=( CArrayP const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] >= t ) ) return false;
		}
		return true;
	}


	/// @brief CArrayP > T
	friend
	inline
	bool
	operator >( CArrayP const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a.array_[ i ] > t ) ) return false;
		}
		return true;
	}


	/// @brief T == CArrayP
	friend
	inline
	bool
	operator ==( T const & t, CArrayP const & a )
	{
		return ( a == t );
	}


	/// @brief T != CArrayP
	friend
	inline
	bool
	operator !=( T const & t, CArrayP const & a )
	{
		return !( t == a );
	}


	/// @brief T < CArrayP
	friend
	inline
	bool
	operator <( T const & t, CArrayP const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t < a.array_[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T <= CArrayP
	friend
	inline
	bool
	operator <=( T const & t, CArrayP const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t <= a.array_[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T >= CArrayP
	friend
	inline
	bool
	operator >=( T const & t, CArrayP const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t >= a.array_[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T > CArrayP
	friend
	inline
	bool
	operator >( T const & t, CArrayP const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t > a.array_[ i ] ) ) return false;
		}
		return true;
	}


public: // Generator


	/// @brief -CArrayP
	friend
	inline
	CArrayP
	operator -( CArrayP const & a )
	{
		CArrayP r( a );
		r *= T( -1 );
		return r;
	}


	/// @brief CArrayP + CArrayP
	friend
	inline
	CArrayP
	operator +( CArrayP const & a, CArrayP const & b )
	{
		CArrayP r( a );
		r += b;
		return r;
	}


	/// @brief CArrayP - CArrayP
	friend
	inline
	CArrayP
	operator -( CArrayP const & a, CArrayP const & b )
	{
		CArrayP r( a );
		r -= b;
		return r;
	}


	/// @brief CArrayP + Value
	friend
	inline
	CArrayP
	operator +( CArrayP const & a, T const & t )
	{
		CArrayP r( a );
		r += t;
		return r;
	}


	/// @brief Value + CArrayP
	friend
	inline
	CArrayP
	operator +( T const & t, CArrayP const & a )
	{
		CArrayP r( a );
		r += t;
		return r;
	}


	/// @brief CArrayP - Value
	friend
	inline
	CArrayP
	operator -( CArrayP const & a, T const & t )
	{
		CArrayP r( a );
		r -= t;
		return r;
	}


	/// @brief Value - CArrayP
	friend
	inline
	CArrayP
	operator -( T const & t, CArrayP const & a )
	{
		CArrayP r( -a );
		r += t;
		return r;
	}


	/// @brief CArrayP * Value
	friend
	inline
	CArrayP
	operator *( CArrayP const & a, T const & t )
	{
		CArrayP r( a );
		r *= t;
		return r;
	}


	/// @brief Value * CArrayP
	friend
	inline
	CArrayP
	operator *( T const & t, CArrayP const & a )
	{
		CArrayP r( a );
		r *= t;
		return r;
	}


	/// @brief CArrayP / Value
	friend
	inline
	CArrayP
	operator /( CArrayP const & a, T const & t )
	{
		CArrayP r( a );
		r /= t;
		return r;
	}


public: // Friend


	/// @brief Dot product
	friend
	inline
	T
	dot_product( CArrayP const & a, CArrayP const & b )
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
	dot( CArrayP const & a, CArrayP const & b )
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
	distance( CArrayP const & a, CArrayP const & b )
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
	distance_squared( CArrayP const & a, CArrayP const & b )
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
	swap( CArrayP & a, CArrayP & b )
	{
		a.swap( b );
	}


private: // Functions


#ifdef OBJEXXFCL_PROXY_CONST_CHECKS

	/// @brief Not a Const Proxy Under Strict Const-Correctness?
	inline
	bool
	not_const_proxy() const
	{
		return ( ! const_proxy_ );
	}

#endif // OBJEXXFCL_PROXY_CONST_CHECKS


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

	/// @brief Owner of the data array or proxy?
	bool owner_;

#ifdef OBJEXXFCL_PROXY_CONST_CHECKS

	/// @brief Proxy for const data array?
	bool const_proxy_;

#endif // OBJEXXFCL_PROXY_CONST_CHECKS


}; // CArrayP


/// @brief Are two CArrayPs comparable?
template< typename T >
bool
comparable( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP == CArrayP
template< typename T >
bool
operator ==( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP != CArrayP
template< typename T >
bool
operator !=( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP < CArrayP
template< typename T >
bool
operator <( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP <= CArrayP
template< typename T >
bool
operator <=( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP >= CArrayP
template< typename T >
bool
operator >=( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP > CArrayP
template< typename T >
bool
operator >( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP == T
template< typename T >
bool
operator ==( CArrayP< T > const & a, T const & t );


/// @brief CArrayP != T
template< typename T >
bool
operator !=( CArrayP< T > const & a, T const & t );


/// @brief CArrayP < T
template< typename T >
bool
operator <( CArrayP< T > const & a, T const & t );


/// @brief CArrayP <= T
template< typename T >
bool
operator <=( CArrayP< T > const & a, T const & t );


/// @brief CArrayP >= T
template< typename T >
bool
operator >=( CArrayP< T > const & a, T const & t );


/// @brief CArrayP > T
template< typename T >
bool
operator >( CArrayP< T > const & a, T const & t );


/// @brief T == CArrayP
template< typename T >
bool
operator ==( T const & t, CArrayP< T > const & a );


/// @brief T != CArrayP
template< typename T >
bool
operator !=( T const & t, CArrayP< T > const & a );


/// @brief T < CArrayP
template< typename T >
bool
operator <( T const & t, CArrayP< T > const & a );


/// @brief T <= CArrayP
template< typename T >
bool
operator <=( T const & t, CArrayP< T > const & a );


/// @brief T >= CArrayP
template< typename T >
bool
operator >=( T const & t, CArrayP< T > const & a );


/// @brief T > CArrayP
template< typename T >
bool
operator >( T const & t, CArrayP< T > const & a );


/// @brief -CArrayP
template< typename T >
CArrayP< T >
operator -( CArrayP< T > const & a );


/// @brief CArrayP + CArrayP
template< typename T >
CArrayP< T >
operator +( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP - CArrayP
template< typename T >
CArrayP< T >
operator -( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief CArrayP + Value
template< typename T >
CArrayP< T >
operator +( CArrayP< T > const & a, T const & t );


/// @brief Value + CArrayP
template< typename T >
CArrayP< T >
operator +( T const & t, CArrayP< T > const & a );


/// @brief CArrayP - Value
template< typename T >
CArrayP< T >
operator -( CArrayP< T > const & a, T const & t );


/// @brief Value - CArrayP
template< typename T >
CArrayP< T >
operator -( T const & t, CArrayP< T > const & a );


/// @brief CArrayP * Value
template< typename T >
CArrayP< T >
operator *( CArrayP< T > const & a, T const & t );


/// @brief Value * CArrayP
template< typename T >
CArrayP< T >
operator *( T const & t, CArrayP< T > const & a );


/// @brief CArrayP / Value
template< typename T >
CArrayP< T >
operator /( CArrayP< T > const & a, T const & t );


/// @brief Dot product
template< typename T >
T
dot_product( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief Dot product
template< typename T >
T
dot( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief Distance
template< typename T >
T
distance( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief Distance squared
template< typename T >
T
distance_squared( CArrayP< T > const & a, CArrayP< T > const & b );


/// @brief Swap
template< typename T >
void
swap( CArrayP< T > & a, CArrayP< T > & b );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.  The legal alternative would be
// to add specializations of swap for each anticipated instantiation.


namespace std {


/// @brief std::swap( CArrayP, CArrayP )
template< typename T >
inline
void
swap( ObjexxFCL::CArrayP< T > & a, ObjexxFCL::CArrayP< T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_CArrayP_HH
