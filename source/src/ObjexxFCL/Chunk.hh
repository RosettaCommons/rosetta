#ifndef INCLUDED_ObjexxFCL_Chunk_hh
#define INCLUDED_ObjexxFCL_Chunk_hh


// Chunk: Contiguous Array for Use in ChunkVector
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


// C++ Headers
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>


namespace ObjexxFCL {


/// @brief Chunk: Contiguous Array for Use in ChunkVector
///
/// @remarks
///  @li size <= capacity
///  @li capacity == size after construction
///  @li capacity == size after assignment if reallocation required
template< typename T >
class Chunk
{


private: // Friend


	template< typename > friend class Chunk;


public: // Types


	// STL style
	typedef  T  value_type;
	typedef  T &  reference;
	typedef  T const &  const_reference;
	typedef  T *  pointer;
	typedef  T const *  const_pointer;
	typedef  std::size_t  size_type;

	// C++ style
	typedef  T  Value;
	typedef  T &  Reference;
	typedef  T const &  ConstReference;
	typedef  T *  Pointer;
	typedef  T const *  ConstPointer;
	typedef  std::size_t  Size;


public: // Creation


	/// @brief Default Constructor
	inline
	Chunk() :
		size_( 0 ),
		capacity_( 0 ),
		array_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	Chunk( Chunk const & c ) :
		size_( c.size_ ),
		capacity_( size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = c.array_[ i ];
		}
	}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	Chunk( Chunk< U > const & c ) :
		size_( c.size_ ),
		capacity_( size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( c.array_[ i ] );
		}
	}


	/// @brief Size Constructor: Built-In Types are Not Initialized!
	inline
	explicit
	Chunk( size_type const size_a ) :
		size_( size_a ),
		capacity_( size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{}


	/// @brief Size + Uniform Value Constructor
	inline
	Chunk(
		size_type const size_a,
		T const & value
	) :
		size_( size_a ),
		capacity_( size_ ),
		array_( size_ > 0 ? new T[ size_ ] : 0 )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = value;
		}
	}


	/// @brief Destructor
	inline
	~Chunk()
	{
		delete[] array_;
	}


public: // Assignment


	/// @brief Copy Assignment
	inline
	Chunk &
	operator =( Chunk const & c )
	{
		if ( this != &c ) {
			if ( size_ != c.size_ ) {
				size_ = c.size_;
				capacity_ = size_;
				delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
			}
			for ( size_type i = 0; i < size_; ++i ) {
				array_[ i ] = c.array_[ i ];
			}
		}
		return *this;
	}


	/// @brief Copy Assignment Template
	template< typename U >
	inline
	Chunk &
	operator =( Chunk< U > const & c )
	{
		if ( size_ != c.size_ ) {
			size_ = c.size_;
			capacity_ = size_;
			delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
		}
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = T( c.array_[ i ] );
		}
		return *this;
	}


	/// @brief Size + Value Assignment
	inline
	Chunk &
	assign(
		size_type const size_a,
		T const & value
	)
	{
		if ( size_ != size_a ) {
			size_ = size_a;
			capacity_ = size_;
			delete[] array_; array_ = ( size_ > 0 ? new T[ size_ ] : 0 );
		}
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = value;
		}
		return *this;
	}


	/// @brief += Chunk
	inline
	Chunk &
	operator +=( Chunk const & c )
	{
		assert( size_ == c.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += c.array_[ i ];
		}
		return *this;
	}


	/// @brief -= Chunk
	inline
	Chunk &
	operator -=( Chunk const & c )
	{
		assert( size_ == c.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= c.array_[ i ];
		}
		return *this;
	}


	/// @brief += Chunk Template
	template< typename U >
	inline
	Chunk &
	operator +=( Chunk< U > const & c )
	{
		assert( size_ == c.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += T( c.array_[ i ] );
		}
		return *this;
	}


	/// @brief -= Chunk Template
	template< typename U >
	inline
	Chunk &
	operator -=( Chunk< U > const & c )
	{
		assert( size_ == c.size_ );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= T( c.array_[ i ] );
		}
		return *this;
	}


	/// @brief = Value
	inline
	Chunk &
	operator =( T const & value )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = value;
		}
		return *this;
	}


	/// @brief += Value
	inline
	Chunk &
	operator +=( T const & value )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] += value;
		}
		return *this;
	}


	/// @brief -= Value
	inline
	Chunk &
	operator -=( T const & value )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] -= value;
		}
		return *this;
	}


	/// @brief *= Value
	inline
	Chunk &
	operator *=( T const & value )
	{
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] *= value;
		}
		return *this;
	}


	/// @brief /= Value
	inline
	Chunk &
	operator /=( T const & value )
	{
		assert( value != T( 0 ) );
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] /= value;
		}
		return *this;
	}


public: // Subscript


	/// @brief Chunk[ i ] const: 0-Based Indexing
	inline
	T const &
	operator []( size_type const i ) const
	{
		assert( i < size_ );
		return array_[ i ];
	}


	/// @brief Chunk[ i ]: 0-Based Indexing
	inline
	T &
	operator []( size_type const i )
	{
		assert( i < size_ );
		return array_[ i ];
	}


public: // Inspector


	/// @brief Size
	inline
	size_type
	size() const
	{
		return size_;
	}


	/// @brief Capacity
	inline
	size_type
	capacity() const
	{
		return capacity_;
	}


	/// @brief Maximum Size
	inline
	size_type
	max_size() const
	{
		return std::numeric_limits< size_type >::max();
	}


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return ( size_ == 0 );
	}


	/// @brief First Element
	inline
	T const &
	front() const
	{
		assert( size_ > 0 );
		return array_[ 0 ];
	}


	/// @brief Last Element
	inline
	T const &
	back() const
	{
		assert( size_ > 0 );
		return array_[ size_ - 1 ];
	}


public: // Modifier


	/// @brief First Element
	inline
	T &
	front()
	{
		assert( size_ > 0 );
		return array_[ 0 ];
	}


	/// @brief Last Element
	inline
	T &
	back()
	{
		assert( size_ > 0 );
		return array_[ size_ - 1 ];
	}


	/// @brief Append an Element
	inline
	Chunk &
	push_back( T const & value )
	{
		assert( size_ < max_size() );
		if ( size_ == capacity_ ) reserve( 2 * capacity_ );
		array_[ size_ ] = value;
		++size_;
		return *this;
	}


	/// @brief Remove the Last Element
	inline
	Chunk &
	pop_back()
	{
		assert( size_ > 0 );
		--size_;
		return *this;
	}


	/// @brief Resize: Values Preserved: Added Built-In Values are Not Initialized!
	inline
	Chunk &
	resize( size_type const size_a )
	{
		if ( size_ != size_a ) {
			if ( size_a > capacity_ ) {
				size_type const size_c( std::min( size_, size_a ) );
				T * const new_array( new T[ size_a ] );
				for ( size_type i = 0; i < size_c; ++i ) {
					new_array[ i ] = array_[ i ];
				}
				delete[] array_; array_ = new_array;
				capacity_ = size_a;
			}
			size_ = size_a;
		}
		return *this;
	}


	/// @brief Resize + Fill Value: Values Preserved
	inline
	Chunk &
	resize(
		size_type const size_a,
		T const & value
	)
	{
		if ( size_ != size_a ) {
			size_type const size_c( std::min( size_, size_a ) );
			if ( size_a > capacity_ ) {
				T * const new_array( new T[ size_a ] );
				for ( size_type i = 0; i < size_c; ++i ) {
					new_array[ i ] = array_[ i ];
				}
				delete[] array_; array_ = new_array;
				capacity_ = size_a;
			}
			for ( size_type i = size_c; i < size_a; ++i ) {
				array_[ i ] = value;
			}
			size_ = size_a;
		}
		return *this;
	}


	/// @brief Resize: Values Not Preserved: Built-In Values are Not Initialized!
	inline
	Chunk &
	non_preserving_resize( size_type const size_a )
	{
		if ( size_ != size_a ) {
			if ( size_a > capacity_ ) {
				delete[] array_; array_ = new T[ size_a ];
				capacity_ = size_a;
			}
			size_ = size_a;
		}
		return *this;
	}


	/// @brief Resize + Fill Value: Values Not Preserved
	inline
	Chunk &
	non_preserving_resize(
		size_type const size_a,
		T const & value
	)
	{
		if ( size_ != size_a ) {
			if ( size_a > capacity_ ) {
				delete[] array_; array_ = new T[ size_a ];
				capacity_ = size_a;
			}
			size_ = size_a;
		}
		for ( size_type i = 0; i < size_; ++i ) {
			array_[ i ] = value;
		}
		return *this;
	}


	/// @brief Reserve: Values Preserved: Added Built-In Values are Not Initialized!
	inline
	Chunk &
	reserve( size_type const capacity_a )
	{
		if ( capacity_ < capacity_a ) {
			T * const new_array( new T[ capacity_a ] );
			for ( size_type i = 0; i < size_; ++i ) {
				new_array[ i ] = array_[ i ];
			}
			delete[] array_; array_ = new_array;
			capacity_ = capacity_a;
		}
		return *this;
	}


	/// @brief Shrink Capacity to Size
	inline
	Chunk &
	shrink()
	{
		if ( size_ < capacity_ ) {
			T * const new_array( ( size_ > 0 ? new T[ size_ ] : 0 ) );
			for ( size_type i = 0; i < size_; ++i ) {
				new_array[ i ] = array_[ i ];
			}
			delete[] array_; array_ = new_array;
			capacity_ = size_;
		}
		return *this;
	}


	/// @brief Swap
	inline
	void
	swap( Chunk & c )
	{
		std::swap( size_, c.size_ );
		std::swap( capacity_, c.capacity_ );
		std::swap( array_, c.array_ );
	}


	/// @brief Clear
	inline
	Chunk &
	clear()
	{
		size_ = 0;
		capacity_ = 0;
		delete[] array_; array_ = 0;
		return *this;
	}


public: // Comparison


	/// @brief Chunk == Chunk
	friend
	inline
	bool
	operator ==( Chunk const & a, Chunk const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( a.array_[ i ] != b.array_[ i ] ) return false; // Elements differ
			}
			return true; // No elements differ
		}
	}


	/// @brief Chunk != Chunk
	friend
	inline
	bool
	operator !=( Chunk const & a, Chunk const & b )
	{
		return !( a == b );
	}


public: // Friend


	/// @brief Swap
	friend
	inline
	void
	swap( Chunk & a, Chunk & b )
	{
		a.swap( b );
	}


private: // Data


	/// @brief Number of elements in use
	size_type size_;

	/// @brief Number of elements it can hold without resizing
	size_type capacity_;

	/// @brief Data array
	T * array_;


}; // Chunk


/// @brief Chunk == Chunk
template< typename T >
bool
operator ==( Chunk< T > const & a, Chunk< T > const & b );


/// @brief Chunk != Chunk
template< typename T >
bool
operator !=( Chunk< T > const & a, Chunk< T > const & b );


/// @brief Swap
template< typename T >
void
swap( Chunk< T > & a, Chunk< T > & b );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.  The legal alternative would be
// to add specializations of swap for each anticipated instantiation.


namespace std {


/// @brief std::swap( Chunk, Chunk )
template< typename T >
inline
void
swap( ObjexxFCL::Chunk< T > & a, ObjexxFCL::Chunk< T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_Chunk_HH
