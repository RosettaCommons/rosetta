#ifndef INCLUDED_ObjexxFCL_ChunkVector_hh
#define INCLUDED_ObjexxFCL_ChunkVector_hh


// ChunkVector: Chunk-Contiguous Vector for Fast Very Large Vectors
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
#include <ObjexxFCL/ChunkVector.fwd.hh>
#include <ObjexxFCL/Chunk.hh>
#include <ObjexxFCL/ChunkExponent.hh>

// C++ Headers
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>


namespace ObjexxFCL {


//// @brief ChunkVector: Chunk-Contiguous Vector for Fast Very Large Vectors
///
/// @remarks
///  @li Chunk allocation avoids large contiguous allocation failures with fragmented memory
///  @li Similar to std::deque but faster accessor performance
///  @li Construction and assignment give right-sized ChunkVector with no reserve capacity
///  @li Exponent is limited to one less than the number of bits in the Chunk size_type
///  @li Resize and push_back can add reserve capacity
///  @li Generators such as ChunkVector + ChunkVector are not provided: Can't specify the chunk exponent
///  @li Double loop operations are used instead of linear indexing for slight efficiency benefit
///
/// @note Invariants:
///  @li chunk_size_ > 0
///  @li Chunks have size == chunk_size_ except last Chunk has size in [1,chunk_size_]
///  @li Chunks have capacity == chunk_size_ except last Chunk has capacity in [chunk.size(),chunk_size_]
template< typename T >
class ChunkVector
{


private: // Friend


	template< typename > friend class ChunkVector;


public: // Types


	typedef  std::vector< Chunk< T > >  Chunks;

	// STL style
	typedef  Chunk< T >  Chunk_type;
	typedef  T  value_type;
	typedef  T &  reference;
	typedef  T const &  const_reference;
	typedef  T *  pointer;
	typedef  T const *  const_pointer;
	typedef  std::size_t  size_type;
	typedef  std::ptrdiff_t  difference_type;
	typedef  typename Chunks::size_type  Chunks_size_type;

	// C++ style
	typedef  Chunk< T >  ChunkType;
	typedef  T  Value;
	typedef  T &  Reference;
	typedef  T const &  ConstReference;
	typedef  T *  Pointer;
	typedef  T const *  ConstPointer;
	typedef  T *  Iterator;
	typedef  T const *  ConstIterator;
	typedef  std::size_t  Size;
	typedef  std::ptrdiff_t  Difference;
	typedef  typename Chunks::size_type  ChunksSize;


public: // Creation


	/// @brief Default Constructor
	inline
	ChunkVector() :
		size_( 0 ),
		chunk_exponent_( 0 ),
		chunk_size_( 1 ),
		chunk_mask_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	ChunkVector( ChunkVector const & v ) :
		size_( v.size_ ),
		chunk_exponent_( v.chunk_exponent_ ),
		chunk_size_( v.chunk_size_ ),
		chunk_mask_( v.chunk_mask_ ),
		chunks_( v.chunks_ )
	{
		assert( v.n_chunk() == computed_n_chunk() );
	}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	ChunkVector( ChunkVector< U > const & v ) :
		size_( v.size_ ),
		chunk_exponent_( v.chunk_exponent_ ),
		chunk_size_( v.chunk_size_ ),
		chunk_mask_( v.chunk_mask_ ),
		chunks_( v.n_chunk() ) // std::vector doesn't have a copy constructor template
	{
		// Size and assign Chunks
		assert( v.n_chunk() == computed_n_chunk() );
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
				Chunk_type & chunk( chunks_[ i ] );
				typename ChunkVector< U >::Chunk_type const & v_chunk( v.chunks_[ i ] );
				assert( v_chunk.size() == ( i < i_last ? chunk_size_ : computed_last_chunk_size() ) );
				chunk = v_chunk;
			}
		}
	}


	/// @brief std::vector + Exponent Constructor Template
	template< typename U, typename L >
	inline
	ChunkVector(
		std::vector< U, L > const & v,
		ChunkExponent const & chunk_exponent_a
	) :
		size_( v.size() ),
		chunk_exponent_( chunk_exponent_a ),
		chunk_size_( size_type( 1u ) << chunk_exponent_ ),
		chunk_mask_( chunk_size_ - size_type( 1u ) ),
		chunks_( computed_n_chunk() )
	{
		// Size and assign Chunks
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			typename std::vector< U, L >::const_iterator k( v.begin() );
			for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
				Chunk_type & chunk( chunks_[ i ] );
				chunk.non_preserving_resize( i < i_last ? chunk_size_ : computed_last_chunk_size() );
				for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
					chunk[ j ] = T( *k );
				}
			}
		}
	}


	/// @brief Iterator Range + Exponent Constructor Template
	template< typename InputIterator >
	inline
	ChunkVector(
		InputIterator const beg,
		InputIterator const end,
		ChunkExponent const & chunk_exponent_a
	) :
		size_( end - beg ),
		chunk_exponent_( chunk_exponent_a ),
		chunk_size_( size_type( 1u ) << chunk_exponent_ ),
		chunk_mask_( chunk_size_ - size_type( 1u ) ),
		chunks_( computed_n_chunk() )
	{
		// Size and assign Chunks
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			InputIterator k( beg );
			for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
				Chunk_type & chunk( chunks_[ i ] );
				chunk.non_preserving_resize( i < i_last ? chunk_size_ : computed_last_chunk_size() );
				for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
					chunk[ j ] = T( *k );
				}
			}
		}
	}


	/// @brief Size + Exponent Constructor: Built-In Values are Not Initialized!
	inline
	ChunkVector(
		size_type const size_a,
		ChunkExponent const & chunk_exponent_a
	) :
		size_( size_a ),
		chunk_exponent_( chunk_exponent_a ),
		chunk_size_( size_type( 1u ) << chunk_exponent_ ),
		chunk_mask_( chunk_size_ - size_type( 1u ) ),
		chunks_( computed_n_chunk() )
	{
		// Size and assign Chunks
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			for ( Chunks_size_type i = 0; i < i_last; ++i ) {
				chunks_[ i ].non_preserving_resize( chunk_size_ );
			}
			chunks_[ i_last ].non_preserving_resize( computed_last_chunk_size() );
		}
	}


	/// @brief Size + Exponent + Uniform Value Constructor
	inline
	ChunkVector(
		size_type const size_a,
		ChunkExponent const & chunk_exponent_a,
		T const & value
	) :
		size_( size_a ),
		chunk_exponent_( chunk_exponent_a ),
		chunk_size_( size_type( 1u ) << chunk_exponent_ ),
		chunk_mask_( chunk_size_ - size_type( 1u ) ),
		chunks_( computed_n_chunk() )
	{
		// Size and assign Chunks
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			for ( Chunks_size_type i = 0; i < i_last; ++i ) {
				chunks_[ i ].non_preserving_resize( chunk_size_, value );
			}
			chunks_[ i_last ].non_preserving_resize( computed_last_chunk_size(), value );
		}
	}


	/// @brief Destructor
	inline
	~ChunkVector()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	ChunkVector &
	operator =( ChunkVector const & v )
	{
		if ( this != &v ) {
			if ( chunk_exponent_ == v.chunk_exponent_ ) { // Resize for efficiency
				non_preserving_resize( v.size_ );
				if ( size_ > 0 ) {
					Chunks_size_type const i_last( i_last_chunk() );
					for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
						Chunk_type & chunk( chunks_[ i ] );
						Chunk_type const & v_chunk( v.chunks_[ i ] );
						assert( ( v_chunk.size() == chunk_size_ ) || ( i == i_last ) );
						assert( v_chunk.size() == chunk.size() );
						chunk = v_chunk;
					}
				}
			} else { // Must reallocate so use member assignment
				size_ = v.size_;
				chunk_exponent( v.chunk_exponent_ );
				chunks_ = v.chunks_;
			}
			assert( n_chunk() == computed_n_chunk() );
		}
		return *this;
	}


	/// @brief Copy Assignment Template
	template< typename U >
	inline
	ChunkVector &
	operator =( ChunkVector< U > const & v )
	{
		if ( chunk_exponent_ == v.chunk_exponent_ ) { // Resize for efficiency
			non_preserving_resize( v.size_ );
			if ( size_ > 0 ) {
				Chunks_size_type const i_last( i_last_chunk() );
				for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
					Chunk_type & chunk( chunks_[ i ] );
					typename ChunkVector< U >::Chunk_type const & v_chunk( v.chunks_[ i ] );
					assert( v_chunk.size() == ( i < i_last ? chunk_size_ : computed_last_chunk_size() ) );
					assert( v_chunk.size() == chunk.size() );
					chunk = v_chunk;
				}
			}
		} else { // Must reallocate so use member assignment
			size_ = v.size_;
			chunk_exponent( v.chunk_exponent_ );
			chunks_.clear();
			chunks_.resize( computed_n_chunk() );
			if ( size_ > 0 ) {
				Chunks_size_type const i_last( i_last_chunk() );
				for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
					Chunk_type & chunk( chunks_[ i ] );
					typename ChunkVector< U >::Chunk_type const & v_chunk( v.chunks_[ i ] );
					assert( v_chunk.size() == ( i < i_last ? chunk_size_ : computed_last_chunk_size() ) );
					chunk = v_chunk;
				}
			}
		}
		assert( n_chunk() == computed_n_chunk() );

		return *this;
	}


	/// @brief std::vector Assignment Template
	template< typename U, typename L >
	inline
	ChunkVector &
	operator =( std::vector< U, L > const & v )
	{
		non_preserving_resize( v.size() );
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			typename std::vector< U, L >::const_iterator k( v.begin() );
			for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
				Chunk_type & chunk( chunks_[ i ] );
				assert( ( chunk.size() == chunk_size_ ) || ( i == i_last ) );
				for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
					chunk[ j ] = T( *k );
				}
			}
		}
		return *this;
	}


	/// @brief std::vector + Exponent Assignment Template
	template< typename U, typename L >
	inline
	ChunkVector &
	assign(
		std::vector< U, L > const & v,
		ChunkExponent const & chunk_exponent_a
	)
	{
		if ( chunk_exponent_ == chunk_exponent_a ) { // Call other assign function for efficiency
			return operator =( v );
		} else { // Must reallocate so use member assignment
			size_ = v.size();
			chunk_exponent( chunk_exponent_a );
			chunks_.clear();
			chunks_.resize( computed_n_chunk() );
			if ( size_ > 0 ) {
				Chunks_size_type const i_last( i_last_chunk() );
				typename std::vector< U, L >::const_iterator k( v.begin() );
				for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
					Chunk_type & chunk( chunks_[ i ] );
					chunk.resize( i < i_last ? chunk_size_ : computed_last_chunk_size() );
					for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
						chunk[ j ] = T( *k );
					}
				}
			}
		}
		return *this;
	}


	/// @brief Iterator Range Assignment Template
	template< typename InputIterator >
	inline
	ChunkVector &
	assign(
		InputIterator const beg,
		InputIterator const end
	)
	{
		non_preserving_resize( end - beg );
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			InputIterator k( beg );
			for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
				Chunk_type & chunk( chunks_[ i ] );
				for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
					chunk[ j ] = T( *k );
				}
			}
		}
		return *this;
	}


	/// @brief Iterator Range + Exponent Assignment Template
	template< typename InputIterator >
	inline
	ChunkVector &
	assign(
		InputIterator const beg,
		InputIterator const end,
		ChunkExponent const & chunk_exponent_a
	)
	{
		if ( chunk_exponent_ == chunk_exponent_a ) { // Call other assign function for efficiency
			return assign( beg, end );
		} else { // Must reallocate so use member assignment
			size_ = end - beg;
			chunk_exponent( chunk_exponent_a );
			chunks_.clear();
			chunks_.resize( computed_n_chunk() );
			if ( size_ > 0 ) {
				Chunks_size_type const i_last( i_last_chunk() );
				InputIterator k( beg );
				for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
					Chunk_type & chunk( chunks_[ i ] );
					chunk.resize( i < i_last ? chunk_size_ : computed_last_chunk_size() );
					for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
						chunk[ j ] = T( *k );
					}
				}
			}
		}
		return *this;
	}


	/// @brief Size + Value Assignment
	inline
	ChunkVector &
	assign(
		size_type const size_a,
		T const & value
	)
	{
		non_preserving_resize( size_a );
		if ( size_ > 0 ) {
			Chunks_size_type const i_last( i_last_chunk() );
			for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
				chunks_[ i ] = value;
			}
		}
		return *this;
	}


	/// @brief Size + Exponent + Value Assignment
	inline
	ChunkVector &
	assign(
		size_type const size_a,
		ChunkExponent const & chunk_exponent_a,
		T const & value
	)
	{
		if ( chunk_exponent_ == chunk_exponent_a ) { // Call other assign function for efficiency
			return assign( size_a, value );
		} else { // Must reallocate so use member assignment
			size_ = size_a;
			chunk_exponent( chunk_exponent_a );
			chunks_.clear();
			chunks_.resize( computed_n_chunk() );
			if ( size_ > 0 ) {
				Chunks_size_type const i_last( i_last_chunk() );
				for ( Chunks_size_type i = 0; i <= i_last; ++i ) {
					chunks_[ i ].assign( i < i_last ? chunk_size_ : computed_last_chunk_size(), value );
				}
			}
		}
		return *this;
	}


	/// @brief += ChunkVector
	inline
	ChunkVector &
	operator +=( ChunkVector const & v )
	{
		assert( size_ == v.size_ );
		size_type k( 0 );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
				chunk[ j ] += v[ k ];
			}
		}
		return *this;
	}


	/// @brief -= ChunkVector
	inline
	ChunkVector &
	operator -=( ChunkVector const & v )
	{
		assert( size_ == v.size_ );
		size_type k( 0 );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
				chunk[ j ] -= v[ k ];
			}
		}
		return *this;
	}


	/// @brief += ChunkVector Template
	template< typename U >
	inline
	ChunkVector &
	operator +=( ChunkVector< U > const & v )
	{
		assert( size_ == v.size_ );
		typename ChunkVector< U >::size_type k( 0 );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
				chunk[ j ] += T( v[ k ] );
			}
		}
		return *this;
	}


	/// @brief -= ChunkVector Template
	template< typename U >
	inline
	ChunkVector &
	operator -=( ChunkVector< U > const & v )
	{
		assert( size_ == v.size_ );
		typename ChunkVector< U >::size_type k( 0 );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
				chunk[ j ] -= T( v[ k ] );
			}
		}
		return *this;
	}


	/// @brief += std::vector Template
	template< typename U, typename L >
	inline
	ChunkVector &
	operator +=( std::vector< U, L > const & v )
	{
		assert( size_ == v.size() );
		typename std::vector< U, L >::size_type k( 0 );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
				chunk[ j ] += T( v[ k ] );
			}
		}
		return *this;
	}


	/// @brief -= std::vector Template
	template< typename U, typename L >
	inline
	ChunkVector &
	operator -=( std::vector< U, L > const & v )
	{
		assert( size_ == v.size() );
		typename std::vector< U, L >::size_type k( 0 );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j, ++k ) {
				chunk[ j ] -= T( v[ k ] );
			}
		}
		return *this;
	}


	/// @brief = Value
	inline
	ChunkVector &
	operator =( T const & value )
	{
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			chunks_[ i ] = value;
		}
		return *this;
	}


	/// @brief += Value
	inline
	ChunkVector &
	operator +=( T const & value )
	{
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			chunks_[ i ] += value;
		}
		return *this;
	}


	/// @brief -= Value
	inline
	ChunkVector &
	operator -=( T const & value )
	{
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			chunks_[ i ] -= value;
		}
		return *this;
	}


	/// @brief *= Value
	inline
	ChunkVector &
	operator *=( T const & value )
	{
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			chunks_[ i ] *= value;
		}
		return *this;
	}


	/// @brief /= Value
	inline
	ChunkVector &
	operator /=( T const & value )
	{
		assert( value != T( 0 ) );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			chunks_[ i ] /= value;
		}
		return *this;
	}


public: // Subscript


	/// @brief ChunkVector[ i ] const: 0-Based Indexing
	inline
	T const &
	operator []( size_type const i ) const
	{
		assert( i < size_ );
		return chunks_[ i >> chunk_exponent_ ][ i & chunk_mask_ ];
	}


	/// @brief ChunkVector[ i ]: 0-Based Indexing
	inline
	T &
	operator []( size_type const i )
	{
		assert( i < size_ );
		return chunks_[ i >> chunk_exponent_ ][ i & chunk_mask_ ];
	}


	/// @brief ChunkVector( i ) const: 1-Based Indexing
	inline
	T const &
	operator ()( size_type const i ) const
	{
		assert( ( i > 0 ) && ( i <= size_ ) );
		return chunks_[ ( i - 1 ) >> chunk_exponent_ ][ ( i - 1 ) & chunk_mask_ ];
	}


	/// @brief ChunkVector( i ): 1-Based Indexing
	inline
	T &
	operator ()( size_type const i )
	{
		assert( ( i > 0 ) && ( i <= size_ ) );
		return chunks_[ ( i - 1 ) >> chunk_exponent_ ][ ( i - 1 ) & chunk_mask_ ];
	}


public: // Inspector


	/// @brief Size
	inline
	size_type
	size() const
	{
		return size_;
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


	/// @brief Chunk Exponent
	inline
	size_type
	chunk_exponent() const
	{
		return chunk_exponent_;
	}


	/// @brief Chunk Size
	inline
	size_type
	chunk_size() const
	{
		return chunk_size_;
	}


	/// @brief Number of Chunks
	inline
	Chunks_size_type
	n_chunk() const
	{
		return chunks_.size();
	}


	/// @brief First Element
	inline
	T const &
	front() const
	{
		assert( size_ > 0 );
		return chunks_[ 0 ][ 0 ];
	}


	/// @brief Last Element
	inline
	T const &
	back() const
	{
		assert( size_ > 0 );
		return operator()( size_ );
	}


	/// @brief Length
	inline
	T
	length() const
	{
		T length_sq( T( 0 ) );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j ) {
				length_sq += square( chunk[ j ] );
			}
		}
		return std::sqrt( length_sq );
	}


	/// @brief Length Squared
	inline
	T
	length_squared() const
	{
		T length_sq( T( 0 ) );
		for ( Chunks_size_type i = 0, ie = chunks_.size(); i < ie; ++i ) {
			Chunk_type & chunk( chunks_[ i ] );
			for ( size_type j = 0, je = chunk.size(); j < je; ++j ) {
				length_sq += square( chunk[ j ] );
			}
		}
		return length_sq;
	}


public: // Modifier


	/// @brief First Element
	inline
	T &
	front()
	{
		assert( size_ > 0 );
		return chunks_[ 0 ][ 0 ];
	}


	/// @brief Last Element
	inline
	T &
	back()
	{
		assert( size_ > 0 );
		return operator()( size_ );
	}


	/// @brief Append an Element
	inline
	ChunkVector &
	push_back( T const & value )
	{
		assert( size_ < max_size() );
		if ( ( size_ == 0 ) || ( last_chunk().size() == chunk_size_ ) ) { // No Chunks or last Chunk is full
			chunks_.push_back( Chunk_type() ); // Add a new Chunk
		}
		Chunk_type & chunk( last_chunk() );
		chunk.reserve( chunk_size_ ); // Reserve full size for efficient appends
		chunk.push_back( value ); // Append the new element
		++size_;
		return *this;
	}


	/// @brief Remove the Last Element
	inline
	ChunkVector &
	pop_back()
	{
		assert( size_ > 0 );
		Chunk_type & chunk( last_chunk() );
		chunk.pop_back();
		if ( chunk.empty() ) chunks_.pop_back();
		--size_;
		return *this;
	}


	/// @brief Append ChunkVector
	inline
	ChunkVector &
	append( ChunkVector const & v )
	{
		if ( v.size_ > 0 ) {
			assert( size_ <= max_size() - v.size_ );
			size_type const size_o( size_ );
			resize( size_ + v.size_ );
			size_type k( 0 );
			for ( size_type i = size_o; i < size_; ++i, ++k ) {
				(*this)[ i ] = v[ k ];
			}
		}
		return *this;
	}


	/// @brief Append ChunkVector Template
	template< typename U >
	inline
	ChunkVector &
	append( ChunkVector< U > const & v )
	{
		if ( v.size_ > 0 ) {
			assert( size_ <= max_size() - v.size_ );
			size_type const size_o( size_ );
			resize( size_ + v.size_ );
			typename ChunkVector< U >::size_type k( 0 );
			for ( size_type i = size_o; i < size_; ++i, ++k ) {
				(*this)[ i ] = T( v[ k ] );
			}
		}
		return *this;
	}


	/// @brief Append std::vector Template
	template< typename U, typename L >
	inline
	ChunkVector &
	append( std::vector< U, L > const & v )
	{
		if ( v.size() > 0 ) {
			assert( size_ <= max_size() - v.size() );
			size_type const size_o( size_ );
			resize( size_ + v.size() );
			typename std::vector< U, L >::size_type k( 0 );
			for ( size_type i = size_o; i < size_; ++i, ++k ) {
				(*this)[ i ] = T( v[ k ] );
			}
		}
		return *this;
	}


	/// @brief Resize with Same Chunk Size + Fill Value: Values Preserved
	inline
	ChunkVector &
	resize(
		size_type const size_a,
		T const & value = T()
	)
	{
		Chunks_size_type const n_chunk_o( n_chunk() );
		Chunks_size_type const n_chunk_a( ( size_a + chunk_size_ - 1 ) / chunk_size_ );
		Chunks_size_type const i_last_chunk_a( n_chunk_a - 1 );
		if ( size_a > size_ ) { // Add more values and maybe Chunks
			if ( n_chunk_a > n_chunk_o ) { // Add more Chunks: Use outer copy + inner swap for speed
				Chunks chunks_a( n_chunk_a ); // Create temporary outer vector with empty Chunks
				for ( Chunks_size_type i = 0; i < n_chunk_o; ++i ) { // Swap the old Chunks to get values
					chunks_[ i ].swap( chunks_a[ i ] );
				}
				chunks_.swap( chunks_a ); // Swap the outer vector
				for ( Chunks_size_type i = n_chunk_o - 1; i < i_last_chunk_a; ++i ) { // Fill out the Chunks
					chunks_[ i ].resize( chunk_size_, value );
				}
			}
		} else if ( size_a < size_ ) { // Remove values and maybe Chunks
			if ( n_chunk_a < n_chunk_o ) { // Remove tail Chunks
				chunks_.resize( n_chunk_a );
			}
		}
		if ( size_a > 0 ) { // Size the last Chunk
			chunks_[ i_last_chunk_a ].resize( size_a - ( i_last_chunk_a * chunk_size_ ), value );
		}
		size_ = size_a;
		return *this;
	}


	/// @brief Resize with Same Chunk Size: Values Not Preserved
	inline
	ChunkVector &
	non_preserving_resize( size_type const size_a )
	{
		Chunks_size_type const n_chunk_o( n_chunk() );
		Chunks_size_type const n_chunk_a( ( size_a + chunk_size_ - 1 ) / chunk_size_ );
		Chunks_size_type const i_last_chunk_a( n_chunk_a - 1 );
		if ( size_a > size_ ) { // Add more values and maybe Chunks
			if ( n_chunk_a > n_chunk_o ) { // Add more Chunks: Use outer copy + inner swap for speed
				Chunks chunks_a( n_chunk_a ); // Create temporary outer vector with empty Chunks
				for ( Chunks_size_type i = 0; i < n_chunk_o; ++i ) { // Swap the old Chunks to save allocation cost
					chunks_[ i ].swap( chunks_a[ i ] );
				}
				chunks_.swap( chunks_a ); // Swap the outer vector
				for ( Chunks_size_type i = n_chunk_o - 1; i < i_last_chunk_a; ++i ) { // Fill out the Chunks
					chunks_[ i ].non_preserving_resize( chunk_size_ );
				}
			}
		} else if ( size_a < size_ ) { // Remove values and maybe Chunks
			if ( n_chunk_a < n_chunk_o ) { // Remove tail Chunks
				chunks_.resize( n_chunk_a );
			}
		}
		if ( size_a > 0 ) { // Size the last Chunk
			chunks_[ i_last_chunk_a ].non_preserving_resize( size_a - ( i_last_chunk_a * chunk_size_ ) );
		}
		size_ = size_a;
		return *this;
	}


	/// @brief Reshape + Fill Value: Values Preserved
	inline
	ChunkVector &
	reshape(
		size_type const size_a,
		ChunkExponent const & chunk_exponent_a,
		T const & value = T()
	)
	{
		ChunkVector v( size_a, chunk_exponent_a, value ); // Temporary with desired shape
		for ( size_type k = 0, ke = std::min( size_, size_a ); k < ke; ++k ) { // Copy values
			v[ k ] = (*this)[ k ];
		}
		swap( v ); // Swap with temporary
		return *this;
	}


	/// @brief Reshape: Values Not Preserved
	inline
	ChunkVector &
	non_preserving_reshape(
		size_type const size_a,
		ChunkExponent const & chunk_exponent_a
	)
	{
		ChunkVector( size_a, chunk_exponent_a ).swap( *this ); // Set to new array
		return *this;
	}


	/// @brief Shrink to Right-Sized
	inline
	ChunkVector &
	shrink()
	{
		if ( size_ > 0 ) {
			Chunk_type & chunk( last_chunk() ); // Only last Chunk can have excess capacity
			if ( chunk.size() < chunk.capacity() ) { // Last Chunk has excess capacity
				chunk.shrink(); // Shrink last chunk
			}
		}
		return *this;
	}


	/// @brief Swap
	inline
	void
	swap( ChunkVector & v )
	{
		std::swap( size_, v.size_ );
		std::swap( chunk_exponent_, v.chunk_exponent_ );
		std::swap( chunk_size_, v.chunk_size_ );
		std::swap( chunk_mask_, v.chunk_mask_ );
		chunks_.swap( v.chunks_ );
	}


	/// @brief Clear
	inline
	ChunkVector &
	clear()
	{
		size_ = 0;
		chunk_exponent_ = 0;
		chunk_size_ = 1;
		chunk_mask_ = 0;
		chunks_.clear();
		return *this;
	}


	/// @brief Normalize to Unit Length
	inline
	ChunkVector &
	normalize()
	{
		T const length_( length() );
		assert( length_ > T( 0 ) );
		operator /=( length_ );
		return *this;
	}


public: // Comparison


	/// @brief Are two ChunkVectors comparable?
	friend
	inline
	bool
	comparable( ChunkVector const & a, ChunkVector const & b )
	{
		return ( a.size_ == b.size_ );
	}


	/// @brief ChunkVector == ChunkVector
	/// @note Value comparison: Chunk size ignored
	friend
	inline
	bool
	operator ==( ChunkVector const & a, ChunkVector const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] == b[ i ] ) ) return false;
			}
			return true; // No elements differ
		}
	}


	/// @brief ChunkVector != ChunkVector
	/// @note Value comparison: Chunk size ignored
	friend
	inline
	bool
	operator !=( ChunkVector const & a, ChunkVector const & b )
	{
		return !( a == b );
	}


	/// @brief ChunkVector < ChunkVector
	/// @note Value comparison: Chunk size ignored
	friend
	inline
	bool
	operator <( ChunkVector const & a, ChunkVector const & b )
	{
		if ( &a == &b ) { // Same objects
			return false;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] < b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector <= ChunkVector
	/// @note Value comparison: Chunk size ignored
	friend
	inline
	bool
	operator <=( ChunkVector const & a, ChunkVector const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] <= b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector >= ChunkVector
	/// @note Value comparison: Chunk size ignored
	friend
	inline
	bool
	operator >=( ChunkVector const & a, ChunkVector const & b )
	{
		if ( &a == &b ) { // Same objects
			return true;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] >= b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector > ChunkVector
	/// @note Value comparison: Chunk size ignored
	friend
	inline
	bool
	operator >( ChunkVector const & a, ChunkVector const & b )
	{
		if ( &a == &b ) { // Same objects
			return false;
		} else if ( a.size_ != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] > b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief Is a ChunkVector comparable to a std::vector?
	template< typename L >
	friend
	inline
	bool
	comparable( ChunkVector const & a, std::vector< T, L > const & b )
	{
		return ( a.size_ == b.size() );
	}


	/// @brief ChunkVector == std::vector Template
	template< typename L >
	friend
	inline
	bool
	operator ==( ChunkVector const & a, std::vector< T, L > const & b )
	{
		if ( a.size_ != b.size() ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] == b[ i ] ) ) return false;
			}
			return true; // No elements differ
		}
	}


	/// @brief ChunkVector != std::vector Template
	template< typename L >
	friend
	inline
	bool
	operator !=( ChunkVector const & a, std::vector< T, L > const & b )
	{
		return !( a == b );
	}


	/// @brief ChunkVector < std::vector
	template< typename L >
	friend
	inline
	bool
	operator <( ChunkVector const & a, std::vector< T, L > const & b )
	{
		if ( a.size_ != b.size() ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] < b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector <= std::vector
	template< typename L >
	friend
	inline
	bool
	operator <=( ChunkVector const & a, std::vector< T, L > const & b )
	{
		if ( a.size_ != b.size() ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] <= b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector >= std::vector
	template< typename L >
	friend
	inline
	bool
	operator >=( ChunkVector const & a, std::vector< T, L > const & b )
	{
		if ( a.size_ != b.size() ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] >= b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector > std::vector
	template< typename L >
	friend
	inline
	bool
	operator >( ChunkVector const & a, std::vector< T, L > const & b )
	{
		if ( a.size_ != b.size() ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
				if ( !( a[ i ] > b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief Is a std::vector comparable to a ChunkVector?
	template< typename L >
	friend
	inline
	bool
	comparable( std::vector< T, L > const & a, ChunkVector const & b )
	{
		return ( a.size() == b.size_ );
	}


	/// @brief std::vector == ChunkVector Template
	template< typename L >
	friend
	inline
	bool
	operator ==( std::vector< T, L > const & a, ChunkVector const & b )
	{
		return ( b == a );
	}


	/// @brief std::vector != ChunkVector Template
	template< typename L >
	friend
	inline
	bool
	operator !=( std::vector< T, L > const & a, ChunkVector const & b )
	{
		return !( b == a );
	}


	/// @brief std::vector < ChunkVector
	template< typename L >
	friend
	inline
	bool
	operator <( std::vector< T, L > const & a, ChunkVector const & b )
	{
		if ( a.size() != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
				if ( !( a[ i ] < b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief std::vector <= ChunkVector
	template< typename L >
	friend
	inline
	bool
	operator <=( std::vector< T, L > const & a, ChunkVector const & b )
	{
		if ( a.size() != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
				if ( !( a[ i ] <= b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief std::vector >= ChunkVector
	template< typename L >
	friend
	inline
	bool
	operator >=( std::vector< T, L > const & a, ChunkVector const & b )
	{
		if ( a.size() != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
				if ( !( a[ i ] >= b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief std::vector > ChunkVector
	template< typename L >
	friend
	inline
	bool
	operator >( std::vector< T, L > const & a, ChunkVector const & b )
	{
		if ( a.size() != b.size_ ) { // Sizes differ
			return false;
		} else { // Compare values
			for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
				if ( !( a[ i ] > b[ i ] ) ) return false;
			}
			return true;
		}
	}


	/// @brief ChunkVector == T
	friend
	inline
	bool
	operator ==( ChunkVector const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( a[ i ] != t ) return false;
		}
		return true;
	}


	/// @brief ChunkVector != T
	friend
	inline
	bool
	operator !=( ChunkVector const & a, T const & t )
	{
		return !( a == t );
	}


	/// @brief ChunkVector < T
	friend
	inline
	bool
	operator <( ChunkVector const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a[ i ] < t ) ) return false;
		}
		return true;
	}


	/// @brief ChunkVector <= T
	friend
	inline
	bool
	operator <=( ChunkVector const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a[ i ] <= t ) ) return false;
		}
		return true;
	}


	/// @brief ChunkVector >= T
	friend
	inline
	bool
	operator >=( ChunkVector const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a[ i ] >= t ) ) return false;
		}
		return true;
	}


	/// @brief ChunkVector > T
	friend
	inline
	bool
	operator >( ChunkVector const & a, T const & t )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( a[ i ] > t ) ) return false;
		}
		return true;
	}


	/// @brief T == ChunkVector
	friend
	inline
	bool
	operator ==( T const & t, ChunkVector const & a )
	{
		return ( a == t );
	}


	/// @brief T != ChunkVector
	friend
	inline
	bool
	operator !=( T const & t, ChunkVector const & a )
	{
		return !( a == t );
	}


	/// @brief T < ChunkVector
	friend
	inline
	bool
	operator <( T const & t, ChunkVector const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t < a[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T <= ChunkVector
	friend
	inline
	bool
	operator <=( T const & t, ChunkVector const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t <= a[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T >= ChunkVector
	friend
	inline
	bool
	operator >=( T const & t, ChunkVector const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t >= a[ i ] ) ) return false;
		}
		return true;
	}


	/// @brief T > ChunkVector
	friend
	inline
	bool
	operator >( T const & t, ChunkVector const & a )
	{
		for ( size_type i = 0, ie = a.size_; i < ie; ++i ) {
			if ( !( t > a[ i ] ) ) return false;
		}
		return true;
	}


public: // Friend


	/// @brief Dot Product
	friend
	inline
	T
	dot_product( ChunkVector const & a, ChunkVector const & b )
	{
		assert( a.size() == b.size() );
		T sum( T( 0 ) );
		for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
			sum += a[ i ] * b[ i ];
		}
		return sum;
	}


	/// @brief Dot Product
	friend
	inline
	T
	dot( ChunkVector const & a, ChunkVector const & b )
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
	distance( ChunkVector const & a, ChunkVector const & b )
	{
		assert( a.size() == b.size() );
		T distance_sq( T( 0 ) );
		for ( size_type i = 0, ie = a.size(); i < ie; ++i ) {
			distance_sq += square( a[ i ] - b[ i ] );
		}
		return std::sqrt( distance_sq );
	}


	/// @brief Distance Squared
	friend
	inline
	T
	distance_squared( ChunkVector const & a, ChunkVector const & b )
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
	swap( ChunkVector & a, ChunkVector & b )
	{
		a.swap( b );
	}


private: // Functions


	/// @brief Index of Last Chunk
	inline
	Chunks_size_type
	i_last_chunk() const
	{
		assert( ! chunks_.empty() );
		return chunks_.size() - 1;
	}


	/// @brief Last Chunk
	inline
	Chunk_type const &
	last_chunk() const
	{
		assert( ! chunks_.empty() );
		return chunks_[ chunks_.size() - 1 ];
	}


	/// @brief Last Chunk
	inline
	Chunk_type &
	last_chunk()
	{
		assert( ! chunks_.empty() );
		return chunks_[ chunks_.size() - 1 ];
	}


	/// @brief Computed Number of Chunks
	inline
	Chunks_size_type
	computed_n_chunk() const
	{
		return ( size_ + chunk_size_ - 1 ) / chunk_size_;
	}


	/// @brief Computed Last Chunk Size
	inline
	size_type
	computed_last_chunk_size() const
	{
		assert( size_ > 0 );
		return size_ - ( ( ( size_ - 1 ) / chunk_size_ ) * chunk_size_ );
	}


	/// @brief Exponent
	inline
	ChunkVector &
	chunk_exponent( ChunkExponent const & chunk_exponent_a )
	{
		chunk_exponent_ = chunk_exponent_a;
		chunk_size_ = size_type( 1u ) << chunk_exponent_;
		chunk_mask_ = chunk_size_ - size_type( 1u );
		return *this;
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


	/// @brief Number of elements
	size_type size_;

	/// @brief Chunk size exponent (< number of bits in size_type)
	size_type chunk_exponent_;

	/// @brief Chunk size (a power of 2) (last Chunk can be smaller)
	size_type chunk_size_;

	/// @brief Chunk index identification mask
	size_type chunk_mask_;

	/// @brief Vector of Chunks
	Chunks chunks_;


}; // ChunkVector


/// @brief Are two ChunkVectors comparable?
template< typename T >
bool
comparable( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector == ChunkVector
/// @note Value comparison: Chunk size ignored
template< typename T >
bool
operator ==( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector != ChunkVector
/// @note Value comparison: Chunk size ignored
template< typename T >
bool
operator !=( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector < ChunkVector
/// @note Value comparison: Chunk size ignored
template< typename T >
bool
operator <( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector <= ChunkVector
/// @note Value comparison: Chunk size ignored
template< typename T >
bool
operator <=( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector >= ChunkVector
/// @note Value comparison: Chunk size ignored
template< typename T >
bool
operator >=( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector > ChunkVector
/// @note Value comparison: Chunk size ignored
template< typename T >
bool
operator >( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief Is a ChunkVector comparable to a std::vector?
template< typename T, typename L >
bool
comparable( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief ChunkVector == std::vector Template
template< typename T, typename L >
bool
operator ==( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief ChunkVector != std::vector Template
template< typename T, typename L >
bool
operator !=( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief ChunkVector < std::vector
template< typename T, typename L >
bool
operator <( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief ChunkVector <= std::vector
template< typename T, typename L >
bool
operator <=( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief ChunkVector >= std::vector
template< typename T, typename L >
bool
operator >=( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief ChunkVector > std::vector
template< typename T, typename L >
bool
operator >( ChunkVector< T > const & a, std::vector< T, L > const & b );


/// @brief Is a std::vector comparable to a ChunkVector?
template< typename T, typename L >
bool
comparable( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief std::vector == ChunkVector Template
template< typename T, typename L >
bool
operator ==( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief std::vector != ChunkVector Template
template< typename T, typename L >
bool
operator !=( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief std::vector < ChunkVector
template< typename T, typename L >
bool
operator <( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief std::vector <= ChunkVector
template< typename T, typename L >
bool
operator <=( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief std::vector >= ChunkVector
template< typename T, typename L >
bool
operator >=( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief std::vector > ChunkVector
template< typename T, typename L >
bool
operator >( std::vector< T, L > const & a, ChunkVector< T > const & b );


/// @brief ChunkVector == T
template< typename T >
bool
operator ==( ChunkVector< T > const & a, T const & t );


/// @brief ChunkVector != T
template< typename T >
bool
operator !=( ChunkVector< T > const & a, T const & t );


/// @brief ChunkVector < T
template< typename T >
bool
operator <( ChunkVector< T > const & a, T const & t );


/// @brief ChunkVector <= T
template< typename T >
bool
operator <=( ChunkVector< T > const & a, T const & t );


/// @brief ChunkVector >= T
template< typename T >
bool
operator >=( ChunkVector< T > const & a, T const & t );


/// @brief ChunkVector > T
template< typename T >
bool
operator >( ChunkVector< T > const & a, T const & t );


/// @brief T == ChunkVector
template< typename T >
bool
operator ==( T const & t, ChunkVector< T > const & a );


/// @brief T != ChunkVector
template< typename T >
bool
operator !=( T const & t, ChunkVector< T > const & a );


/// @brief T < ChunkVector
template< typename T >
bool
operator <( T const & t, ChunkVector< T > const & a );


/// @brief T <= ChunkVector
template< typename T >
bool
operator <=( T const & t, ChunkVector< T > const & a );


/// @brief T >= ChunkVector
template< typename T >
bool
operator >=( T const & t, ChunkVector< T > const & a );


/// @brief T > ChunkVector
template< typename T >
bool
operator >( T const & t, ChunkVector< T > const & a );


/// @brief Dot Product
template< typename T >
T
dot_product( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief Dot Product
template< typename T >
T
dot( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief Distance
template< typename T >
T
distance( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief Distance Squared
template< typename T >
T
distance_squared( ChunkVector< T > const & a, ChunkVector< T > const & b );


/// @brief Swap
template< typename T >
void
swap( ChunkVector< T > & a, ChunkVector< T > & b );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.  The legal alternative would be
// to add specializations of swap for each anticipated instantiation.


namespace std {


/// @brief std::swap( ChunkVector, ChunkVector )
template< typename T >
inline
void
swap( ObjexxFCL::ChunkVector< T > & a, ObjexxFCL::ChunkVector< T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_ChunkVector_HH
