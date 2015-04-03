// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/ArrayPool.hh
/// @brief  ArrayPool class declaration and implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_graph_ArrayPool_hh
#define INCLUDED_core_graph_ArrayPool_hh

/// Project Headers
#include <platform/types.hh>

/// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

/// C++ headers
#include <list>

#include <utility/vector1.hh>


namespace core {
namespace graph {

/// @brief Class Array0 is a c-style array wrapper that does bounds checking in debug mode.
/// It indexes from 0 just like regular c-arrays.  Class Array0 does not manage it's own
/// memory.  It does not allocate memory if you want to make it larger, nor does it
/// deallocate memory when you destroy it.  Bounds checking only ensures that the user
/// does not go outside of the memory Array0 thinks it's in charge of.  If the user
/// should happen to point the array0 at memory that has not been allocated, Array0
/// is not responsible for segmentation fault that will likely occur.  Garbage in,
/// garbage out.
template< class T >
class Array0 {
public:
	/// @brief Default ctor points at null.
	Array0() :
		array_( 0 ),
		size_( 0 )
	{}

	/// @brief Array and size constructor -- point this Array0 at a block of memory
	Array0( T * mem_begin, platform::Size size ) :
		array_( mem_begin ),
		size_( size )
	{}

	/// @brief Copy constructor -- point this Array0 at a block of memory
	Array0( Array0< T > const & other ) :
		array_( other.array_ ),
		size_( other.size_ )
	{}

	/// @brief Assignment operator -- point this Array0 at a different block of memory.
	Array0< T > const &
	operator = ( Array0< T > const & rhs ) {
		if ( &rhs != this ) {
			array_ = rhs.array_;
			size_ = rhs.size_;
		}
		return *this;
	}

	/// @brief The destructor does not deallocate the memory that this Array0 points at.
	/// That is the responsibility of some other class.  Array0 is for bounds checking only.
	~Array0() {}

public:
	/// Accessors and mutators.

	T & operator [] ( platform::Size index ) {
	debug_assert( bounds_check( index ) );
		return array_[ index ];
	}

	T const & operator [] ( platform::Size index ) const {
	debug_assert( bounds_check( index ) );
		return array_[ index ];
	}

	T & operator [] ( int index ) {
	debug_assert( bounds_check( index ) );
		return array_[ index ];
	}

	T const & operator [] ( int index ) const {
	debug_assert( bounds_check( index ) );
		return array_[ index ];
	}

	platform::Size size() const { return size_; }

private:

	bool
	bounds_check( platform::Size index ) const {
		return index < size_;
	}

	bool
	bounds_check( int index ) const {
		return index < (int) size_ && index >= 0;
	}

private:
	T * array_;
	platform::Size size_;

};


/// @brief NegSpaceElement represents a single element in the singly-linked
/// list of negative space in an array pool.
template < class T >
class NegSpaceElement
{
public:
	NegSpaceElement() : next_( 0 ), array_( 0 ), allocated_( false )
	{}

	NegSpaceElement( NegSpaceElement * next,  T * array ) :
		next_( next ),
		array_( array ),
		allocated_( false )
	{}

	void set_next( NegSpaceElement * next ) {
		next_ = next;
	}

	void set_array( T * array ) {
		array_ = array;
	}

	void set_allocated( bool setting ) {
		allocated_ = setting;
	}


	NegSpaceElement * next() const { return next_; }
	T * array() const { return array_; }
	bool allocated() const { return allocated_; }

	/// @brief Remove my next element from negative space,
	/// and return a pointer to it.  Maintain negative space
	/// integrity by pointing my next_ pointer at the next_
	/// pointer of the removed element
	NegSpaceElement * pop() {
		NegSpaceElement< T > * removed = next_;
		next_ = removed->next_;
		return removed;
	}

	/// @brief Add an element to negative space by inserting it
	/// behind this element.
	void insert_after( NegSpaceElement * element ) {
		element->next_ = next_;
		next_ = element;
	}

private:

	NegSpaceElement * next_;
	T * array_;
	bool allocated_;

};

template < class T >
class ArrayPoolElement
{
public:

	ArrayPoolElement( platform::Size size, NegSpaceElement< T > * neg_ptr ) :
		array_( neg_ptr->array(), size ),
		neg_ptr_( neg_ptr )
	{}

	ArrayPoolElement( ArrayPoolElement< T > const & other ) :
		array_( other.array_ ),
		neg_ptr_( other.neg_ptr_ )
	{}

	~ArrayPoolElement() {}

	T const & operator [] ( platform::Size index ) const {
	debug_assert( neg_ptr_->allocated() );
		return array_[ index ];
	}

	T & operator [] ( platform::Size index ) {
	debug_assert( neg_ptr_->allocated() );
		return array_[ index ];
	}

	T const & operator[] ( int index ) const {
	debug_assert( neg_ptr_->allocated() );
		return array_[ index ];
	}

	T & operator[] ( int index ) {
	debug_assert( neg_ptr_->allocated() );
		return array_[ index ];
	}

	bool valid() const {
		return neg_ptr_->allocated();
	}

	void copy_array_contents( ArrayPoolElement< T > const & other ) {
	debug_assert( array_.size() == other.array_.size() );
		for ( platform::Size ii = 0; ii < array_.size(); ++ii ) {
			array_[ ii ] = other.array_[ ii ];
		}
	}

	template< typename > friend class ArrayPool;

private:
	Array0< T > array_;
	NegSpaceElement< T > * neg_ptr_;

};

template < class T >
class ArrayPool : public utility::pointer::ReferenceCount
{
public:
	/// Creation and destruction

	/// @brief Default constructor uses a block size of 32
	ArrayPool() :
		block_size_( 32 ),
		array_size_( 0 ),
		nblocks_( 0 ),
		nnegative_( 0 )
	{}

	/// @brief Constructor with block-size specification.
	ArrayPool( platform::Size block_size ) :
		block_size_( block_size ),
		array_size_( 0 ),
		nblocks_( 0 ),
		nnegative_( 0 )
	{}

	/// @brief Constructor with block-size and array-size specification.
	ArrayPool( platform::Size block_size, platform::Size array_size ) :
		block_size_( block_size ),
		array_size_( array_size ),
		nblocks_( 0 ),
		nnegative_( 0 )
	{}

	~ArrayPool() {
		if ( !empty() ) {
			utility_exit_with_message( "Error in ArrayPool destructor: cannot free a non-empty ArrayPool" );
		}
		clear();
	}

public:
	/// Methods to read information about the pool size

	/// @brief Returns the size of each array to be allocated.
	platform::Size array_size() const {
		return array_size_;
	}

	/// @brief Returns the number of arrays to allocate in each block
	platform::Size block_size() const {
		return block_size_;
	}

	/// @brief Returns the number of bytes occupied by all allocated and
	/// not-yet allocated arrays.
	//
	/// @details Approximate the cost of each list element in the two std::lists
	/// as 4 times the cost of a pointer.
	platform::Size footprint() const {
		return nblocks_ * ( block_size_ * (sizeof( NegSpaceElement< T > * ) + array_size_ * sizeof( T )) + 8 * sizeof( platform::Size ))
			+ sizeof( ArrayPool< T > );
	}

public:
	/// Methods to modify the size of the pool

	/// @brief Set the size of the arrays that the pool is meant to allocate.
	///
	/// @details The
	void set_array_size( platform::Size size ) {
		if ( !empty() ) {
			utility_exit_with_message( "ERROR: ArrayPool array size cannot be changed unless the ArrayPool is empty" );
		}
		if ( array_size_ != size ) {
			clear();
			array_size_ = size;
		}
	}

	void set_block_size( platform::Size size ) {
		if ( !empty() ) {
			utility_exit_with_message( "ERROR: ArrayPool block size cannot be changed unless the ArrayPool is empty" );
		}
		if ( block_size_ != size ) {
			clear();
			block_size_ = size;
		}
	}

	/// @brief Return the number of allocated array blocks
	platform::Size nblocks() const {
		return nblocks_;
	}

	/// @brief Return the number of allocated arrays.  If this were multiplied by the
	/// array_size_, then you'd have the amount of space spent representing the T arrays
	/// but you wouldn't have the amount of space spent representing the negative space
	/// linked list.
	platform::Size nallocated() const {
		return block_size_ * nblocks_;
	}

	/// @brief Return the number of allocated arrays that have been requested by the
	///  user and are now "outstanding" in that they have to be returned before the
	///  array pool may be deleted.
	platform::Size noutstanding() const {
		return nallocated() - nnegative_;
	}

	bool empty() const {
		return nnegative_ == nallocated();
	}

public:
	/// Methods for allocating and deallocating arrays from this pool

	/// @brief Request a new ArrayPoolElement from the pool.  This method
	/// will enlarge the pool if it has grown to its capacity.
	ArrayPoolElement< T >
	new_array() {
		if ( neg_begin_.next() == 0 ) create_new_block();
		NegSpaceElement< T > * first = neg_begin_.pop();
	debug_assert( first != 0 );
	debug_assert( first->allocated() == false );

		first->set_allocated( true );
		--nnegative_;

		return ArrayPoolElement< T >( array_size_, first );
	}

	/// @brief Free an ArrayPoolElement to the pool.  The element is invalid
	/// after this call completes.
	void
	deallocate_array(
		ArrayPoolElement< T > const & element
	)
	{
	debug_assert( mine( element.neg_ptr_ ) );

		NegSpaceElement< T > * neg_element = element.neg_ptr_;
	debug_assert( neg_element != 0 );
	debug_assert( neg_element->allocated() );

		neg_begin_.insert_after( neg_element );
		neg_element->set_allocated( false );
		++nnegative_;

	debug_assert( nnegative_ <= nallocated() );
	}

private:

	/// @brief Add a new block to the pool and add the elements of this block to the
	/// negative space singly linked list.  The block sizes and the array sizes must
	/// be set before this method is called and, so long as the ArrayPool is not empty,
	/// these sizes may not be altered.
	///
	/// @details This array pool supports, but does not recommend using, array sizes of 0.
	/// The resulting ArrayPoolElements have Array0 objecst that point at null, and therefore,
	/// should never be dereferenced.  However, there is nothing inherrently wrong with pointing
	/// at null, or with having an array that is never examined, and so it is supported.
	void create_new_block()
	{
	debug_assert( neg_begin_.next() == 0 );
		if ( array_size_ != 0 ) {

			NegSpaceElement< T > * const neg_block = new NegSpaceElement< T >[ block_size_ ];

			if ( neg_block == 0 ) {
				utility_exit_with_message( "ERROR: new failed in ArrayPool when requesting neg space block of size " +
					utility::to_string( block_size_ ) + " (" +
					utility::to_string( block_size_ * sizeof( NegSpaceElement< T > )) +
					" bytes)" );
			}

			T * const t_block = new T[ block_size_ * array_size_ ];

			if ( t_block == 0 ) {
				utility_exit_with_message( "ERROR: new failed in ArrayPool when requesting Array block of size " +
					utility::to_string( block_size_ ) + "x" + utility::to_string( array_size_ ) + " (" +
					utility::to_string( block_size_ * array_size_ * sizeof( T )) +
					" bytes)" );
			}

			NegSpaceElement< T > * neg_iter  = neg_block;
			T * t_iter                       = t_block;

			for ( platform::Size ii = 1; ii <= block_size_; ++ii ) {
				neg_iter->set_array( t_iter );
				t_iter += array_size_;
				neg_iter += 1;
			}
			neg_iter = neg_block;
			for ( platform::Size ii = 1; ii < block_size_; ++ii ) {
				NegSpaceElement< T > * last = neg_iter;
				++neg_iter;
				last->set_next( neg_iter );
			}
			neg_iter->set_next( 0 );
			neg_begin_.set_next( neg_block );

			++nblocks_;
			nnegative_ += block_size_;
			neg_space_blocks_.push_back( neg_block );
			array_blocks_.push_back( t_block );
		} else {
			/// 0-size array support.
			NegSpaceElement< T > * const neg_block = new NegSpaceElement< T >[ block_size_ ];

			if ( neg_block == 0 ) {
				utility_exit_with_message( "ERROR: new failed in ArrayPool when requesting neg space block of size " +
					utility::to_string( block_size_ ) + " (" +
					utility::to_string( block_size_ * sizeof( NegSpaceElement< T > )) +
					" bytes)" );
			}

			T * const t_block = 0;
			NegSpaceElement< T > * neg_iter  = neg_block;

			for ( platform::Size ii = 1; ii <= block_size_; ++ii ) {
				neg_iter->set_array( t_block );
				neg_iter += 1;
			}
			neg_iter = neg_block;
			for ( platform::Size ii = 1; ii < block_size_; ++ii ) {
				NegSpaceElement< T > * last = neg_iter;
				++neg_iter;
				last->set_next( neg_iter );
			}
			neg_iter->set_next( 0 );
			neg_begin_.set_next( neg_block );

			++nblocks_;
			nnegative_ += block_size_;
			neg_space_blocks_.push_back( neg_block );
			array_blocks_.push_back( t_block );
		}
	}

	/// @brief Deallocate all the allocated blocks.  The pool must be empty before
	/// this should be called or dangling references are likely.
	void clear()
	{
	debug_assert( empty() );
		for ( typename std::list< NegSpaceElement< T > * >::const_iterator
				iter = neg_space_blocks_.begin(),
				end_iter = neg_space_blocks_.end();
				iter != end_iter; ++iter ) {
			delete [] (*iter );
		}
		neg_space_blocks_.clear();

		for ( typename std::list< T * >::const_iterator
				iter = array_blocks_.begin(),
				end_iter = array_blocks_.end();
				iter != end_iter; ++iter ) {
			delete [] (*iter );
		}
		array_blocks_.clear();

		nblocks_ = 0;
		nnegative_ = 0;
		neg_begin_.set_next( 0 );
	}


	/// @brief Determine if a given pointer to a negative-space element belongs to this
	/// pool.  The compiler cannot ensure that when a ArrayPoolElement is handed to
	/// the deallocate_array method that it belongs to the pool being invoked upon.
	///
	/// @details This method is too slow to include in release builds, but is here to
	/// ensure in debug mode that two pools don't accidentally start crossing their
	/// pointers and "sharing" data, since that would be disasterous.
	bool mine( NegSpaceElement< T > const * neg_element ) const {
		for ( typename std::list< NegSpaceElement< T > * >::const_iterator
				iter = neg_space_blocks_.begin(),
				end_iter = neg_space_blocks_.end();
				iter != end_iter; ++iter ) {
			if ( neg_element >= (*iter ) && neg_element < (*iter) + block_size_ ) {
				return true;
			}
		}
		/// This neg space element must belong to some other pool.
		return false;
	}

private:

	platform::Size block_size_;
	platform::Size array_size_;
	platform::Size nblocks_;
	platform::Size nnegative_; // the number of arrays which have been allocated but not requested by the user
	NegSpaceElement< T > neg_begin_;

	std::list< NegSpaceElement< T > * > neg_space_blocks_;
	std::list< T * > array_blocks_;

};


}
}

#endif
