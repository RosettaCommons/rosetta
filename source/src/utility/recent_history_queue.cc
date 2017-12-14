// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/recent_history_queue.cc
/// @brief  A queue for holding the history in which certain members of a set
///         are promoted to the front of the queue before eventually falling
///         off the end of the queue.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/recent_history_queue.hh>


namespace utility {

recent_history_queue::recent_history_queue() :
	num_elements_( 0 ),
	history_size_( 0 ),
	curr_num_in_recent_history_( 0 ),
	head_of_rh_queue_ptr_( 0 ),
	end_of_rh_queue_ptr_( 0 )
{}

recent_history_queue::recent_history_queue( int num_elements, int history_size ) :
	num_elements_( num_elements ),
	history_size_( history_size ),
	curr_num_in_recent_history_( 0 ),
	recent_history_queue_( history_size ),
	head_of_rh_queue_ptr_( 0 ),
	end_of_rh_queue_ptr_( 0 ),
	elements_2_recent_history_( num_elements, 0 )
{
	history_queue_struct zero;
	zero.more_recent_ptr = 0;
	zero.element_in_rh = 0;
	zero.more_ancient_ptr = 0;
	for ( int ii = 1; ii <= history_size_; ++ii ) {
		recent_history_queue_[ ii ] = zero;
	}
}

recent_history_queue::~recent_history_queue() = default;

void recent_history_queue::clear() {
	curr_num_in_recent_history_ = 0;
	history_queue_struct zero;
	zero.more_recent_ptr = 0;
	zero.element_in_rh = 0;
	zero.more_ancient_ptr = 0;
	for ( int ii = 1; ii <= history_size_; ++ii ) {
		recent_history_queue_[ ii ] = zero;
	}
	head_of_rh_queue_ptr_ = 0;
	end_of_rh_queue_ptr_ = 0;
	for ( int ii = 1; ii <= num_elements_; ++ii ) {
		elements_2_recent_history_[ ii ] = 0;
	}
}

void recent_history_queue::num_elements( int num_elements )
{
	num_elements_ = num_elements;
	elements_2_recent_history_.resize( num_elements_ );
	clear();
}

void recent_history_queue::history_size( int history_size )
{
	history_size_ = history_size;
	recent_history_queue_.resize( history_size_ );
	clear();
}

int recent_history_queue::num_elements() const { return num_elements_; }
int recent_history_queue::history_size() const { return history_size_; }

int recent_history_queue::curr_num_in_recent_history() const
{
	return curr_num_in_recent_history_;
}

int recent_history_queue::head_of_queue() const { return head_of_rh_queue_ptr_; }
int recent_history_queue::end_of_queue() const { return end_of_rh_queue_ptr_; }


int recent_history_queue::pos_in_history_queue( int element ) const {
	return elements_2_recent_history_[ element ];
}

int recent_history_queue::push_to_front_of_history_queue( int element )
{
	if ( elements_2_recent_history_[ element ] != 0 ) {
		//already stored in recent history -- nothing gets bumped
		int const element_rh_id = elements_2_recent_history_[ element ];
		if ( head_of_rh_queue_ptr_ == element_rh_id ) return 0;

		int const anc_id = recent_history_queue_[ element_rh_id ].more_ancient_ptr;
		int const rec_id = recent_history_queue_[ element_rh_id ].more_recent_ptr;

		recent_history_queue_[ rec_id ].more_ancient_ptr = anc_id;
		if ( element_rh_id != end_of_rh_queue_ptr_ ) {
			recent_history_queue_[ anc_id ].more_recent_ptr = rec_id;
		} else {
			end_of_rh_queue_ptr_ = rec_id;
		}

		recent_history_queue_[ head_of_rh_queue_ptr_ ].more_recent_ptr = element_rh_id;
		recent_history_queue_[ element_rh_id ].more_ancient_ptr = head_of_rh_queue_ptr_;
		recent_history_queue_[ element_rh_id ].more_recent_ptr = 0;
		head_of_rh_queue_ptr_ = element_rh_id;

		return 0;
	} else if ( curr_num_in_recent_history_ < history_size_ ) {
		++curr_num_in_recent_history_;
		if ( curr_num_in_recent_history_ == 1 ) end_of_rh_queue_ptr_ = 1;

		elements_2_recent_history_[ element ] = curr_num_in_recent_history_;
		recent_history_queue_[ curr_num_in_recent_history_ ].element_in_rh = element;

		if ( head_of_rh_queue_ptr_ != 0 ) {
			recent_history_queue_[ head_of_rh_queue_ptr_ ].more_recent_ptr = curr_num_in_recent_history_;
		}
		recent_history_queue_[ curr_num_in_recent_history_ ].more_ancient_ptr = head_of_rh_queue_ptr_;
		recent_history_queue_[ curr_num_in_recent_history_ ].more_recent_ptr = 0;
		head_of_rh_queue_ptr_ = curr_num_in_recent_history_;

	} else {
		//not in recent history, something gets bumped

		if ( history_size_ == 1 ) {
			elements_2_recent_history_[ recent_history_queue_[ 1 ].element_in_rh ] = 0;
			elements_2_recent_history_[ element ] = 1;
			recent_history_queue_[ 1 ].element_in_rh = element;
			return 1;
		} else {
			int const prev_end = end_of_rh_queue_ptr_;
			int const one_before_end = recent_history_queue_[ prev_end ].more_recent_ptr;
			recent_history_queue_[ one_before_end ].more_ancient_ptr = 0;
			end_of_rh_queue_ptr_ = one_before_end;

			elements_2_recent_history_[ recent_history_queue_[ prev_end ].element_in_rh ] = 0;

			elements_2_recent_history_[ element ] = prev_end;
			recent_history_queue_[ prev_end ].element_in_rh = element;

			recent_history_queue_[ prev_end ].more_recent_ptr = 0;
			recent_history_queue_[ prev_end ].more_ancient_ptr = head_of_rh_queue_ptr_;
			recent_history_queue_[ head_of_rh_queue_ptr_ ].more_recent_ptr = prev_end;
			head_of_rh_queue_ptr_ = prev_end;

			return prev_end;
		}
	}
	return 0; //control of flow never reaches here; makes compiler happy

}

unsigned int
recent_history_queue::dynamic_memory_usage() const {
	unsigned int total( 0 );
	total += recent_history_queue_.size() * sizeof( history_queue_struct );
	total += elements_2_recent_history_.size() * sizeof( int );
	return total;
}

}

