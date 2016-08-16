// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/recent_history_queue.hh
/// @brief  A queue for holding the history in which certain members of a set
///         are promoted to the front of the queue before eventually falling
///         off the end of the queue.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_recent_history_queue_hh
#define INCLUDED_utility_recent_history_queue_hh

// Unit headers
//#include <utility/recent_history_queue.fwd.hh>

// ObjexxFCL Headers
#include <utility/vector1.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace utility {


struct history_queue_struct
{
	int more_recent_ptr;
	int element_in_rh;
	int more_ancient_ptr;

	/// WARNING WARNING WARNING:
	/// if you going to modify this struct make sure to adjust definition of comparison operator below.
};

// Have to define it by hand because compiler does not automatically generate comparison for structs (and we need it for vector1::has_one implementaion)
inline bool operator== (history_queue_struct const &a, history_queue_struct const &b)
{
	return (a.more_recent_ptr  == b.more_recent_ptr &&
		a.element_in_rh    == b.element_in_rh &&
		a.more_ancient_ptr == b.more_ancient_ptr);
}


/// @brief A class for keeping track of a subset of elements in a set that are
/// pushed into a queue in a certain order, and which fall off the end of the queue
/// in ther order in which they arrive.  Elements in the set can be bumped to the
/// front of the queue.
///
/// The queue is "in place", so there are no calls to new or delete with repeated
/// calls to push_to_front_of_history_queue().
///
/// The position in queue can be used to keep track of data for elements, where
/// the position is understood not as the number of elements between the element
/// and the front of the queue, but rather, an index that the object has in the
/// "recent_history_queue_" data member -- an array.  If an element is in the
/// queue with index X, and then it is pushed to the front of the history queue,
/// its index will still be X.
class recent_history_queue : public utility::pointer::ReferenceCount {
public:
	recent_history_queue();
	recent_history_queue( int num_elements, int history_size );
	virtual ~recent_history_queue();

	void clear();
	void num_elements( int num_elements );
	void history_size( int history_size );

	int num_elements() const;
	int history_size() const;
	int curr_num_in_recent_history() const;

	/// @brief For unit-testing purposes
	int head_of_queue() const;
	/// @brief For unit-testing purposes
	int end_of_queue() const;

	/// @brief Returns the position in the recent history for a given element in the set.  Returns 0
	/// if the element is not part of the recent history.
	int pos_in_history_queue( int element ) const;

	/// @brief Push an element to the front of the recent history queue.  This will likely bump
	/// an element that had been in the recent history queue; in that event, this
	/// function returns the in-place index for that bumped element.  If the new
	/// element doesn't displace some other element (i.e. if the queue is either not
	/// yet full, or if the element was already in the queue), then this function returns
	/// a fictional index of "0".
	int push_to_front_of_history_queue( int element );

	unsigned int dynamic_memory_usage() const;


private:
	int num_elements_;
	int history_size_;
	int curr_num_in_recent_history_;
	utility::vector1< history_queue_struct > recent_history_queue_;
	int head_of_rh_queue_ptr_;
	int end_of_rh_queue_ptr_;
	utility::vector1< int > elements_2_recent_history_;


}; // class recent_history_queue

}

#endif
