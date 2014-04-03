// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_BoundedPriorityQueue_hh
#define INCLUDED_protocols_frag_picker_BoundedPriorityQueue_hh

// package headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <algorithm>

//Auto Headers


namespace protocols {
namespace frag_picker {

template<class T, class StrictWeakOrdering>
class BoundedPriorityQueue : public utility::pointer::ReferenceCount {
public:
	BoundedPriorityQueue(StrictWeakOrdering cmp, Size max_capacity) :
		comp(cmp) {
		max_capacity_ = max_capacity;
		sorted = false;
		initialized_ = false;
		last_ = 0;

		n_sorts = 0;
		n_denied = 0;
	}

	~BoundedPriorityQueue() {
	    std::cerr<<"sorts, denied: "<<n_sorts<<" "<<n_denied<<"\n";
	}

	inline const T& top() {
		BoundedPriorityQueue<T, StrictWeakOrdering>::lazy_sort();
		return data_.front();
	}

	inline bool push(const T& x) {

		if ((last_ == max_capacity_) && (comp(worst_, x))) {
			n_denied++;
			return false;
		}
		if (!initialized_) {
			worst_ = x;
			initialized_ = true;
			last_ = 1;
			data_.push_back(x);
			sorted = false;
			return true;
		}

		if (last_ < max_capacity_) {
			last_++;
			data_.push_back(x);
			sorted = false;
			return true;
		}
		data_[last_] = x;
		sorted = false;
		BoundedPriorityQueue<T, StrictWeakOrdering>::lazy_sort();

		worst_ = data_.back();
		sorted = true;

		return true;
	}

	inline Size count_inserted() { return last_; }

	/// @brief Removes and returns a reference to the best element in the queue.
	/// @detailed This element, when compared to any other element in the queue will give true.
	inline T& pop() {

		T & t = BoundedPriorityQueue<T, StrictWeakOrdering>::top();
		data_.pop_back();
		return t;
	}

	inline void lazy_sort() {
		if (!sorted) {
			std::sort(data_.begin(), data_.end(), comp);
			sorted = true;
			n_sorts++;
		}
	}

	/// @brief sets new capacity for the container
	inline void set_boundary(Size max_capacity) {
		max_capacity_ = max_capacity;
	}

	/// @brief Returns a reference to the worst element in the queue.
	/// @detailed This element when compared to any other element in the queue will give false.
	inline T& peek_back() {
		BoundedPriorityQueue<T, StrictWeakOrdering>::lazy_sort();
		return data_.back();
	}

	inline T& peek_front() {
		BoundedPriorityQueue<T, StrictWeakOrdering>::lazy_sort();
		return data_.front();
	}

	inline T& at(Size index) {
		return data_.at(index);
	}

	inline T& operator[](Size index) {
		return data_[index];
	}

	inline Size size() {
		return data_.size();
	}

	inline utility::vector1<T>& expose_data() {
		return data_;
	}

	inline void clear() {
		data_.clear();
	}
private:
	StrictWeakOrdering comp;
	bool sorted;
	Size max_capacity_;
	T worst_;
	Size last_;
	bool initialized_;
	utility::vector1<T> data_;

	// debug info
	Size n_sorts;
	Size n_denied;
};

} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_BoundedPriorityQueue_hh
