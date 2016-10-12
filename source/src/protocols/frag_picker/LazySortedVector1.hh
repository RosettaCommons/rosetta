// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_LazySortedVector1_hh
#define INCLUDED_protocols_frag_picker_LazySortedVector1_hh

// package headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <algorithm>

namespace protocols {
namespace frag_picker {

template<class T, class StrictWeakOrdering>
class LazySortedVector1 : public utility::pointer::ReferenceCount {
public:
	LazySortedVector1(StrictWeakOrdering cmp, core::Size sorted_capacity,core::Size max_capacity) :
		comp(cmp) {
		max_capacity_ = max_capacity;
		sorted_ = false;
		sorted_capacity_ = sorted_capacity;
		last_ = 1;
		data_.resize( max_capacity );
		n_inserts_ = 0;
		n_denied_ = 0;
		n_sorts_ = 0;
	}

	~LazySortedVector1() override {

		//     std::cerr<<n_inserts_<<" "<<n_denied_<<" "<<n_sorts_<<'\n';
	}

	void resize(core::Size sorted_capacity,core::Size max_capacity) {

		max_capacity_ = max_capacity;
		sorted_capacity_ = sorted_capacity;
		data_.resize( max_capacity );
	}

	inline const T& top() {

		if ( !sorted_ ) {
			std::sort(data_.begin(), data_.begin() + last_ - 1, comp);
		}
		return data_[1];
	}

	inline bool push(const T& x) {

		if ( comp(worst_,x) ) {
			n_denied_++;
			return false;
		}

		if ( last_ > max_capacity_ ) {
			std::sort(data_.begin(), data_.begin() + last_ - 1, comp);
			last_ = sorted_capacity_ + 1;
			worst_ = data_[sorted_capacity_];

			n_sorts_++;
		}
		n_inserts_++;
		data_[last_] = x;
		last_++;
		return true;
	}

	inline void set_worst( T new_worst) { worst_ = new_worst; }
	inline T peek_back() const { return data_[sorted_capacity_]; }
	inline T peek_front() const { return data_[1]; }

	inline core::Size count_inserted() const { return data_.size(); }


	/// @brief sets new capacity for the container
	inline void set_boundary(core::Size sorted_capacity) {
		sorted_capacity_ = sorted_capacity;
	}

	inline T& at(core::Size index) {
		if ( !sorted_ ) {
			std::sort(data_.begin(), data_.begin() + last_ - 1, comp);
		}
		return data_.at(index);
	}

	inline T& operator[](core::Size index) {
		if ( !sorted_ ) {
			std::sort(data_.begin(), data_.begin() + last_ - 1, comp);
		}
		return data_[index];
	}

	inline core::Size size() const {
		return std::min(sorted_capacity_,last_);
	}

	utility::vector1<T>& expose_data() {

		good_data_.clear();
		if ( last_<=2 ) {
			return good_data_;
		}
		good_data_.reserve(size());
		if ( !sorted_ ) {
			std::sort(data_.begin(), data_.begin() + last_ - 1, comp);
		}

		//  std::copy(data_.begin(), data_.begin() + size() - 1, good_data_.begin() );
		for ( core::Size i=1; i<=size(); i++ ) {
			good_data_.push_back(data_[i]);
		}

		return good_data_;
	}

	inline void clear() {
		data_.clear();
	}
private:
	core::Size last_;
	T worst_;
	StrictWeakOrdering comp;
	bool sorted_;
	core::Size max_capacity_;
	core::Size sorted_capacity_;
	utility::vector1<T> data_;
	utility::vector1<T> good_data_;
	// debug info:
	core::Size n_inserts_;
	core::Size n_denied_;
	core::Size n_sorts_;
};

} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_LazySortedVector1_hh
