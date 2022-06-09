// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/VBWrapper.hh
/// @brief  A class derived from VirtualBase for wrapping a primative
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_VBWrapper_hh
#define INCLUDED_utility_VBWrapper_hh

// Unit headers
#include <utility/VBWrapper.hh>

// Package Headers
#include <utility/VirtualBase.hh>
#include <utility/pointer/owning_ptr.hh>

namespace utility {

template <class T>
class VBWrapper: public VirtualBase
{
public:
	VBWrapper(T val) : payload_(val) {}
	VBWrapper(VBWrapper const & ) = default;

	T const &
	data() const {
		return payload_;
	}

	void data(T const & new_val) {
		payload_ = new_val;
	}


private:
	T payload_;

}; // VirtualBase

template <class T>
utility::pointer::shared_ptr<VBWrapper<T>>
make_shared_vb_wrapper(T val) {
	return utility::pointer::make_shared<utility::VBWrapper<T>>(val);
}

} // namespace utility

#endif // INCLUDED_utility_VirtualBase_hh
