// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/AdaptiveScoreHistogram.hh
/// @brief A very basic histogram to keep statistics of scores.
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_AdaptiveScoreHistogram_hh
#define INCLUDED_protocols_frag_picker_scores_AdaptiveScoreHistogram_hh

// package headers
#include <protocols/frag_picker/scores/AdaptiveScoreHistogram.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>

#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

class AdaptiveScoreHistogram : public utility::pointer::ReferenceCount {
public:
	AdaptiveScoreHistogram(Real,Real);
	virtual ~AdaptiveScoreHistogram();

	void insert(Real);

	inline Size at(Size index) {
		return data_.at(index);
	}

	inline Size operator[](Size index) {
		return data_[index];
	}

	inline Size size() {
		return data_.size();
	}

	inline utility::vector1<Size>& expose_counts() {
		return data_;
	}

	inline void clear() {
		data_.clear();
	}

	Size sum();

	Real p_value(Real);
private:
	Real bin_size_;
	bool is_up_to_date_;
	utility::vector1<Size> data_;
	utility::vector1<Size> cumulative_sums_;
};
} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_AdaptiveScoreHistogram_HH */
