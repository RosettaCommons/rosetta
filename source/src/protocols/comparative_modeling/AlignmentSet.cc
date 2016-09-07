// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/AlignmentSet.hh
/// @brief
/// @author James Thompson

#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <protocols/comparative_modeling/AlignmentSet.hh>

#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

/// @brief Gets the next alignment from the stream. Might be the same
/// alignment every time, or might be a different alignment after each
/// call to this function.
utility::vector1< core::sequence::SequenceAlignment >
AlignmentSet::alignments() const {
	utility::vector1< core::sequence::SequenceAlignment > alignments_;
	using std::set;
	using core::sequence::SequenceAlignment;
	for ( auto copy : align_set_ ) {
		alignments_.push_back( copy );
	}
	return alignments_;
}

void AlignmentSet::insert(
	core::sequence::SequenceAlignment aln
) {
	align_set_.insert( aln );
}

core::Size AlignmentSet::size() const {
	return align_set_.size();
}

} // comparative_modeling
} // protocols
