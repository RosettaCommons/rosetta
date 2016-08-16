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

#ifndef INCLUDED_protocols_comparative_modeling_AlignmentSet_HH
#define INCLUDED_protocols_comparative_modeling_AlignmentSet_HH

#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

class AlignmentSet {

public:
	AlignmentSet () {}
	virtual ~AlignmentSet() {}

	/// @brief Gets the next alignment from the stream. Might be the same
	/// alignment every time, or might be a different alignment after each
	/// call to this function.
	utility::vector1< core::sequence::SequenceAlignment >
	alignments() const;

	void insert( core::sequence::SequenceAlignment aln );

	core::Size size() const;

private:
	std::set< core::sequence::SequenceAlignment > align_set_;

}; // AlignmentSet

} // comparative_modeling
} // protocols

#endif
