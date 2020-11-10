// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationCollectionBase.cc
/// @brief Base structure for storing citations.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <basic/citation_manager/CitationCollectionBase.hh>

namespace basic {
namespace citation_manager {

utility::vector1< CitationCollectionBaseCOP > const &
CitationCollectionList::citations() const {
	return entries_;
}

bool
CitationCollectionList::empty() const {
	return entries_.empty();
}

void
CitationCollectionList::add( CitationCollectionBaseCOP const & cite ) {
	for ( CitationCollectionBaseCOP const & entry: entries_ ) {
		if ( *entry == *cite ) {
			return; // We already have it.
		}
	}
	// We didn't find it - add it.
	entries_.push_back( cite );
}

void
CitationCollectionList::add( CitationCollectionList const & other ) {
	add( other.entries_ );
}

} //basic
} //citation_manager

