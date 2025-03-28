// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationCollection.fwd.hh
/// @brief A class for keeping track of a collection of citations for a particular Rosetta module.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_basic_citation_manager_CitationCollection_fwd_hh
#define INCLUDED_basic_citation_manager_CitationCollection_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace citation_manager {

class CitationCollection;

typedef utility::pointer::shared_ptr< CitationCollection > CitationCollectionOP;
typedef utility::pointer::shared_ptr< CitationCollection const > CitationCollectionCOP;

} //basic
} //citation_manager

#endif //INCLUDED_basic_citation_manager_CitationCollection_fwd_hh
