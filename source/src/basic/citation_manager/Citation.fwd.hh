// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/Citation.fwd.hh
/// @brief Data structure for storing a citation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_basic_citation_manager_Citation_fwd_hh
#define INCLUDED_basic_citation_manager_Citation_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace citation_manager {

/// @brief The format for writing a citation.
/// @details If this is altered, update Citation::get_formatted_citation().
enum class CitationFormat {
	DefaultStyle = 1, //Keep first
	NatureStyle, //Keep second-to-last
	end_of_list = NatureStyle //Keep last
};

class AuthorNames;

typedef utility::pointer::shared_ptr< AuthorNames > AuthorNamesOP;
typedef utility::pointer::shared_ptr< AuthorNames const > AuthorNamesCOP;

class Citation;

typedef utility::pointer::shared_ptr< Citation > CitationOP;
typedef utility::pointer::shared_ptr< Citation const > CitationCOP;

} //basic
} //citation_manager

#endif //INCLUDED_basic_citation_manager_Citation_fwd_hh
