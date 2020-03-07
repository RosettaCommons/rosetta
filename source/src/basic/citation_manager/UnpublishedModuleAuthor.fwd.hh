// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/UnpublishedModuleAuthor.fwd.hh
/// @brief Authorship information for a single author in an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_citation_manager_UnpublishedModuleAuthor_fwd_hh
#define INCLUDED_basic_citation_manager_UnpublishedModuleAuthor_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace citation_manager {

class UnpublishedModuleAuthor;

using UnpublishedModuleAuthorOP = utility::pointer::shared_ptr< UnpublishedModuleAuthor >;
using UnpublishedModuleAuthorCOP = utility::pointer::shared_ptr< UnpublishedModuleAuthor const >;

} //citation_manager
} //basic

#endif //INCLUDED_basic_citation_manager_UnpublishedModuleAuthor_fwd_hh
