// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationManager.fwd.hh
/// @brief A class that receives lists of works to cite from Rosetta modules, then returns a list of all
/// works to cite on demand.  Threadsafe.
/// @details For works with publications that have DOIs, this loads a list of references from the database,
/// indexed by DOI, once in a threadsafe manner on object creation.  This allows modules to only specify
/// the DOI.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_citation_manager_CitationManager_fwd_hh
#define INCLUDED_basic_citation_manager_CitationManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace citation_manager {

class CitationManager;

using CitationManagerOP = utility::pointer::shared_ptr< CitationManager >;
using CitationManagerCOP = utility::pointer::shared_ptr< CitationManager const >;

} //citation_manager
} //basic

#endif //INCLUDED_basic_citation_manager_CitationManager_fwd_hh
