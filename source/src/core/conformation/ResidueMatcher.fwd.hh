// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Residue.fwd.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_ResidueMatcher_fwd_hh
#define INCLUDED_core_conformation_ResidueMatcher_fwd_hh


// Project headers

// Utility headers
//#include <utility/vector1.fwd.hh>
//#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers


namespace core {
namespace conformation {

class ResidueMatcher;

typedef  utility::pointer::shared_ptr< ResidueMatcher >  ResidueMatcherOP;
typedef  utility::pointer::shared_ptr< ResidueMatcher const >  ResidueMatcherCOP;

} // namespace conformation
} // namespace core


#endif // INCLUDED_core_conformation_Residue_FWD_HH
