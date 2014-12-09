// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/ScoringSecMatchRPE.fwd.hh
/// @brief
/// @author Kui Chan, kuichan@u.washington.edu, Oct 2009

#ifndef INCLUDED_protocols_match_downstream_ScoringSecMatchRPE_fwd_hh
#define INCLUDED_protocols_match_downstream_ScoringSecMatchRPE_fwd_hh

// Unit headers
// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace downstream {

class ScoringSecMatchRPE;
typedef utility::pointer::shared_ptr< ScoringSecMatchRPE > ScoringSecMatchRPEOP;
typedef utility::pointer::shared_ptr< ScoringSecMatchRPE const > ScoringSecMatchRPECOP;


}
}
}

#endif
