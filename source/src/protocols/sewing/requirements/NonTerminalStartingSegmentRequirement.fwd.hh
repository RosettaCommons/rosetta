// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/NonTerminalStartingSegmentRequirement.fwd.hh
/// @brief a Requirement that an Assembly have less than a certain number of clashes
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirement_fwd_hh
#define INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirement_fwd_hh

// Utility headers
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace sewing {
namespace requirements {

class NonTerminalStartingSegmentRequirement;

typedef utility::pointer::shared_ptr< NonTerminalStartingSegmentRequirement > NonTerminalStartingSegmentRequirementOP;
typedef utility::pointer::shared_ptr< NonTerminalStartingSegmentRequirement const > NonTerminalStartingSegmentRequirementCOP;



} //protocols
} //sewing
} //requirements


#endif //INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirement_fwd_hh





