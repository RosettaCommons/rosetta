// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_noesy_assign_ResonanceList_fwd_hh
#define INCLUDED_protocols_noesy_assign_ResonanceList_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace noesy_assign {

// Forward
class ResonanceList;

// Types
typedef  utility::pointer::shared_ptr< ResonanceList >  ResonanceListOP;
typedef  utility::pointer::shared_ptr< ResonanceList const >  ResonanceListCOP;

typedef  utility::pointer::weak_ptr< ResonanceList >  ResonanceListAP;
typedef  utility::pointer::weak_ptr< ResonanceList const >  ResonanceListCAP;

}
}

#endif
