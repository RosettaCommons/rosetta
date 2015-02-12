// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/symmetry/SymmetricAddMembraneMover.fwd.hh
///
/// @brief      Add Membrane Representation to a Symmetric starting pose
/// @details    Given a symmetrized pose, add the membrane components,
///             ensuring all data descriptions are for the asymmetric subunit
///             (in terms of scoring) and the membrane is the root of the
///             whole system. After applying SymmetricAddMembraneMover
///             pose.conformaiton().is_membrane() AND is_symmetric( pose )
///             should both return true
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/9/15)

#ifndef INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_fwd_hh
#define INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace symmetry {
    
class SymmetricAddMembraneMover;
typedef utility::pointer::shared_ptr< SymmetricAddMembraneMover > SymmetricAddMembraneMoverOP;
typedef utility::pointer::shared_ptr< SymmetricAddMembraneMover const > SymmetricAddMembraneMoverCOP;

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_fwd_hh
