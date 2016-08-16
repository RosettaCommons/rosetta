// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_noesy_assign_Resonance_fwd_hh
#define INCLUDED_protocols_noesy_assign_Resonance_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace noesy_assign {

// Forward
class Resonance;

// Types
typedef  utility::pointer::shared_ptr< Resonance >  ResonanceOP;
typedef  utility::pointer::shared_ptr< Resonance const >  ResonanceCOP;

typedef  utility::pointer::weak_ptr< Resonance >  ResonanceAP;
typedef  utility::pointer::weak_ptr< Resonance const >  ResonanceCAP;

}
}

#endif
