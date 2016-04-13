// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/power_diagram/PowerDiagram.fwd.hh
/// @brief  core::scoring::power_diagram::PowerDiagram forward declaration
/// @author Jim Havranek

#ifndef INCLUDED_core_scoring_power_diagram_PowerDiagram_fwd_hh
#define INCLUDED_core_scoring_power_diagram_PowerDiagram_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace power_diagram {

class PDsphere;
typedef utility::pointer::shared_ptr< PDsphere > PDsphereOP;
typedef utility::pointer::shared_ptr< PDsphere const > PDsphereCOP;

class PDvertex;
typedef utility::pointer::shared_ptr< PDvertex > PDvertexOP;
typedef utility::pointer::shared_ptr< PDvertex const > PDvertexCOP;

class PDinter;
typedef utility::pointer::shared_ptr< PDinter > PDinterOP;
typedef utility::pointer::shared_ptr< PDinter const > PDinterCOP;

// forward declaration for PowerDiagram
class PowerDiagram;
typedef utility::pointer::shared_ptr< PowerDiagram > PowerDiagramOP;
typedef utility::pointer::shared_ptr< PowerDiagram const > PowerDiagramCOP;

}
}
}


#endif /*INCLUDED_core_scoring_power_diagram_PowerDiagram_fwd_hh*/
