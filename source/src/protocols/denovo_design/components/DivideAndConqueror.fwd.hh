// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/DivideAndConqueror.fwd.hh
/// @brief Splits a denovo structure into pieces, and devises a strategy for folding the structure piece-by-piece
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_DivideAndConqueror_fwd_hh
#define INCLUDED_protocols_denovo_design_components_DivideAndConqueror_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <set>

// Forward
namespace protocols {
namespace denovo_design {
namespace components {

class DivideAndConqueror;

typedef utility::pointer::shared_ptr< DivideAndConqueror > DivideAndConquerorOP;
typedef utility::pointer::shared_ptr< DivideAndConqueror const > DivideAndConquerorCOP;

class BuildPhases;

typedef utility::pointer::shared_ptr< BuildPhases > BuildPhasesOP;
typedef utility::pointer::shared_ptr< BuildPhases const > BuildPhasesCOP;

} //protocols
} //denovo_design
} //components


#endif //INCLUDED_protocols_denovo_design_components_DivideAndConqueror_fwd_hh


