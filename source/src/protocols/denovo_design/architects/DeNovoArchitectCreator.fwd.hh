// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/denovo_design/architects/DeNovoArchitectCreator.fwd.hh
/// @brief  Forward declaration of a class that instantiates a particular DeNovoArchitect
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_protocols_denovo_design_architects_DeNovoArchitectCreator_FWD_HH
#define INCLUDED_protocols_denovo_design_architects_DeNovoArchitectCreator_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace denovo_design {
namespace architects {

class DeNovoArchitectCreator;

typedef utility::pointer::shared_ptr< DeNovoArchitectCreator > DeNovoArchitectCreatorOP;
typedef utility::pointer::shared_ptr< DeNovoArchitectCreator const > DeNovoArchitectCreatorCOP;

} //namespace architects
} //namespace denovo_design
} //namespace protocols

#endif
