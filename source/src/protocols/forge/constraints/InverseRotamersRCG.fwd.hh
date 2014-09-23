// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/constraints/InverseRotamersRCG.fwd.hh
///
/// @author Florian Richter, floric@u.washington.edu, april 2010

#ifndef INCLUDED_protocols_forge_constraints_InverseRotamersRCG_fwd_hh
#define INCLUDED_protocols_forge_constraints_InverseRotamersRCG_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace forge{
namespace constraints{

class InverseRotamersRCG;

typedef utility::pointer::shared_ptr< InverseRotamersRCG > InverseRotamersRCGOP;
typedef utility::pointer::weak_ptr< InverseRotamersRCG const > InverseRotamersRCGCAP;

}
}
}
#endif
