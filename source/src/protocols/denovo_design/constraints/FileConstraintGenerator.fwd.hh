// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/constraints/FileConstraintGenerator.fwd.hh
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @author Tom Linsky ( tlinsky at uw dot edu ), Mar 2016

#ifndef INCLUDED_protocols_forge_constraints_FileConstraintGenerator_fwd_hh
#define INCLUDED_protocols_forge_constraints_FileConstraintGenerator_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace denovo_design {
namespace constraints {

class FileConstraintGenerator;

typedef utility::pointer::shared_ptr< FileConstraintGenerator > FileConstraintGeneratorOP;
typedef utility::pointer::weak_ptr< FileConstraintGenerator const > FileConstraintGeneratorCAP;

} //namespace constraints
} //namespace denovo_design
} //namespace protocols

#endif // INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_fwd_hh
