// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/constraint_generator/HydrogenBondConstraintGenerator.fwd.hh
/// @brief
/// @author Tom Linsky ( tlinsky at uw dot edu )


#ifndef INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGenerator_fwd_hh
#define INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGenerator_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace constraint_generator {

class HydrogenBondConstraintGenerator;

typedef utility::pointer::shared_ptr< HydrogenBondConstraintGenerator > HydrogenBondConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< HydrogenBondConstraintGenerator const > HydrogenBondConstraintGeneratorCOP;

typedef utility::vector1< HydrogenBondConstraintGeneratorOP > HydrogenBondConstraintGeneratorOPs;
typedef utility::vector1< HydrogenBondConstraintGeneratorCOP > HydrogenBondConstraintGeneratorCOPs;

class HBondData;

typedef utility::pointer::shared_ptr< HBondData > HBondDataOP;
typedef utility::pointer::shared_ptr< HBondData const > HBondDataCOP;

typedef utility::vector1< HBondDataOP > HBondDataOPs;
typedef utility::vector1< HBondDataCOP > HBondDataCOPs;

}
}

#endif
