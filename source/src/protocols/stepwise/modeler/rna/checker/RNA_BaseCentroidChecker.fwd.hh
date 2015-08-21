// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
//  vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/StepWise_RNA_BaseCentroidChecker
/// @brief  Various OPS for various classes.
/// @author Rhiju Das

#ifndef INCLUDED_protocols_stepwise_rna_RNA_BaseCentroidChecker_FWD_HH
#define INCLUDED_protocols_stepwise_rna_RNA_BaseCentroidChecker_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

class RNA_BaseCentroidChecker;
typedef utility::pointer::shared_ptr< RNA_BaseCentroidChecker > RNA_BaseCentroidCheckerOP;
typedef utility::pointer::shared_ptr< RNA_BaseCentroidChecker const > RNA_BaseCentroidCheckerCOP;

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#endif
