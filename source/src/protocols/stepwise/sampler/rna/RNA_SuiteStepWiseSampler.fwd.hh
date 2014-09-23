// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_SuiteStepWiseSampler.fwd.hh
/// @brief Generate rotamers for one RNA suite (from residue i to i+1).
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_sampler_rna_RNA_SuiteStepWiseSampler_fwd_HH
#define INCLUDED_protocols_sampler_rna_RNA_SuiteStepWiseSampler_fwd_HH


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

class RNA_SuiteStepWiseSampler;
typedef utility::pointer::shared_ptr< RNA_SuiteStepWiseSampler > RNA_SuiteStepWiseSamplerOP;

} //rna
} //sampler
} //stepwise
} //protocols

#endif
