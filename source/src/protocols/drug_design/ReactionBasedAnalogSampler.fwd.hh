// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionBasedAnalogSampler.hh
/// @brief apply RDKit's reaction mechanism to react reagents based on a given reaction
/// @author Tracy Tang (yidan.tang@vanderbilt.edu)

#ifndef INCLUDED_protocols_drug_design_ReactionBasedAnalogSampler_fwd_hh
#define INCLUDED_protocols_drug_design_ReactionBasedAnalogSampler_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace drug_design {


class ReactionBasedAnalogSampler; // fwd declaration
typedef utility::pointer::shared_ptr< ReactionBasedAnalogSampler > ReactionBasedAnalogSamplerOP;
typedef utility::pointer::shared_ptr< ReactionBasedAnalogSampler const > ReactionBasedAnalogSamplerCOP;


} // namespace drug_design
} // namespace protocols

#endif // INCLUDED_protocols_drug_design_ReactionBasedAnalogSampler_fwd_hh



