// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.cc
/// @brief A class for predicting amino acid probabilities P(Sequence|Sequence) using the ESM language model.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

#ifndef INCLUDED_protocols_esm_perplexity_EsmPerplexityTensorflowProtocol_fwd_hh
#define INCLUDED_protocols_esm_perplexity_EsmPerplexityTensorflowProtocol_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace esm_perplexity {

class EsmPerplexityTensorflowProtocol;

using EsmPerplexityTensorflowProtocolOP = utility::pointer::shared_ptr< EsmPerplexityTensorflowProtocol >;
using EsmPerplexityTensorflowProtocolCOP = utility::pointer::shared_ptr< EsmPerplexityTensorflowProtocol const >;

} //esm_perplexity
} //protocols

#endif //INCLUDED_protocols_esm_perplexity_EsmPerplexityTensorflowProtocol_fwd_hh
