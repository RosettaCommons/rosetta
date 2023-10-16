// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/PseudoPerplexityMetric.fwd.hh
/// @brief A class for calculating the pseudo-perplexity from a given PerResidueProbabilitiesMetric.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)


#ifndef INCLUDED_protocols_esm_perplexity_PseudoPerplexityMetric_fwd_hh
#define INCLUDED_protocols_esm_perplexity_PseudoPerplexityMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace esm_perplexity {

class PseudoPerplexityMetric;

using PseudoPerplexityMetricOP = utility::pointer::shared_ptr< PseudoPerplexityMetric >;
using PseudoPerplexityMetricCOP = utility::pointer::shared_ptr< PseudoPerplexityMetric const >;

} //esm_perplexity
} //protocols

#endif //INCLUDED_protocols_esm_perplexity_PseudoPerplexityMetric_fwd_hh
