// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeRelax.fwd.hh
/// @brief A protocol for optimizing glycan trees using GlycanRelax from the base of the tree out to the leaves.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_carbohydrates_GlycanTreeRelax_fwd_hh
#define INCLUDED_protocols_carbohydrates_GlycanTreeRelax_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace carbohydrates {

class GlycanTreeRelax;

typedef utility::pointer::shared_ptr< GlycanTreeRelax > GlycanTreeRelaxOP;
typedef utility::pointer::shared_ptr< GlycanTreeRelax const > GlycanTreeRelaxCOP;

} //protocols
} //carbohydrates

#endif //INCLUDED_protocols_carbohydrates_GlycanTreeRelax_fwd_hh
