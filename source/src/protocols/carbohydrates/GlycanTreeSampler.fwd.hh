// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeSampler.fwd.hh
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_carbohydrates_GlycanTreeSampler_fwd_hh
#define INCLUDED_protocols_carbohydrates_GlycanTreeSampler_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace carbohydrates {

class GlycanTreeSampler;

typedef utility::pointer::shared_ptr< GlycanTreeSampler > GlycanTreeSamplerOP;
typedef utility::pointer::shared_ptr< GlycanTreeSampler const > GlycanTreeSamplerCOP;



} //protocols
} //carbohydrates


#endif //INCLUDED_protocols/carbohydrates_GlycanTreeSampler_fwd_hh





