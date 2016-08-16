// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/bb_sampler/SmallBBSampler.fwd.hh
/// @brief A bb sampler that samples within a range of a starting angle.  Similar to small mover.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_bb_sampler_SmallBBSampler_fwd_hh
#define INCLUDED_protocols_simple_moves_bb_sampler_SmallBBSampler_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace simple_moves {
namespace bb_sampler {

class SmallBBSampler;

typedef utility::pointer::shared_ptr< SmallBBSampler > SmallBBSamplerOP;
typedef utility::pointer::shared_ptr< SmallBBSampler const > SmallBBSamplerCOP;



} //protocols
} //simple_moves
} //bb_sampler


#endif //INCLUDED_protocols_simple_moves_bb_sampler_SmallBBSampler_fwd_hh





