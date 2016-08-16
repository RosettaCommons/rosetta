// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/SegmentSwap.fwd.hh
/// @brief forward declaration for SegmentSwap
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_SegmentSwap_fwd_hh
#define INCLUDED_protocols_forge_build_SegmentSwap_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for SegmentSwap
class SegmentSwap;


/// @brief owning pointer for SegmentSwap
typedef utility::pointer::shared_ptr< SegmentSwap > SegmentSwapOP;


/// @brief owning pointer for const SegmentSwap
typedef utility::pointer::shared_ptr< SegmentSwap const > SegmentSwapCOP;


/// @brief access pointer for SegmentSwap
typedef utility::pointer::weak_ptr< SegmentSwap > SegmentSwapAP;


/// @brief access pointer for const SegmentSwap
typedef utility::pointer::weak_ptr< SegmentSwap const > SegmentSwapCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_SegmentSwap_FWD_HH */
