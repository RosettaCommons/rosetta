// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/SegmentRebuild.fwd.hh
/// @brief forward declaration for SegmentRebuild
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_SegmentRebuild_fwd_hh
#define INCLUDED_protocols_forge_build_SegmentRebuild_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for SegmentRebuild
class SegmentRebuild;


/// @brief owning pointer for SegmentRebuild
typedef utility::pointer::shared_ptr< SegmentRebuild > SegmentRebuildOP;


/// @brief owning pointer for const SegmentRebuild
typedef utility::pointer::shared_ptr< SegmentRebuild const > SegmentRebuildCOP;


/// @brief access pointer for SegmentRebuild
typedef utility::pointer::weak_ptr< SegmentRebuild > SegmentRebuildAP;


/// @brief access pointer for const SegmentRebuild
typedef utility::pointer::weak_ptr< SegmentRebuild const > SegmentRebuildCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_SegmentRebuild_FWD_HH */
