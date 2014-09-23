// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/build/RelativeSequencePosition.fwd.hh
/// @brief forward declaration for RelativeSequencePosition
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_RelativeSequencePosition_fwd_hh
#define INCLUDED_protocols_forge_build_RelativeSequencePosition_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for RelativeSequencePosition
struct RelativeSequencePosition;
/// @brief owning pointer for RelativeSequencePosition
typedef utility::pointer::shared_ptr< RelativeSequencePosition > RelativeSequencePositionOP;
/// @brief owning pointer for const RelativeSequencePosition
typedef utility::pointer::shared_ptr< RelativeSequencePosition const > RelativeSequencePositionCOP;
/// @brief access pointer for RelativeSequencePosition
typedef utility::pointer::weak_ptr< RelativeSequencePosition > RelativeSequencePositionAP;
/// @brief access pointer for const RelativeSequencePosition
typedef utility::pointer::weak_ptr< RelativeSequencePosition const > RelativeSequencePositionCAP;


/// @brief forward declaration for CountFromLeft
struct CountFromLeft;
/// @brief owning pointer for CountFromLeft
typedef utility::pointer::shared_ptr< CountFromLeft > CountFromLeftOP;
/// @brief owning pointer for const CountFromLeft
typedef utility::pointer::shared_ptr< CountFromLeft const > CountFromLeftCOP;
/// @brief access pointer for CountFromLeft
typedef utility::pointer::weak_ptr< CountFromLeft > CountFromLeftAP;
/// @brief access pointer for const CountFromLeft
typedef utility::pointer::weak_ptr< CountFromLeft const > CountFromLeftCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_RelativeSequencePosition_FWD_HH */
