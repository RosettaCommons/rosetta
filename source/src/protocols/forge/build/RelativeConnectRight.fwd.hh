// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/RelativeConnectRight.fwd.hh
/// @brief forward declaration for RelativeConnectRight
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_RelativeConnectRight_fwd_hh
#define INCLUDED_protocols_forge_build_RelativeConnectRight_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for RelativeConnectRight
class RelativeConnectRight;


/// @brief owning pointer for RelativeConnectRight
typedef utility::pointer::shared_ptr< RelativeConnectRight > RelativeConnectRightOP;


/// @brief owning pointer for const RelativeConnectRight
typedef utility::pointer::shared_ptr< RelativeConnectRight const > RelativeConnectRightCOP;


/// @brief access pointer for RelativeConnectRight
typedef utility::pointer::weak_ptr< RelativeConnectRight > RelativeConnectRightAP;


/// @brief access pointer for const RelativeConnectRight
typedef utility::pointer::weak_ptr< RelativeConnectRight const > RelativeConnectRightCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_RelativeConnectRight_FWD_HH */
