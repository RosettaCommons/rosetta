// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/ConnectRight.fwd.hh
/// @brief forward declaration for ConnectRight
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_ConnectRight_fwd_hh
#define INCLUDED_protocols_forge_build_ConnectRight_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for ConnectRight
class ConnectRight;


/// @brief owning pointer for ConnectRight
typedef utility::pointer::shared_ptr< ConnectRight > ConnectRightOP;


/// @brief owning pointer for const ConnectRight
typedef utility::pointer::shared_ptr< ConnectRight const > ConnectRightCOP;


/// @brief access pointer for ConnectRight
typedef utility::pointer::weak_ptr< ConnectRight > ConnectRightAP;


/// @brief access pointer for const ConnectRight
typedef utility::pointer::weak_ptr< ConnectRight const > ConnectRightCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_ConnectRight_FWD_HH */
