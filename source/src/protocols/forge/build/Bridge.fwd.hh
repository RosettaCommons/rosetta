// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/build/Bridge.fwd.hh
/// @brief  forward declaration for Bridge
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_Bridge_fwd_hh
#define INCLUDED_protocols_forge_build_Bridge_fwd_hh


// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for Bridge
class Bridge;


/// @brief owning pointer for Bridge
typedef utility::pointer::shared_ptr< Bridge > BridgeOP;


/// @brief owning pointer for const Bridge
typedef utility::pointer::shared_ptr< Bridge const > BridgeCOP;


/// @brief access pointer for Bridge
typedef utility::pointer::weak_ptr< Bridge > BridgeAP;


/// @brief access pointer for const Bridge
typedef utility::pointer::weak_ptr< Bridge const > BridgeCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_Bridge_FWD_HH */
