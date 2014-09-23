// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/components/BDR.fwd.hh
/// @brief  forward declaration for BDR
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_components_BDR_fwd_hh
#define INCLUDED_protocols_forge_components_BDR_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace components {


/// @brief forward declaration for BDR
class BDR;


/// @brief BDR owning pointer
typedef utility::pointer::shared_ptr< BDR > BDR_OP;


/// @brief BDR const owning pointer
typedef utility::pointer::shared_ptr< BDR const > BDR_COP;


/// @brief BDR access pointer
typedef utility::pointer::weak_ptr< BDR > BDR_AP;


/// @brief BDR const access pointer
typedef utility::pointer::weak_ptr< BDR const > BDR_CAP;


} // namespace components
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_BDR_FWD_HH */
