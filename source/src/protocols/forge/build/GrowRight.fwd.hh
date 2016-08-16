// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/GrowRight.fwd.hh
/// @brief forward declaration for GrowRight
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_GrowRight_fwd_hh
#define INCLUDED_protocols_forge_build_GrowRight_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for GrowRight
class GrowRight;


/// @brief owning pointer for GrowRight
typedef utility::pointer::shared_ptr< GrowRight > GrowRightOP;


/// @brief owning pointer for const GrowRight
typedef utility::pointer::shared_ptr< GrowRight const > GrowRightCOP;


/// @brief access pointer for GrowRight
typedef utility::pointer::weak_ptr< GrowRight > GrowRightAP;


/// @brief access pointer for const GrowRight
typedef utility::pointer::weak_ptr< GrowRight const > GrowRightCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_GrowRight_FWD_HH */
