// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/AnchoredDesign/AnchorMoversData.fwd.hh
/// @brief  Anchor classes forward declarations header
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_anchored_design_AnchorMoversData_fwd_hh
#define INCLUDED_protocols_anchored_design_AnchorMoversData_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace anchored_design {

//Forwards and OP typedefs
class AnchorMoversData;
typedef utility::pointer::shared_ptr< AnchorMoversData > AnchorMoversDataOP;
typedef utility::pointer::shared_ptr< AnchorMoversData const > AnchorMoversDataCOP;

}//AnchoredDesign
}//protocols

#endif //INCLUDED_protocols_AnchoredDesign_AnchorMoversData_FWD_HH

