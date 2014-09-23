// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/AntibodyInfo.fwd.hh
/// @brief  AntibodyInfo class forward declarations header
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyInfo_fwd_hh
#define INCLUDED_protocols_antibody_AntibodyInfo_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

// C++ Headers
namespace protocols {
namespace antibody {

// Forward
class AntibodyInfo;

typedef utility::pointer::shared_ptr< AntibodyInfo > AntibodyInfoOP;
typedef utility::pointer::shared_ptr< AntibodyInfo const > AntibodyInfoCOP;
typedef  utility::pointer::weak_ptr< AntibodyInfo >  AntibodyInfoAP;
typedef  utility::pointer::weak_ptr< AntibodyInfo const >  AntibodyInfoCAP;



} //namespace antibody
} //namespace protocols

#endif //INCLUDED_protocols_AntibodyInfo_FWD_HH
