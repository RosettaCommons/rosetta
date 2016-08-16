// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/AntibodyEnumManager.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyEnumManager_fwd_hh
#define INCLUDED_protocols_antibody_AntibodyEnumManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace antibody {

// Forward
class AntibodyEnumManager;

typedef utility::pointer::shared_ptr< AntibodyEnumManager > AntibodyEnumManagerOP;
typedef utility::pointer::shared_ptr< AntibodyEnumManager const > AntibodyEnumManagerCOP;


} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_antibody_AntibodyEnumManager.fwd.hh
