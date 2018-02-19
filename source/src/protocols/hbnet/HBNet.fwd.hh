// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/hbnet/HBNet.fwd.hh
/// @brief
/// @author Scott Boyken (sboyken@gmail.com)

#ifndef INCLUDED_protocols_hbnet_HBNet_fwd_hh
#define INCLUDED_protocols_hbnet_HBNet_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace hbnet {

class HBNet;
using HBNetOP = utility::pointer::shared_ptr< HBNet >;
using HBNetCOP = utility::pointer::shared_ptr< HBNet const >;

struct HBondResStruct;
using HBondResStructOP = utility::pointer::shared_ptr< HBondResStruct >;
using HBondResStructCOP = utility::pointer::shared_ptr< HBondResStruct const >;

struct HBondNetStruct;
using HBondNetStructOP = utility::pointer::shared_ptr< HBondNetStruct >;
using HBondNetStructCOP = utility::pointer::shared_ptr< HBondNetStruct const >;

}  //namespace hbnet
}  //namespace protocols

#endif
