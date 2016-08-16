// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/enzdes/EnzdesBaseProtocol.fwd.hh
///
/// @brief
/// @author Florian Richter

#ifndef INCLUDED_protocols_enzdes_EnzdesBaseProtocol_fwd_hh
#define INCLUDED_protocols_enzdes_EnzdesBaseProtocol_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace enzdes {

class EnzdesBaseProtocol;
typedef utility::pointer::shared_ptr< EnzdesBaseProtocol > EnzdesBaseProtocolOP;
typedef utility::pointer::shared_ptr< EnzdesBaseProtocol const > EnzdesBaseProtocolCOP;

typedef utility::pointer::weak_ptr< EnzdesBaseProtocol > EnzdesBaseProtocolAP;
typedef utility::pointer::weak_ptr< EnzdesBaseProtocol const > EnzdesBaseProtocolCAP;


}  //namespace enzdes
}  //namespace protocols

#endif
