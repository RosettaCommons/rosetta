// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/HBondGraph.fwd.hh
/// @brief HBondGraph (derived from core::graph::Graph) forward declarations
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_protocols_hbnet_HBondGraph_FWD_HH
#define INCLUDED_protocols_hbnet_HBondGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace hbnet {

class HBondGraph;
typedef utility::pointer::shared_ptr< HBondGraph > HBondGraphOP;
typedef utility::pointer::shared_ptr< HBondGraph const > HBondGraphCOP;

} //hbnet
} //core

#endif
