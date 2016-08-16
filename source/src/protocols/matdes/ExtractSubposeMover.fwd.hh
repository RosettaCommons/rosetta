// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file src/protocols/matdes/ExtractSubposeMover.cc
/// @brief  Extract primary component associated with symdofs and all neighboring components.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_matdes_ExtractSubposeMover_fwd_HH
#define INCLUDED_protocols_matdes_ExtractSubposeMover_fwd_HH
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace matdes {

class ExtractSubposeMover;

typedef utility::pointer::shared_ptr< ExtractSubposeMover > ExtractSubposeMoverOP;
typedef utility::pointer::shared_ptr< ExtractSubposeMover const > ExtractSubposeMoverCOP;


} // matdes
} // protocols

#endif
