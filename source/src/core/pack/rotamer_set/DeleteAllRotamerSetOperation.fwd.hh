// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /rosetta/rosetta_source/src/protocols/sewing/DeleteAllRotamerSetOperation.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_core_pack_rotamer_set_DeleteAllRotamerSetOperation_FWD_HH
#define INCLUDED_core_pack_rotamer_set_DeleteAllRotamerSetOperation_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace rotamer_set {

class DeleteAllRotamerSetOperation;

typedef utility::pointer::shared_ptr< DeleteAllRotamerSetOperation > DeleteAllRotamerSetOperationOP;
typedef utility::pointer::shared_ptr< DeleteAllRotamerSetOperation const > DeleteAllRotamerSetOperationCOP;


} //namespace rotamer_set
} //namespace pack
} //namespace core

#endif
