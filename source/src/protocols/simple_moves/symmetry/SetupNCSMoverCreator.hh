// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   SetupNCSMoverCreator.hh
/// @brief  protocols::moves::MoverCreator for SetupNCSMover
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_simple_moves_symmetry_SetupNCSMoverCreator_hh
#define INCLUDED_protocols_simple_moves_symmetry_SetupNCSMoverCreator_hh


// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace simple_moves {
namespace symmetry {

class SetupNCSMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static  std::string mover_name();
};

}
}
}

#endif

