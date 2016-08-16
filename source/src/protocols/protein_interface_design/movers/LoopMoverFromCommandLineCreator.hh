// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/LoopMoverFromCommandLineCreator.hh
/// @brief  Declaration of the MoverCreator class for the LoopRemodelFromCommandLine
/// @author Jordan Willis(jordan.r.willis@vanderbilt.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_LoopMoverFromCommandLineCreator_HH
#define INCLUDED_protocols_protein_interface_design_movers_LoopMoverFromCommandLineCreator_HH

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {
class LoopMoverFromCommandLineCreator : public moves::MoverCreator
{
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static  std::string mover_name();

};

}
}
}

#endif
