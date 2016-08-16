// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/movers/DockSetupMover.hh
/// @brief  Declaration of the MoverCreator class for the DockAndRetrieveSidechains
/// @author Oliver Lange (oliver.lange@tum.de)

#ifndef INCLUDED_protocols_docking_DockSetupMoverCreator_hh
#define INCLUDED_protocols_docking_DockSetupMoverCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace docking {

class DockSetupMoverCreator : public moves::MoverCreator
{
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static  std::string mover_name();

};

}
}

#endif
