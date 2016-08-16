// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/movers/TaskAwareCstsCreator.hh
/// @brief  Declaration of the MoverCreator class for the TaskAwareCsts
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_TaskAwareCstsCreator_hh
#define INCLUDED_protocols_protein_interface_design_movers_TaskAwareCstsCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

class TaskAwareCstsCreator : public moves::MoverCreator
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
