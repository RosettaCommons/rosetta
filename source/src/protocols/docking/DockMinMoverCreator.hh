// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/docking/DockMinMoverCreator.hh
/// @brief  Declaration of the MoverCreator class for the DockMinMover
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_docking_DockMinMoverCreator_HH
#define INCLUDED_protocols_docking_DockMinMoverCreator_HH

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace docking {

class DockMinMoverCreator : public moves::MoverCreator
{
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();

};

} // namespace docking
} // namespace protocols

#endif // INCLUDED_protocols_docking_DockMinMoverCreator_HH