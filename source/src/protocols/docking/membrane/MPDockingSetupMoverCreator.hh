// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/MPDockingSetupMoverCreator.hh
/// @brief      Setup MPDock protocol
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @author     Rebecca Alford (rfalford12@gmail.com) - Added RosettaScripts hook
/// @note       Last Modified (2/9/15)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingSetupMoverCreator_hh
#define INCLUDED_protocols_docking_membrane_MPDockingSetupMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace docking {
namespace membrane {

/// @brief Mover Creator
class MPDockingSetupMoverCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingSetupMoverCreator_hh
