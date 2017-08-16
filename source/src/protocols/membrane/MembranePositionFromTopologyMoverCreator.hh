// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/membrane/MembranePositionFromTopologyMoverCreator.hh
///
/// @brief      Compute the Initial Position of the membrane
/// @details Compute the initial position of the membrane from
///    sequence or structure (can be specified by the user at construction
///    or as a seutp cmd flag)
/// @file  protocols/membrane/MembranePositionRTMover.hh
///
/// @brief  Sets the membrane position normal and center
/// @details Sets the membrane position normal and center
///    CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///    Last Modified: 6/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)
///
///    NOTE: Requires a membrane pose!
///    NOTE: sequence not yet implemented
///    Last Modified: 6/21/14
///
/// @author  Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionFromTopologyMoverCreator_hh
#define INCLUDED_protocols_membrane_MembranePositionFromTopologyMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator
class MembranePositionFromTopologyMoverCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionFromTopologyMoverCreator_hh
