// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/metal_interface/RemoveMetalConnectionsMoverCreator.hh
/// @brief Creator for a mover that removes the connections to metals that were added by the SetupMetalsMover
/// or by the -auto_setup_metals flag.
/// @details This mover:
///     - Removes the bonds between metals and metal-binding residues.
///     - Reverts metal-liganding residues back to their pre-bonded types.
///     - Reverts metals back to their pre-bonded types.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMoverCreator_HH
#define INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace metal_interface {

class RemoveMetalConnectionsMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //metal_interface
} //protocols

#endif //INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMoverCreator_HH
