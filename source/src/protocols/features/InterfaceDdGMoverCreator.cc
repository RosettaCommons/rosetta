// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/InterfaceDdGMoverCreator.cc
/// @brief Calculates ddG of binding
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_InterfaceDdGMoverCreator_cc
#define INCLUDED_protocols_features_InterfaceDdGMoverCreator_cc

#include <protocols/moves/MoverCreator.hh>
#include <protocols/features/InterfaceDdGMover.hh>
#include <protocols/features/InterfaceDdGMoverCreator.hh>

namespace protocols {
namespace features {

/////////////// Creator ///////////////

protocols::moves::MoverOP
InterfaceDdGMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new InterfaceDdGMover );
}

std::string
InterfaceDdGMoverCreator::keyname() const
{
	return InterfaceDdGMoverCreator::mover_name();
}

std::string
InterfaceDdGMoverCreator::mover_name()
{
	return InterfaceDdGMover::mover_name();
}

void
InterfaceDdGMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceDdGMover::provide_xml_schema( xsd );
}

} //protocols
} //features

#endif //INCLUDED_protocols_features_InterfaceDdGMoverCreator_cc
