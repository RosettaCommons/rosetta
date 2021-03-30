// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file trRosetta_protocols/movers/trRosettaProtocolMoverCreator.hh
/// @brief The full trRosetta structure prediction protocol from Yang et al, converted to C++ and implemented as a mover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMoverCreator_HH
#define INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace trRosetta_protocols {
namespace movers {

class trRosettaProtocolMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //movers
} //trRosetta_protocols
} //protocols

#endif //INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMoverCreator_HH
