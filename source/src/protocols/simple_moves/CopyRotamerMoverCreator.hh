// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CopyRotamerMoverCreator.hh
/// @brief A mover to copy a rotamer (residue identity and conformation) from one position in a pose to another.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

#ifndef INCLUDED_protocols_simple_moves_CopyRotamerMoverCreator_hh
#define INCLUDED_protocols_simple_moves_CopyRotamerMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace simple_moves {

class CopyRotamerMoverCreator : public protocols::moves::MoverCreator {

public:


	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} //protocols
} //simple_moves

#endif //INCLUDED_protocols/simple_moves_CopyRotamerMover_fwd_hh
