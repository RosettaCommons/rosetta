// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddHelixSequenceConstraintsMoverCreator.hh
/// @brief This mover adds sequence constraints to the ends of each helix, requiring at least one positively-charged residue in the three C-terminal residues, and at least one negatively-charged resiude in the three N-terminal residues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMoverCreator_HH
#define INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace aa_composition {

class AddHelixSequenceConstraintsMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //protocols
} //aa_composition

#endif //INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMoverCreator_HH
