// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
///
/// @file   protocols/cyclic_peptide/CreateTorsionConstraint.hh
/// @brief  Add torsion constraints to the current pose conformation.  This class will create instances of mover CreateTorsionConstraint for the MoverFactory.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraintCreator_hh
#define INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraintCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace cyclic_peptide {

class CreateTorsionConstraintCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}

#endif //INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraintCreator_h

