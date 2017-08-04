// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MetalContactsConstraintGeneratorCreator.hh
/// @brief This constraint generator takes residue selectors for a residue/residues containing metal ion(s) and for residue(s) for which to set up contacts. It allows users to specify which base atoms will be used to define angles/dihedrals to constrain; ideal values for angles/dihedrals/distances; and an option to constrain to native values.
/// @author guffysl (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_constraint_generator_MetalContactsConstraintGeneratorCreator_hh
#define INCLUDED_protocols_constraint_generator_MetalContactsConstraintGeneratorCreator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>

namespace protocols {
namespace constraint_generator {

class MetalContactsConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	protocols::constraint_generator::ConstraintGeneratorOP
	create_constraint_generator() const override;

	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols/constraint_generator_MetalContactsConstraintGenerator_fwd_hh
