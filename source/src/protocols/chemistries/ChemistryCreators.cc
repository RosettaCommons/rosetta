// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryCreators.hh
/// @brief  Class declarations for the ChemistryCreators for Chemistry objects in core/chemical/modifications
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/chemistries/ChemistryCreators.hh>
#include <protocols/chemistries/Chemistry.hh>

#include <protocols/chemistries/WrappedChemistries.hh>

namespace protocols {
namespace chemistries {

ChemistryOP
ReprotonateCreator::create_chemistry() const {
	return ChemistryOP( new ReprotonateChemistry );
}

std::string
ReprotonateCreator::keyname() const {
	return ReprotonateChemistry::class_name();
}

void
ReprotonateCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ReprotonateChemistry::provide_xml_schema( xsd );
}

} //namespace chemistries
} //namespace protocols

