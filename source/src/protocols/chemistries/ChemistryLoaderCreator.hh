// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryLoaderCreator.hh
/// @brief  Creator classes for the ChemistryLoader
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_ChemistryLoaderCreator_hh
#define INCLUDED_protocols_chemistries_ChemistryLoaderCreator_hh

// Package headers
#include <protocols/parser/DataLoaderCreator.hh>

// Utility Headers

namespace protocols {
namespace chemistries {

class ChemistryLoaderCreator : public parser::DataLoaderCreator
{
public:
	parser::DataLoaderOP create_loader() const override;
	std::string keyname() const override;
	DerivedNameFunction schema_ct_naming_function() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace chemistries
} //namespace protocols

#endif
