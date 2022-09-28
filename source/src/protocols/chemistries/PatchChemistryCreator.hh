// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/PatchChemistryCreator.hh
/// @brief  A chemistry which applies a patch
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_PatchChemistryCreator_HH
#define INCLUDED_protocols_chemistries_PatchChemistryCreator_HH

// Package headers
#include <protocols/chemistries/Chemistry.fwd.hh>
#include <protocols/chemistries/ChemistryCreator.hh>

namespace protocols {
namespace chemistries {

class PatchChemistryCreator : public ChemistryCreator {
public:
	ChemistryOP create_chemistry() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

} //namespace chemistries
} //namespace protocols

#endif
