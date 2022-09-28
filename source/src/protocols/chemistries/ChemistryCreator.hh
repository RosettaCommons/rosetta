// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryCreator.hh
/// @brief  Class for instantiating a particular Chemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_ChemistryCreator_HH
#define INCLUDED_protocols_chemistries_ChemistryCreator_HH

// Package headers
#include <protocols/chemistries/Chemistry.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace chemistries {

class ChemistryCreator : public utility::VirtualBase {
public:
	virtual ChemistryOP create_chemistry() const = 0;
	virtual std::string keyname() const = 0;

	/// @brief Define the structure of the XML file for the Chemistry that this
	/// ChemistryCreator instantiates using the XML Schema language.
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const = 0;
};


} //namespace chemistries
} //namespace protocols


#endif
