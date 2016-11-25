// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/rdf/RDFFunctionCreator.fwd.hh
/// @brief  Header for base class for RDFFunctionCreator
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_ligand_docking_rdf_RDFFunctionCreator_hh
#define INCLUDED_protocols_ligand_docking_rdf_RDFFunctionCreator_hh

// Unit Headers
#include <protocols/ligand_docking/rdf/RDFFunctionCreator.fwd.hh>

// Package Headers
#include <protocols/ligand_docking/rdf/RDFBase.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// c++ headers
#include <string>

namespace protocols {
namespace ligand_docking {
namespace rdf {

/// @brief The Creator class is responsible for creating a particular
/// RDF Function class.
class RDFFunctionCreator : public utility::pointer::ReferenceCount
{
public:
	RDFFunctionCreator() {}
	virtual ~RDFFunctionCreator() {}

	virtual RDFBaseOP create_rdf_function() const = 0;
	virtual std::string type_name() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

}
} //namespace
} //namespace

#endif
