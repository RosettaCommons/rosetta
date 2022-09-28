// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryFactory.hh
/// @brief  Class for instantiating arbitrary Chemistry objects from a string
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_ChemistryFactory_HH
#define INCLUDED_protocols_chemistries_ChemistryFactory_HH

// Package headers
#include <protocols/chemistries/Chemistry.fwd.hh>
#include <protocols/chemistries/ChemistryCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/SingletonBase.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>
#include <string>

namespace protocols {
namespace chemistries {

class ChemistryFactory : public utility::SingletonBase< ChemistryFactory > {
public:
	friend class utility::SingletonBase< ChemistryFactory >;

	void factory_register( ChemistryCreatorOP creator );
	bool has_type( std::string const & ) const;

	ChemistryOP new_chemistry(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;

	void set_throw_on_double_registration();

	void define_chemistry_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	static std::string chemistry_xml_schema_group_name();

private:
	ChemistryFactory();

	// Unimplemented -- uncopyable
	ChemistryFactory( ChemistryFactory const & ) = delete;
	ChemistryFactory const & operator = ( ChemistryFactory const & ) = delete;

private:
	std::map< std::string, ChemistryCreatorOP > creator_map_;
	bool throw_on_double_registration_;
};


} //namespace chemistries
} //namespace protocols


#endif
