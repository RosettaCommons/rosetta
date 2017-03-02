// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpSelectorFactory.hh
/// @brief  Class for instantiating arbitrary JumpSelectors from a string --> JumpSelectorCreator map
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_JumpSelectorFactory_HH
#define INCLUDED_core_select_jump_selector_JumpSelectorFactory_HH

// Package headers
#include <core/select/jump_selector/JumpSelector.fwd.hh>
#include <core/select/jump_selector/JumpSelectorCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>
#include <string>

namespace core {
namespace select {
namespace jump_selector {

class JumpSelectorFactory : public utility::SingletonBase< JumpSelectorFactory > {
public:
	friend class utility::SingletonBase< JumpSelectorFactory >;

	void factory_register( JumpSelectorCreatorOP creator );
	bool has_type( std::string const & ) const;

	/// @brief Get the XML schema for a given jump selector.
	/// @details Throws an error if the jump selector is unknown to Rosetta.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void provide_xml_schema( std::string const &selector_name, utility::tag::XMLSchemaDefinition & xsd ) const;

	JumpSelectorOP new_jump_selector(
		std::string const & selector_name,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;

	/// @brief Should the Factory throw an exception or call utility::exit when it encounters the
	/// second of two JumpSelctorCreators with the same keyname?  It's default behavior is to
	/// call utility::exit, but this method allows you to set it so that it will throw an
	/// exception instead (which is unit testable).
	void set_throw_on_double_registration();

	/// @brief The %JumpSelectorFactory is the point of entry for the definition of the XML Schemas
	/// for every JumpSelector that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "jump_selector" listing each of the jump-selector-complex types that may
	/// be initialized using the %JumpSelectorFactory and to iterate across each of the
	/// JumpSelectorCreators it contains asking them for the XML schema of the JumpSelector they
	/// are responsible for creating.
	void define_jump_selector_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Read access to the map of creator names to creators -- for unit testing purposes only
	std::map< std::string, JumpSelectorCreatorOP > const & creator_map() const;

	static std::string jump_selector_xml_schema_group_name();

private:
	JumpSelectorFactory();

	// Unimplemented -- uncopyable
	JumpSelectorFactory( JumpSelectorFactory const & ) = delete;
	JumpSelectorFactory const & operator = ( JumpSelectorFactory const & ) = delete;

private:
	std::map< std::string, JumpSelectorCreatorOP > creator_map_;
	bool throw_on_double_registration_;
};


} //namespace jump_selector
} //namespace select
} //namespace core


#endif
