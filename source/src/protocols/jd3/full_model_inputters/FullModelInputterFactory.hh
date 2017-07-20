// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model_inputters/FullModelInputterFactory.hh
/// @brief  FullModelInputterFactory class that holds the list of classes able to create Poses from inputs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_full_model_inputters_FullModelInputterFactory_hh
#define INCLUDED_protocols_jd3_full_model_inputters_FullModelInputterFactory_hh

// Unit headers
#include <protocols/jd3/full_model_inputters/FullModelInputterFactory.fwd.hh>

// Package headers
#include <protocols/jd3/full_model_inputters/FullModelInputter.fwd.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputterCreator.fwd.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputSource.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/SingletonBase.hh>

// c++ headers
#include <list>
#include <map>

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

/// @brief This templated class will register an instance of an
/// FullModelInputterCreator (class T) with the FullModelInputterFactory.  It will ensure
/// that no FullModelInputterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class FullModelInputterRegistrator : public utility::factory::WidgetRegistrator< FullModelInputterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< FullModelInputterFactory, T > parent;
public:
	FullModelInputterRegistrator() : parent() {}
};



class FullModelInputterFactory : public utility::SingletonBase< FullModelInputterFactory > {
public:

	typedef std::map< std::string, FullModelInputterCreatorOP > FullModelInputterMap;
	typedef std::list< FullModelInputterCreatorCOP > CreatorList;
	typedef utility::vector1< std::pair< FullModelInputSourceOP, FullModelInputterOP > > FullModelInputSourcesAndInputters;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

	friend class utility::SingletonBase< FullModelInputterFactory >;

	void factory_register( FullModelInputterCreatorOP creator );

	/// @brief Create a pose inputter given its identifying string, e.g. the name of an XML element.
	FullModelInputterOP new_full_model_inputter( std::string const & ) const;

	/// @brief Read the command line (i.e. in the absence of a job-definition file) to determine the set
	/// of input poses.  Note that this reads from the global options system and does not read per-job
	/// options -- this is necessary because there isn't a job-definition file from which to read per-job
	/// options. Note also that this will read all options on the command line, possibly from different
	/// FullModelInputters.  The order in which the FullModelInputSources are reported are sorted by alphabetically
	/// by the keys of the FullModelInputters (e.g. "PDB" will precede "Silent").

	FullModelInputSourcesAndInputters
	full_model_inputs_from_command_line() const;

	/// @brief Should the Factory throw an exception or call utility::exit when it encounters the
	/// second of two FullModelInputterCreators with the same keyname?  It's default behavior is to
	/// call utility::exit, but this method allows you to set it so that it will throw an
	/// exception instead (which is unit testable).
	void set_throw_on_double_registration();

	/// @brief The %ResidueSelectorFactory is the point of entry for the definition of the XML Schemas
	/// for every ResidueSelector that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "residue_selector" listing each of the residue-selector-complex types that may
	/// be initialized using the %ResidueSelectorFactory and to iterate across each of the
	/// ResidueSelectorCreators it contains asking them for the XML schema of the ResidueSelector they
	/// are responsible for creating.
	void define_full_model_inputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	void list_options_read( utility::options::OptionKeyList & read_options ) const;

	CreatorList const & full_model_inputter_creators() const;

	static std::string full_model_inputter_xml_schema_group_name();
	static std::string complex_type_name_for_full_model_inputter( std::string const & inputter_key );

private:
	FullModelInputterFactory();

	// Unimplemented -- uncopyable
	FullModelInputterFactory( FullModelInputterFactory const & ) = delete;
	FullModelInputterFactory const & operator = ( FullModelInputterFactory const & ) = delete;

private:

	FullModelInputterMap full_model_inputter_creator_map_;
	CreatorList creator_list_;
	bool throw_on_double_registration_;

};

} // namespace full_model_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_FullModelInputterFactory_HH
