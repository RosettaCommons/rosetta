// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/PoseOutputterFactory.hh
/// @brief  PoseOutputterFactory class that holds the list of classes able to create Poses from inputs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_pose_outputters_PoseOutputterFactory_hh
#define INCLUDED_protocols_jd3_pose_outputters_PoseOutputterFactory_hh

// Unit headers
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.fwd.hh>

// Package headers
#include <protocols/jd3/pose_outputters/PoseOutputter.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// c++ headers
#include <list>
#include <map>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief This templated class will register an instance of an
/// PoseOutputterCreator (class T) with the PoseOutputterFactory.  It will ensure
/// that no PoseOutputterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class PoseOutputterRegistrator : public utility::factory::WidgetRegistrator< PoseOutputterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PoseOutputterFactory, T > parent;
public:
	PoseOutputterRegistrator() : parent() {}
};



class PoseOutputterFactory : public utility::SingletonBase< PoseOutputterFactory >
{
public:
	friend class utility::SingletonBase< PoseOutputterFactory >;

	typedef std::map< std::string, PoseOutputterCreatorOP > PoseOutputterMap;
	typedef std::list< PoseOutputterCreatorCOP > CreatorList;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

	virtual ~PoseOutputterFactory();

	void factory_register( PoseOutputterCreatorOP creator );

	/// @brief Create a pose outputter given its identifying string, e.g. the name of an XML element.
	PoseOutputterOP new_pose_outputter( std::string const & ) const;

	/// @brief Construct a PoseOutputter for a job relying on the command line in the absence of a specified
	/// outputter class in the job-definition file.
	PoseOutputterOP
	pose_outputter_from_command_line() const;

	/// @brief Should the Factory throw an exception or call utility::exit when it encounters the
	/// second of two PoseOutputterCreators with the same keyname?  It's default behavior is to
	/// call utility::exit, but this method allows you to set it so that it will throw an
	/// exception instead (which is unit testable).
	void set_throw_on_double_registration();

	/// @brief The %ResidueSelectorFactory is the point of entry for the definition of the XML Schemas
	/// for every ResidueSelector that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "residue_selector" listing each of the residue-selector-complex types that may
	/// be initialized using the %ResidueSelectorFactory and to iterate across each of the
	/// ResidueSelectorCreators it contains asking them for the XML schema of the ResidueSelector they
	/// are responsible for creating.
	void define_pose_outputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	void list_options_read( utility::options::OptionKeyList & read_options ) const;

	static std::string pose_outputter_xml_schema_group_name();
	static std::string complex_type_name_for_pose_outputter( std::string const & outputter_key );

private:
	PoseOutputterFactory();

	// Unimplemented -- uncopyable
	PoseOutputterFactory( PoseOutputterFactory const & ) = delete;
	PoseOutputterFactory const & operator = ( PoseOutputterFactory const & ) = delete;

private:

	PoseOutputterMap pose_outputter_creator_map_;
	CreatorList creator_list_;
	bool throw_on_double_registration_;

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PoseOutputterFactory_HH
