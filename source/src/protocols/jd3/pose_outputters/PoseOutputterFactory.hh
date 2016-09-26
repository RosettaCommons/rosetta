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
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.fwd.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputterCreator.fwd.hh>

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

/// @brief This templated class will register an instance of an
/// SecondaryPoseOutputterCreator (class T) with the PoseOutputterFactory.  It will ensure
/// that no SecondaryPoseOutputterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class SecondaryPoseOutputterRegistrator : public utility::factory::WidgetRegistrator< PoseOutputterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PoseOutputterFactory, T > parent;
public:
	SecondaryPoseOutputterRegistrator() : parent() {}
};



class PoseOutputterFactory : public utility::SingletonBase< PoseOutputterFactory >
{
public:
	typedef std::map< std::string, PoseOutputterCreatorOP > PoseOutputterMap;
	typedef std::map< std::string, SecondaryPoseOutputterCreatorOP > SecondaryOutputterMap;
	typedef std::list< PoseOutputterCreatorCOP > CreatorList;
	typedef std::list< SecondaryPoseOutputterCreatorCOP > SecondaryCreatorList;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

	friend class utility::SingletonBase< PoseOutputterFactory >;

	void factory_register( PoseOutputterCreatorOP creator );

	void factory_register( SecondaryPoseOutputterCreatorOP creator );

	/// @brief Create a pose outputter given its identifying string, e.g. the name of an XML element.
	PoseOutputterOP new_pose_outputter( std::string const & ) const;

	/// @brief Create a SecondaryPoseOutputter given its identifying string, e.g. the name of an XML element.
	SecondaryPoseOutputterOP new_secondary_outputter( std::string const & ) const;

	/// @brief Construct a PoseOutputter for a job relying on the command line in the absence of a specified
	/// outputter class in the job-definition file.
	PoseOutputterOP
	pose_outputter_from_command_line() const;

	/// @brief Construct all SecondaryPoseOutputters for a job relying on the command line in the absence of a specified
	/// list of outputter classes in the job-definition file.
	std::list< SecondaryPoseOutputterOP >
	secondary_pose_outputters_from_command_line() const;

	/// @brief Should the Factory throw an exception or call utility::exit when it encounters the
	/// second of two PoseOutputterCreators or two SecondaryPoseOutputters with the same keyname?
	/// It's default behavior is to call utility::exit, but this method allows you to set it so
	/// that it will throw an exception instead (which is unit testable).
	void set_throw_on_double_registration();

	/// @brief The %PoseOutputterFactory is the point of entry for the definition of the XML Schemas
	/// for every PoseOutputter that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "pose_outputter" listing each of the pose-outputter-complex types that may
	/// be initialized using the %PoseOutputterFactory and to iterate across each of the
	/// PoseOutputterCreators it contains asking them for the XML schema of the PoseOutputter they
	/// are responsible for creating.
	void define_pose_outputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief The %PoseOutputterFactory is the point of entry for the definition of the XML Schemas
	/// for every SecondaryPoseOutputter that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "secondary_pose_outputter" listing each of the pose-outputter-complex types that may
	/// be initialized using the %PoseOutputterFactory and to iterate across each of the SecondaryPoseOutputterCreators
	/// it contains asking them for the XML schema of the SecondaryPoseOutputter they are responsible for creating.
	void define_secondary_outputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	void list_options_read( utility::options::OptionKeyList & read_options ) const;

	static std::string pose_outputter_xml_schema_group_name();
	static std::string complex_type_name_for_pose_outputter( std::string const & outputter_key );

	static std::string secondary_pose_outputter_xml_schema_group_name();
	static std::string complex_type_name_for_secondary_pose_outputter( std::string const & outputter_key );

private:
	PoseOutputterFactory();

	// Unimplemented -- uncopyable
	PoseOutputterFactory( PoseOutputterFactory const & ) = delete;
	PoseOutputterFactory const & operator = ( PoseOutputterFactory const & ) = delete;

private:

	PoseOutputterMap pose_outputter_creator_map_;
	SecondaryOutputterMap secondary_pose_outputter_creator_map_;
	CreatorList creator_list_;
	SecondaryCreatorList secondary_creator_list_;
	bool throw_on_double_registration_;

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PoseOutputterFactory_HH
