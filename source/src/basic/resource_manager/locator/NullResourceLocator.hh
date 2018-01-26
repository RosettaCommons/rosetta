// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/NullResourceLocater.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_locator_NullResourceLocater_hh
#define INCLUDED_basic_resource_manager_locator_NullResourceLocater_hh

//unit headers
#include <basic/resource_manager/locator/NullResourceLocator.fwd.hh>

// Package headers
#include <basic/resource_manager/ResourceLocator.hh>

//utility headers
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {
namespace locator {

/// @brief The %NullStream acts as an empty stream object that may be
/// returned by the NullResourceLocator.  It does not open any files.
class NullStream : public basic::resource_manager::ResourceStream
{
public:
	NullStream();

public:

	virtual
	~NullStream();

	virtual
	std::istream &
	stream();

private: // members
	utility::io::izstream stream_;

};


/// @brief The %NullResourceLocator is meant for cases where a resource can
/// be created without reading from an input file.  It goes through the motions
/// of returning a ResourceStream (an empty NullStream) as is required of all
/// ResourceLocators, but the stream that it creates will not be used.
class NullResourceLocator : public basic::resource_manager::ResourceLocator
{
public:
	NullResourceLocator();

	NullResourceLocator(NullResourceLocator const & src);

	virtual ~NullResourceLocator();

	/// @brief Describe this instance to the given output stream; since there is no
	/// data in this class, merely print the name of this class.
	void
	show(
		std::ostream & out) const override;

	/// @brief Return the name of this class: "NullResourceLocator"
	std::string
	type() const override;

	/// @brief Create an empty NullResource object that will not be used
	ResourceStreamOP
	locate_resource_stream(
		std::string const & input_id
	) const override;

	/// @brief Noop, since there is no data in this class.
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	) override;

	/// @brief Describe the schema for this resource locator to the XSD.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static std::string classname();
};

} // namespace locator
} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_locator_NullResourceLocator_hh
