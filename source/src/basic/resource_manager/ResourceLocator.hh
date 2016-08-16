// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_magnager/ResourceLocater.hh
/// @brief  Locate a resource in a data store and make it available to a ResourceLoader in the form of an ResourceStream
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLocater_hh
#define INCLUDED_basic_resource_manager_ResourceLocater_hh

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//unit headers
#include <basic/resource_manager/ResourceLocator.fwd.hh>

//project headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {

/// @brief The %ResourceStream represents an abstract class for packaging
/// up a standard istream so that data that the ResourceLocator needs to
/// deliver to a ResourceLoader can come from arbitrary sources (e.g.
/// from either a file or from a database).
class ResourceStream : public utility::pointer::ReferenceCount
{
public:
	virtual
	~ResourceStream();

	/// @brief Return an istream reference so that the ResourceLoader can access arbitrary data
	/// returned by the %ResourceLocator
	virtual
	std::istream &
	stream() = 0;

};

/// @brief %ResourceLocator classes are responsible for retrieving data
/// from a data store that will be used to construct a Resource.  This
/// data store could be a file system or a database or any other
/// place where data is stored.
///
/// The ResourceManager asks the ResourceLocator to produce a
/// ResourceStream object when given a "locator_id."  A "locator_id"
/// is what's needed to identify a data source from a data store:
/// for example, for the FileSystemResourceLocator, the locator id
/// is a file name; for the DatabaseResourceLocator, it would be
/// a database query.
class ResourceLocator : public utility::pointer::ReferenceCount
{
public:
	/// @brief Construct a %ResourceLocator and initialize its name
	/// (its locator_tag) to the empty string.
	ResourceLocator();

	/// @brief Construct a %ResourceLocator while setting its name to the input
	/// locater_tag
	ResourceLocator(
		std::string const & locator_tag);

	/// @brief Copy construct a %ResourceLocator from an example locator
	ResourceLocator( ResourceLocator const & src );

	virtual ~ResourceLocator();


	/// @brief Create and return a ResourceStream object from the given locator_id
	/// so that its stream can be passed to the ResourceLoader
	virtual
	ResourceStreamOP
	locate_resource_stream(
		std::string const & locator_id
	) const = 0;

	/// @brief Initialize the parameters for this %ResourceLocator from the contents of
	/// an XML file.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	) = 0;

	/// @brief Set the name for this %ResourceLocator
	virtual
	void
	locator_tag( std::string const & locator_tag );

	/// @brief Return the name for this %ResourceLocator
	virtual
	std::string
	locator_tag() const;

	/// @brief Write a description of this %ResourceLocator to an out stream
	virtual
	void
	show(
		std::ostream & out
	) const = 0;

	/// @brief Return the class name for this %ResourceLocator instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const = 0;

private:
	std::string locator_tag_;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLoader_hh
