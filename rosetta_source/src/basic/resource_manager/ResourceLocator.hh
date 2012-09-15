// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <utility/tag/Tag.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {

class ResourceStream : public utility::pointer::ReferenceCount
{
public:
	virtual
	~ResourceStream();

	virtual
	std::istream &
	stream() = 0;
};

class ResourceLocator : public utility::pointer::ReferenceCount
{
public:
	ResourceLocator();

	ResourceLocator(
		std::string const & locator_tag);

	ResourceLocator( ResourceLocator const & src );

	virtual ~ResourceLocator();


	/// @brief Create a ResourceStream object from the given resource
	/// source, so that its stream can be passed to the ResourceLoader
	virtual
	ResourceStreamOP
	locate_resource_stream(
		std::string const & input_tag
	) const = 0;

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr tag
	) = 0;

	virtual
	void
	locator_tag( std::string const & locator_tag );

	virtual
	std::string
	locator_tag() const;

	virtual
	void
	show(
		std::ostream & out
	) const = 0;

	/// @brief The class name for a particular ResourceLocator instance.
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
