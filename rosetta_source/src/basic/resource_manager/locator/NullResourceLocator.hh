// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/NullResourceLocater.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_locator_NullResourceLocater_hh
#define INCLUDED_basic_resource_manager_locator_NullResourceLocater_hh

//unit headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/locator/NullResourceLocator.fwd.hh>


//project headers
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {
namespace locator {

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



class NullResourceLocator : public basic::resource_manager::ResourceLocator
{
public:
	NullResourceLocator();

	// Undefined, commenting out to fix PyRosetta build  NullResourceLocator(NullResourceLocator const & src);

	virtual ~NullResourceLocator();

	virtual
	void
	show(
		std::ostream & out) const;

	virtual
	std::string
	type() const;

	/// @brief Create a ResourceStream object from the given resource
	/// source, so that its stream can be passed to the ResourceLoader
	virtual
	ResourceStreamOP
	locate_resource_stream(
		std::string const & locator_tag
	) const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr tag
	);

};

} // namespace locator
} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_locator_NullResourceLocator_hh
