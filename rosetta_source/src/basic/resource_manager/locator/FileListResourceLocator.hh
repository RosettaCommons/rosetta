// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/FileSystemResourceLocater.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_locator_FileListResourceLocator_hh
#define INCLUDED_basic_resource_manager_locator_FileListResourceLocator_hh


#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/locator/FileListResourceLocator.fwd.hh>


namespace basic {
namespace resource_manager {
namespace locator {

class FileListResourceLocator : public basic::resource_manager::ResourceLocator
{
public:
	FileListResourceLocator();

	FileListResourceLocator(
		FileListResourceLocator const & src);

	virtual ~FileListResourceLocator();

	virtual
	void
	show(
		std::ostream & out) const;

	//friend
	//std::ostream &
	//operator<<(
	//	std::ostream & out,
	//	const FileSystemResourceLocator & file_system_resource_locator);

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

private:
	std::ios_base::openmode open_mode_;

};

}
}
}

#endif // INCLUDED_basic_resource_manager_locator_FileListResourceLocator_hh
