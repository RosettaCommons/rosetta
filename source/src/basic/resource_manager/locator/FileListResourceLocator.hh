// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief The %FileListResourceLocator concatenates a set of listed files; e.g. useful
/// for constructing a pose from two separate PDB files.
class FileListResourceLocator : public basic::resource_manager::ResourceLocator
{
public:
	FileListResourceLocator();

	FileListResourceLocator(
		FileListResourceLocator const & src);

	virtual ~FileListResourceLocator();

	/// @brief Describe the %FileListResourceLocator to the output stringstream; since
	/// this class has no data, merely prints the name of the class.
	virtual
	void
	show(
		std::ostream & out
	) const;

	/// @brief Return the name for this class: "FileListResourceLocator"
	virtual
	std::string
	type() const;

	void
	set_open_mode(
		std::ios_base::openmode open_mode
	);

	std::ios_base::openmode
	get_open_mode() const;

	/// @brief Take the input locator tag and split it by whitespace, interpret each substring
	/// as the name of a file, and concatenate each file into a single stringstream to be
	/// returned.
	virtual
	ResourceStreamOP
	locate_resource_stream(
		std::string const & locator_tag
	) const;

	/// @brief Do nothing, since there is no data that the %FileListResourceLocator needs.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	);

private:
	std::ios_base::openmode open_mode_;

};

}
}
}

#endif // INCLUDED_basic_resource_manager_locator_FileListResourceLocator_hh
