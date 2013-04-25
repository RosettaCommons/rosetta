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

#ifndef INCLUDED_basic_resource_manager_locator_FileSystemResourceLocater_hh
#define INCLUDED_basic_resource_manager_locator_FileSystemResourceLocater_hh

//unit headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/locator/FileSystemResourceLocator.fwd.hh>


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

class FileStream : public basic::resource_manager::ResourceStream
{
public:
	FileStream();

	FileStream(
		std::string const & filename,
		std::ios_base::openmode open_mode = std::ios_base::in
	);

private:
	///@brief This is private. The FileStream shouldn't be copied
	FileStream(
		FileStream const & src);

public:

	virtual
	~FileStream();

	void
	open(
		std::string const & filename,
		std::ios_base::openmode open_mode = std::ios_base::in
	);


	virtual
	std::istream &
	stream();

private: // members
	utility::io::izstream stream_;

};



class FileSystemResourceLocator : public basic::resource_manager::ResourceLocator
{
public:
	FileSystemResourceLocator(
		std::ios_base::openmode open_mode = std::ios_base::in);

	FileSystemResourceLocator(
		FileSystemResourceLocator const & src);

	virtual ~FileSystemResourceLocator();

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


	void
	set_open_mode(
		std::ios_base::openmode open_mode);

	std::ios_base::openmode
	get_open_mode() const;

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

} // namespace locator
} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_locator_FileSystemResourceLocator_hh
