// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/FileContentsMap.hh
/// @brief  Declaration for a file-driven definition for the active states in a multistate design
///         and the set of mathematical expressions that together define the fitness function
///         for the sequence being designed.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_io_FileContentsMap_hh
#define INCLUDED_utility_io_FileContentsMap_hh

// Unit headers
#include <utility/io/FileContentsMap.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>
#include <string>

// Platform headers
#include <platform/types.hh>

namespace utility {
namespace io {

/// @brief The %FileContentsMap is a class that will minimize the number of times
/// files are accessed from disk.  The first time the contents of a file are
/// requested, it will cache the contents the file in memory as a strings.  All
/// subsequent requests for the contents of that file will then be delivered
/// without having to go to disk.  WARNING: The %FileContentsMap will not update
/// the contents that it holds for a file after the first time it loads it, so if
/// the contents change on disk, the %FileContentsMap will be out of date.
///
/// The %FileContentsMap offers two settings that are both off by default.
///
/// 1) You can set whether it should delete files when they reach a "read limit,"
/// which in turn requires that you first define a read limit for the file.
/// This can prevent the %FileContentsMap from loading too much into memory as
/// long as the utility of a file is expected to be short (i.e. it's read several
/// times in a short period and then never again).  Setting a file's read limit to
/// zero signals that the contents of the file should never be deleted from memory
/// i.e. that it has an infinite read limit. If the contents of a file is requested
/// before read limit has been set, then the read limit is set to zero
/// (i.e. infinity).
///
/// 2) You can set whether it should refuse to load unexpected files.  Before
/// you request the contents of a file, you must provide a read limit (zero is
/// an acceptible value here; see above) for that file, or the %FileContentsMap
/// will throw an exception when you request that file.
class FileContentsMap : public utility::pointer::ReferenceCount
{
public:

	/// @brief Construct an empty %FileContentsMap with both boolean options set to false
	FileContentsMap();

	/// @brief Construct and populate a %FileContentsMap from an input std::map between
	/// file names and file contents.  Both boolean options are set to false.
	FileContentsMap( std::map< std::string, std::string > const & fcontents );

	/// @brief Enable or disable the behavior that the %FileContentsMap will delete
	/// the contents of files from memory after they have been read as many times
	/// as they were expected to be read (their read limit).
	void delete_contents_at_nread_limit( bool setting );

	/// @brief Return the boolean value describing whether the %FileContentsMap will
	/// delete the contents of files from memory after they have reached their
	/// read limit.
	bool delete_contents_at_nread_limit() const;

	/// @brief Enable or disable the behavior that the FCM will refuse to load a file if
	/// a non-zero expected number of reads for that file has not been previously provided.
	void refuse_unexpected_files( bool setting );

	/// @brief Return the boolean value describing whether the %FileContentsMap will refuse
	/// to load files that it does not expect to load.
	bool refuse_unexpected_files() const;

	/// @brief Set the maximum number of reads for a file; if this is zero, then that
	/// file can be read an indefinite number of times.
	void set_nread_limit_for_file( std::string const & fname, platform::Size limit );

	/// @brief Increment the internal nread limit for a file
	void increment_nread_limit( std::string const & fname );

	/// @brief Return the number of times a particular file has been read
	platform::Size nreads_for_file( std::string const & fname ) const;

	/// @brief Return the read limit that has been set for this file; this will return
	/// a value of zero even if the file has not had a read limit set.  Use
	/// has_read_limit to determine if a read limit has been set for a particular file.
	platform::Size nread_limit_for_file( std::string const & fname ) const;

	/// @brief Return whether a read limit (even one of zero) has been set for a file
	bool has_read_limit_for_file( std::string const & fname ) const;

	/// @brief Return whether or not the %FileContentsMap has the contents of a file in memory
	bool has_file_contents( std::string const & fname ) const;

	/// @brief Set the contents of a file from an input string thereby avoiding a read from
	/// disk.
	void set_file_contents( std::string const & fname, std::string const & contents );

	/// @brief Delete the contents of a particular file from memory
	void delete_file_contents( std::string const & fname );

	/// @brief Delete the contents of all files and file-read counts from memory
	void clear();

	/// @brief Return the string contents of a file, possibly loading the contents of that file
	/// into memory for the first time, always updating internal data, and possibly
	/// deleting the contents of a file if it is expected to no longer be needed.
	std::string
	get_file_contents( std::string const & fname );

	/// @brief Return a reference to the string contents of a file, possibly loading the contents
	/// of that file into memory for the first time, but not deleting the contents of that file.
	std::string const &
	get_file_contents_ref( std::string const & fname );

private:
	bool delete_contents_at_nread_limit_;
	bool refuse_unexpected_files_;
	std::map< std::string, platform::Size > read_limit_;
	std::map< std::string, platform::Size > read_counts_;
	std::map< std::string, std::string > file_contents_;


};


}
}

#endif
