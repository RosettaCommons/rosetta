// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/FileContentsMap.cc
/// @brief  Implementation of the FileContentsMap class, which loads and then maintains the contents of
///         files as strings in memory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/io/FileContentsMap.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>


namespace utility {
namespace io {


FileContentsMap::FileContentsMap() :
	delete_contents_at_nread_limit_( false ),
	refuse_unexpected_files_( false )
{}

FileContentsMap::FileContentsMap(
	std::map< std::string, std::string > const & fcontents
) :
	delete_contents_at_nread_limit_( false ),
	refuse_unexpected_files_( false ),
	file_contents_( fcontents )
{}

void FileContentsMap::delete_contents_at_nread_limit( bool setting )
{
	delete_contents_at_nread_limit_ = setting;
}

bool FileContentsMap::delete_contents_at_nread_limit() const
{
	return delete_contents_at_nread_limit_;
}

void FileContentsMap::refuse_unexpected_files( bool setting )
{
	refuse_unexpected_files_ = setting;
}

bool FileContentsMap::refuse_unexpected_files() const
{
	return refuse_unexpected_files_;
}


void FileContentsMap::set_nread_limit_for_file( std::string const & fname, platform::Size limit )
{
	read_limit_[ fname ] = limit;
}

void FileContentsMap::increment_nread_limit( std::string const & fname )
{
	std::map< std::string, platform::Size >::iterator iter = read_limit_.find( fname );
	if ( iter == read_limit_.end() ) {
		read_limit_[ fname ] = 1;
		read_counts_[ fname ] = 0;
	}
	++(iter->second);
}

platform::Size
FileContentsMap::nreads_for_file( std::string const & fname ) const
{
	std::map< std::string, platform::Size >::const_iterator iter = read_counts_.find( fname );
	return iter == read_counts_.end() ? 0 : iter->second;
}

platform::Size
FileContentsMap::nread_limit_for_file( std::string const & fname ) const
{
	std::map< std::string, platform::Size >::const_iterator iter = read_limit_.find( fname );
	return iter == read_limit_.end() ? 0 : iter->second;
}

bool
FileContentsMap::has_read_limit_for_file( std::string const & fname ) const
{
	std::map< std::string, platform::Size >::const_iterator iter = read_limit_.find( fname );
	return iter != read_limit_.end();
}

bool
FileContentsMap::has_file_contents( std::string const & fname ) const
{
	std::map< std::string, std::string >::const_iterator iter = file_contents_.find( fname );
	return iter != file_contents_.end();
}


void FileContentsMap::set_file_contents( std::string const & fname, std::string const & contents )
{
	file_contents_[ fname ] = contents;
}

void FileContentsMap::delete_file_contents( std::string const & fname )
{
	std::map< std::string, std::string >::iterator iter = file_contents_.find( fname );
	if ( iter == file_contents_.end() ) return;
	file_contents_.erase( iter );
}

void FileContentsMap::clear()
{
	read_limit_.clear();
	read_counts_.clear();
	file_contents_.clear();
}

/// @throws Relies on utility::file_contents which throws a utility::excn::EXCN_Msg_Exception
/// if the desired file cannot be found.  It will also throw an exception if the %FileContentsMap
/// has the "refuse unexpected files" setting as true and the file that's requested does not
/// have a read limit that's been set.
std::string
FileContentsMap::get_file_contents( std::string const & filename )
{
	std::map< std::string, platform::Size >::const_iterator rl_iter = read_limit_.find( filename );
	std::map< std::string, platform::Size >::iterator rc_iter = read_counts_.find( filename );
	std::map< std::string, std::string >::iterator fc_iter = file_contents_.find( filename );

	if ( refuse_unexpected_files_ ) {
		if ( rl_iter == read_limit_.end() ) {
			throw utility::excn::EXCN_Msg_Exception( "Unexpected file-read requested: " + filename );
		}
		if ( rc_iter == read_counts_.end() ) {
			throw utility::excn::EXCN_Msg_Exception( "read_counts_ map does not contain an entry for file: " + filename );
		}
	} else {
		if ( rl_iter == read_limit_.end() ) {
			read_limit_[ filename ] = 0;
			read_counts_[ filename ] = 0;
			rl_iter = read_limit_.find( filename );
			rc_iter = read_counts_.find( filename );
		}
	}

	++(rc_iter->second);

	if ( delete_contents_at_nread_limit_ ) {

		if ( rl_iter->second == 1 && fc_iter == file_contents_.end() ) {
			return utility::file_contents( filename );
		}

		if ( rc_iter->second == 1 && fc_iter == file_contents_.end() ) {
			// this is the first time we're reading the file, so we have to load it.
			file_contents_[ filename ] = utility::file_contents( filename );
			fc_iter = file_contents_.find( filename );
		}

		if ( rl_iter->second != 0 && rl_iter->second == rc_iter->second ) {
			if ( fc_iter == file_contents_.end() ) {
				throw utility::excn::EXCN_Msg_Exception( "file-contents map does not contain an entry for file that has been read before: " + filename );
			}
			std::string fc = fc_iter->second;
			if ( delete_contents_at_nread_limit_ ) file_contents_.erase( fc_iter );
			return fc;
		}
	}

	// otherwise, we're not deleting the contents of the file

	if ( fc_iter == file_contents_.end() ) {
		std::string fc = utility::file_contents( filename );
		file_contents_[ filename ] = fc;
		return fc;
	} else {
		return fc_iter->second;
	}
}

/// @throws Relies on utility::file_contents which throws a utility::excn::EXCN_Msg_Exception
/// if the desired file cannot be found.
std::string const &
FileContentsMap::get_file_contents_ref( std::string const & filename )
{
	std::map< std::string, platform::Size >::iterator rc_iter = read_counts_.find( filename );
	std::map< std::string, platform::Size >::iterator rl_iter = read_limit_.find( filename );
	std::map< std::string, std::string >::iterator fc_iter = file_contents_.find( filename );

	if ( rc_iter == read_counts_.end() ) {
		read_counts_[ filename ] = 1;
	} else {
		++(rc_iter->second);
	}

	if ( rl_iter == read_limit_.end() ) {
		// mark an indefinite read limit for this file.
		read_limit_[ filename ] = 0;
	}

	if ( fc_iter == file_contents_.end() ) {
		file_contents_[ filename ] = utility::file_contents( filename );
		fc_iter = file_contents_.find( filename );
	}

	return fc_iter->second;

}


}
}
