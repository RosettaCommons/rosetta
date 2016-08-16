// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/file/PathName.hh
/// @brief  Path name class supporting Windows and UN*X/Linux format names
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit headers
#include <utility/file/PathName.hh>
#include <utility/exit.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <algorithm>
#include <cctype>
#include <iostream> // PB fix gcc 4.2.1 error?

namespace utility {
namespace file {


/// @brief Volume assignment
PathName &
PathName::vol( std::string const & vol_a )
{
	if ( platform::file::VOLUME_USED ) {
		vol_.clear();
		runtime_assert( ( vol_a.length() == 2 ) && ( std::isalpha( vol_a[ 0 ] ) ) && ( vol_a[ 1 ] == ':' ) );
		vol_.assign( vol_a, 0, 2 );
	}
	return *this;
}


/// @brief Path assignment
PathName &
PathName::path( std::string const & path_a )
{
	path_ = path_a;
	char const not_sep( platform::file::PATH_SEPARATOR == '/' ? '\\' : '/' );
	std::replace( path_.begin(), path_.end(), not_sep, platform::file::PATH_SEPARATOR ); // Convert to platform separators
	if ( ( ! path_.empty() ) && ( path_[ path_.length() - 1 ] != platform::file::PATH_SEPARATOR ) ) { // Append separator
		path_ += platform::file::PATH_SEPARATOR;
	}
	return *this;
}


/// @brief Absolute path?
bool
PathName::absolute() const
{
	return ( ( ! path_.empty() ) && ( path_[ 0 ] == platform::file::PATH_SEPARATOR ) );
}


/// @brief Relative path?
bool
PathName::relative() const
{
	return !absolute();
}


/// @details For relative paths, this will stop when it reaches the empty string --
/// no "../" will be prepended to possibly go further up the directory tree.
PathName
PathName::parent() const
{
	PathName newpath( name() );
	if ( path_.length() >= 2 ) {
		// could be "/a/" or "a/"
		std::string::size_type second_to_last_sep = path_.rfind( platform::file::PATH_SEPARATOR, path_.length() - 2 );
		if ( second_to_last_sep == std::string::npos ) newpath.path( "" );
		else newpath.path( path_.substr(0, second_to_last_sep+1) );
	}
	// else the path is "" or "/", neither of which needs modification
	return newpath;
}

/// @brief Input from stream
std::istream &
operator >>( std::istream & stream, PathName & name )
{
	std::string name_string;
	stream >> name_string;
	name.assign( name_string );
	return stream;
}


/// @brief Output to stream
std::ostream &
operator <<( std::ostream & stream, PathName const & name )
{
	stream << name.name();
	return stream;
}


/// @brief Parse and assign a path name string
void
PathName::assign( std::string const & name_string )
{
	using std::string;

	// Local copy to cut up as we parse
	string cut( name_string );

	// Volume
	vol_.clear();
	if ( ( std::isalpha( cut[ 0 ] ) ) && ( cut[ 1 ] == ':' ) ) {
		if ( platform::file::VOLUME_USED ) vol_.assign( cut, 0, 2 );
		cut.erase( 0, 2 ); // Remove volume part
	}

	// Path
	path_ = cut; // Path is whatever is left
	char const not_sep( platform::file::PATH_SEPARATOR == '/' ? '\\' : '/' );
	std::replace( path_.begin(), path_.end(), not_sep, platform::file::PATH_SEPARATOR ); // Convert to platform separators
	if ( ( ! path_.empty() ) && ( path_[ path_.length() - 1 ] != platform::file::PATH_SEPARATOR ) ) { // Append separator
		path_ += platform::file::PATH_SEPARATOR;
	}

} // assign


/// @brief Characters equal case-insensitively?
inline
bool
char_equali( char c, char d )
{
	return ( std::tolower( c ) == std::tolower( d ) );
}


/// @brief Strings equal case-insensitively?
inline
bool
equali( std::string const & s, std::string const & t )
{
	return (
		( s.size() == t.size() ) &&
		( std::equal( s.begin(), s.end(), t.begin(), char_equali ) ) );
}


/// @brief PathNames equal on this platform?
bool
PathName::equal( PathName const & name1, PathName const & name2 )
{
	if ( platform::file::CASE_SENSITIVE ) { // Case-sensitive
		return (
			( name1.vol_  == name2.vol_  ) &&
			( name1.path_ == name2.path_ ) );
	} else { // Case-insensitive
		return (
			equali( name1.vol_ , name2.vol_  ) &&
			equali( name1.path_, name2.path_ ) );
	}
}


} // namespace file
} // namespace utility
