// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/file/FileName.hh
/// @brief  File name class supporting Windows and UN*X/Linux format names
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit headers
#include <utility/file/FileName.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <algorithm>
#include <cctype>
#include <iostream> // PB fix gcc 4.2.1 error?

namespace utility {
namespace file {


/// @brief Volume assignment
FileName &
FileName::vol( std::string const & vol_a )
{
	if ( platform::file::VOLUME_USED ) {
		vol_ = vol_a;
		//vol_.clear();
		//runtime_assert( ( vol_a.length() == 2 ) && ( std::isalpha( vol_a[ 0 ] ) ) && ( vol_a[ 1 ] == ':' ) );
		//vol_.assign( vol_a, 0, 2 );
	}
	return *this;
}


/// @brief Path assignment
FileName &
FileName::path( std::string const & path_a )
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
FileName::absolute() const
{
	return ( ( ! path_.empty() ) && ( path_[ 0 ] == platform::file::PATH_SEPARATOR ) );
}


/// @brief Relative path?
bool
FileName::relative() const
{
	return !absolute();
}


/// @brief Input from stream
std::istream &
operator >>( std::istream & stream, FileName & name )
{
	std::string name_string;
	stream >> name_string;
	name.assign( name_string );
	return stream;
}


/// @brief Output to stream
std::ostream &
operator <<( std::ostream & stream, FileName const & name )
{
	stream << name.name();
	return stream;
}


/// @brief Parse and assign a file name string
void
FileName::assign( std::string const & name_string )
{
	using std::string;
	typedef  std::string::size_type  size_type;

	// Local copy to cut up as we parse
	string cut( name_string );

	// Volume
	vol_.clear();
	if ( ( std::isalpha( cut[ 0 ] ) ) && ( cut[ 1 ] == ':' ) ) {
		if ( platform::file::VOLUME_USED ) vol_.assign( cut, 0, 2 );
		cut.erase( 0, 2 ); // Remove volume part
	}

	// Path
	path_.clear();
	size_type const i_path( cut.find_last_of( "/\\" ) ); // Accept Windows or UN*X separators
	if ( i_path != string::npos ) { // Path present
		path_.assign( cut, 0, i_path + 1 );
		char const not_sep( platform::file::PATH_SEPARATOR == '/' ? '\\' : '/' );
		std::replace( path_.begin(), path_.end(), not_sep, platform::file::PATH_SEPARATOR ); // Convert to platform separators
		cut.erase( 0, i_path + 1 ); // Remove path part
	}

	// Version
	ver_.clear();
	size_type const i_ver = cut.find_last_of( '~' );
	if ( i_ver != string::npos ) { // Version present
		if ( i_ver + 1 < cut.length() ) ver_.assign( cut, i_ver + 1, string::npos );
		cut.erase( i_ver );
	}

	// Extension
	ext_.clear();
	size_type i_ext( cut.find_last_of( '.' ) );
	if ( i_ext != string::npos ) { // Extension present
		if ( i_ext + 1 < cut.length() ) ext_.assign( cut, i_ext + 1, string::npos );
		cut.erase( i_ext );
		if ( ext_=="gz" ) { ///special case for dealing with files that have 2 extensions, ".pdb.gz"
			i_ext= cut.find_last_of( '.' ) ;
			if ( i_ext != string::npos ) { // Extension present
				if ( i_ext + 1 < cut.length() ) {
					ext_.insert( 0, ".");
					ext_.insert( 0, cut, i_ext+1, string::npos);
				}
				cut.erase( i_ext );
			}
		}
	}

	// Base
	base_ = cut; // Base part is whatever is left

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


/// @brief FileNames equal on this platform?
bool
FileName::equal( FileName const & name1, FileName const & name2 )
{
	if ( platform::file::CASE_SENSITIVE ) { // Case-sensitive
		return (
			( name1.vol_  == name2.vol_  ) &&
			( name1.path_ == name2.path_ ) &&
			( name1.base_ == name2.base_ ) &&
			( name1.ext_  == name2.ext_  ) &&
			( name1.ver_  == name2.ver_  ) );
	} else { // Case-insensitive
		return (
			equali( name1.vol_ , name2.vol_  ) &&
			equali( name1.path_, name2.path_ ) &&
			equali( name1.base_, name2.base_ ) &&
			equali( name1.ext_ , name2.ext_  ) &&
			equali( name1.ver_ , name2.ver_  ) );
	}
}


} // namespace file
} // namespace utility
