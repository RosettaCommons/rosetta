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
///
/// @remarks
///  @li Converts to local platform format


#ifndef INCLUDED_utility_file_PathName_hh
#define INCLUDED_utility_file_PathName_hh


// Unit headers
#include <utility/file/PathName.fwd.hh>

// C++ headers
#include <iosfwd>
#include <string>


namespace utility {
namespace file {


/// @brief Path name class supporting Windows and UN*X/Linux format names
class PathName
{


public: // Creation


	/// @brief Default constructor
	inline
	PathName()
	{}


	/// @brief String constructor
	inline
	PathName( std::string const & name_string )
	{
		assign( name_string );
	}


	/// @brief Destructor
	inline
	~PathName()
	{}


public: // Assignment


	/// @brief String assignment
	inline
	PathName &
	operator =( std::string const & name_string )
	{
		assign( name_string );
		return *this;
	}


public: // Conversion


	/// @brief String conversion
	inline
	operator std::string() const
	{
		return name();
	}


public: // Methods


	/// @brief Functor copy assignment
	inline
	PathName &
	operator ()( PathName const & name )
	{
		operator =( name );
		return *this;
	}


	/// @brief Functor string assignment
	inline
	PathName &
	operator ()( std::string const & name_string )
	{
		assign( name_string );
		return *this;
	}


	/// @brief Volume assignment
	PathName &
	vol( std::string const & vol_a );


	/// @brief Volume assignment
	inline
	PathName &
	volume( std::string const & vol_a )
	{
		return vol( vol_a );
	}


	/// @brief Path assignment
	PathName &
	path( std::string const & path_a );


	/// @brief Erase the path name
	inline
	PathName &
	erase()
	{
		vol_.clear();
		path_.clear();
		return *this;
	}


	/// @brief Clear the path name
	inline
	PathName &
	clear()
	{
		vol_.clear();
		path_.clear();
		return *this;
	}


	/// @brief Returns the parent of this directory, or itself if no parent is available.
	PathName
	parent() const;


public: // Properties


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return (
			( vol_.empty() ) &&
			( path_.empty() ) );
	}


	/// @brief Absolute path?
	bool
	absolute() const;


	/// @brief Relative path?
	bool
	relative() const;


	/// @brief Volume
	inline
	std::string const &
	vol() const
	{
		return vol_;
	}


	/// @brief Volume
	inline
	std::string const &
	volume() const
	{
		return vol_;
	}


	/// @brief Path
	inline
	std::string const &
	path() const
	{
		return path_;
	}


	/// @brief Path name string
	inline
	std::string
	name() const
	{
		return
			vol_ +
			path_;
	}


	/// @brief Functor
	inline
	std::string
	operator ()() const
	{
		return name();
	}


public: // Comparison


	/// @brief PathName == PathName
	friend
	inline
	bool
	operator ==( PathName const & name1, PathName const & name2 )
	{
		return equal( name1, name2 );
	}


	/// @brief PathName != PathName
	friend
	inline
	bool
	operator !=( PathName const & name1, PathName const & name2 )
	{
		return ( ! equal( name1, name2 ) );
	}


	/// @brief PathName < PathName
	friend
	inline
	bool
	operator <( PathName const & name1, PathName const & name2 )
	{
		return ( name1() < name2() );
	}


	/// @brief PathName <= PathName
	friend
	inline
	bool
	operator <=( PathName const & name1, PathName const & name2 )
	{
		return ( name1() <= name2() );
	}


	/// @brief PathName >= PathName
	friend
	inline
	bool
	operator >=( PathName const & name1, PathName const & name2 )
	{
		return ( name1() >= name2() );
	}


	/// @brief PathName > PathName
	friend
	inline
	bool
	operator >( PathName const & name1, PathName const & name2 )
	{
		return ( name1() > name2() );
	}


public: // I/O


	/// @brief Input from stream
	friend
	std::istream &
	operator >>( std::istream & stream, PathName & name );


	/// @brief Output to stream
	friend
	std::ostream &
	operator <<( std::ostream & stream, PathName const & name );


private: // Functions


	/// @brief Parse and assign a path name string
	void
	assign( std::string const & name_string );


private: // Static functions


	/// @brief PathNames equal on this platform?
	static
	bool
	equal( PathName const & name1, PathName const & name2 );


private: // Fields


	/// @brief Volume (X:) (Windows only)
	std::string vol_;

	/// @brief Path (with trailing separator)
	std::string path_;


}; // PathName


// Friend function namespace declarations


/// @brief PathName == PathName
#ifndef __clang__
bool
operator ==( PathName const & name1, PathName const & name2 );
#endif

/// @brief PathName != PathName
#ifndef __clang__
bool
operator !=( PathName const & name1, PathName const & name2 );
#endif

/// @brief PathName < PathName
#ifndef __clang__
bool
operator <( PathName const & name1, PathName const & name2 );
#endif

/// @brief PathName <= PathName
#ifndef __clang__
bool
operator <=( PathName const & name1, PathName const & name2 );
#endif

/// @brief PathName >= PathName
#ifndef __clang__
bool
operator >=( PathName const & name1, PathName const & name2 );
#endif

/// @brief PathName > PathName
#ifndef __clang__
bool
operator >( PathName const & name1, PathName const & name2 );
#endif

} // namespace file
} // namespace utility


#endif // INCLUDED_utility_file_PathName_HH
