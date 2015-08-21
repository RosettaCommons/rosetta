// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/FileName.hh
/// @brief  File name class supporting Windows and UN*X/Linux format names
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Converts to local platform format


#ifndef INCLUDED_utility_file_FileName_hh
#define INCLUDED_utility_file_FileName_hh


// Unit headers
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.hh>

// C++ headers
#include <iosfwd>
#include <string>
#include <vector>


namespace utility {
namespace file {


/// @brief File name class supporting Windows and UN*X/Linux format names
class FileName
{


public: // Creation


	/// @brief Default constructor
	inline
	FileName()
	{}


	/// @brief String constructor
	inline
	FileName( std::string const & name_string )
	{
		assign( name_string );
	}


	/// @brief Uses filename from string but path and vol from PathName.
	inline
	FileName( std::string const & name_string, PathName const & use_path )
	{
		assign( name_string );
		path( use_path.path() );
		vol( use_path.vol() );
	}

	/// @brief FileName vector constructor
	inline
	FileName( std::vector<FileName> const & file_names)
	{
		std::vector<FileName>::const_iterator begin= file_names.begin();

		erase();
		if ( begin == file_names.end() ) return;

		vol_.assign(begin->vol_);
		path_.assign(begin->path_);
		ver_.assign(begin->ver_);
		ext_.assign(begin->ext_);
		base_.append(begin->base_);

		for ( ++begin ; begin != file_names.end(); ++begin ) {
			base_.append("_");
			base_.append(begin->base_);
		}
	}

	/// @brief Destructor
	inline
	~FileName()
	{}


public: // Assignment


	/// @brief String assignment
	inline
	FileName &
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
	FileName &
	operator ()( FileName const & name )
	{
		operator =( name );
		return *this;
	}


	/// @brief Functor string assignment
	inline
	FileName &
	operator ()( std::string const & name_string )
	{
		assign( name_string );
		return *this;
	}


	/// @brief Volume assignment
	FileName &
	vol( std::string const & vol_a );


	/// @brief Volume assignment
	inline
	FileName &
	volume( std::string const & vol_a )
	{
		return vol( vol_a );
	}


	/// @brief Path assignment
	FileName &
	path( std::string const & path_a );


	/// @brief Base assignment
	inline
	FileName &
	base( std::string const & base_a )
	{
		base_ = base_a;
		return *this;
	}


	/// @brief Extension assignment
	inline
	FileName &
	ext( std::string const & ext_a )
	{
		ext_ = ( ( ! ext_a.empty() ) && ( ext_a[ 0 ] == '.' ) ? ext_a.substr( 1 ) : ext_a );
		return *this;
	}


	/// @brief Extension assignment
	inline
	FileName &
	extension( std::string const & ext_a )
	{
		ext_ = ( ( ! ext_a.empty() ) && ( ext_a[ 0 ] == '.' ) ? ext_a.substr( 1 ) : ext_a );
		return *this;
	}


	/// @brief Version assignment
	inline
	FileName &
	version( std::string const & ver_a )
	{
		ver_ = ( ( ! ver_a.empty() ) && ( ver_a[ 0 ] == '~' ) ? ver_a.substr( 1 ) : ver_a );
		return *this;
	}


	/// @brief Version assignment
	inline
	FileName &
	ver( std::string const & ver_a )
	{
		ver_ = ( ( ! ver_a.empty() ) && ( ver_a[ 0 ] == '~' ) ? ver_a.substr( 1 ) : ver_a );
		return *this;
	}


	/// @brief Change to local name (without volume or path)
	inline
	FileName &
	to_local_name()
	{
		vol_.clear();
		path_.clear();
		return *this;
	}


	/// @brief Change to bare name (without volume or path or version)
	inline
	FileName &
	to_bare_name()
	{
		vol_.clear();
		path_.clear();
		ver_.clear();
		return *this;
	}


	/// @brief Erase the file name
	inline
	FileName &
	erase()
	{
		vol_.clear();
		path_.clear();
		base_.clear();
		ext_.clear();
		ver_.clear();
		return *this;
	}


	/// @brief Clear the file name
	inline
	FileName &
	clear()
	{
		vol_.clear();
		path_.clear();
		base_.clear();
		ext_.clear();
		ver_.clear();
		return *this;
	}


public: // Properties


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return (
			( vol_.empty() ) &&
			( path_.empty() ) &&
			( base_.empty() ) &&
			( ext_.empty() ) &&
			( ver_.empty() ) );
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


	/// @brief Base
	inline
	std::string const &
	base() const
	{
		return base_;
	}


	/// @brief Extension
	inline
	std::string const &
	ext() const
	{
		return ext_;
	}


	/// @brief Extension with separator
	inline
	std::string
	extension() const
	{
		return ( ext_.empty() ? std::string() : '.' + ext_ );
	}


	/// @brief Version
	inline
	std::string const &
	ver() const
	{
		return ver_;
	}


	/// @brief Version with separator
	inline
	std::string
	version() const
	{
		return ( ver_.empty() ? std::string() : '~' + ver_ );
	}


	/// @brief File name string
	inline
	std::string
	name() const
	{
		return
			vol_ +
			path_ +
			base_ +
			extension() +
			version();
	}


	/// @brief Functor
	inline
	std::string
	operator ()() const
	{
		return name();
	}


	/// @brief Local name (without volume or path)
	inline
	std::string
	local_name() const
	{
		return
			base_ +
			extension() +
			version();
	}


	/// @brief Bare name (without volume or path or version)
	inline
	std::string
	bare_name() const
	{
		return
			base_ +
			extension();
	}


public: // Comparison


	/// @brief FileName == FileName
	friend
	inline
	bool
	operator ==( FileName const & name1, FileName const & name2 )
	{
		return equal( name1, name2 );
	}


	/// @brief FileName != FileName
	friend
	inline
	bool
	operator !=( FileName const & name1, FileName const & name2 )
	{
		return ( ! equal( name1, name2 ) );
	}


	/// @brief FileName < FileName
	friend
	inline
	bool
	operator <( FileName const & name1, FileName const & name2 )
	{
		return ( name1() < name2() );
	}


	/// @brief FileName <= FileName
	friend
	inline
	bool
	operator <=( FileName const & name1, FileName const & name2 )
	{
		return ( name1() <= name2() );
	}


	/// @brief FileName >= FileName
	friend
	inline
	bool
	operator >=( FileName const & name1, FileName const & name2 )
	{
		return ( name1() >= name2() );
	}


	/// @brief FileName > FileName
	friend
	inline
	bool
	operator >( FileName const & name1, FileName const & name2 )
	{
		return ( name1() > name2() );
	}


public: // I/O


	/// @brief Input from stream
	friend
	std::istream &
	operator >>( std::istream & stream, FileName & name );


	/// @brief Output to stream
	friend
	std::ostream &
	operator <<( std::ostream & stream, FileName const & name );


private: // Functions


	/// @brief Parse and assign a file name string
	void
	assign( std::string const & name_string );


private: // Static functions


	/// @brief FileNames equal on this platform?
	static
	bool
	equal( FileName const & name1, FileName const & name2 );


private: // Fields


	/// @brief Volume (X:) (Windows only)
	std::string vol_;

	/// @brief Path (with trailing separator)
	std::string path_;

	/// @brief Base name (without path or extension or version)
	std::string base_;

	/// @brief Extension (without separator)
	std::string ext_;

	/// @brief Version (without separator)
	std::string ver_;


}; // FileName


// Friend function namespace declarations


/// @brief FileName == FileName
#ifndef __clang__
bool
operator ==( FileName const & name1, FileName const & name2 );
#endif

/// @brief FileName != FileName
#ifndef __clang__
bool
operator !=( FileName const & name1, FileName const & name2 );
#endif

/// @brief FileName < FileName
#ifndef __clang__
bool
operator <( FileName const & name1, FileName const & name2 );
#endif

/// @brief FileName <= FileName
#ifndef __clang__
bool
operator <=( FileName const & name1, FileName const & name2 );
#endif


/// @brief FileName >= FileName
#ifndef __clang__
bool
operator >=( FileName const & name1, FileName const & name2 );
#endif

/// @brief FileName > FileName
#ifndef __clang__
bool
operator >( FileName const & name1, FileName const & name2 );
#endif

} // namespace file
} // namespace utility


#endif // INCLUDED_utility_file_FileName_HH
