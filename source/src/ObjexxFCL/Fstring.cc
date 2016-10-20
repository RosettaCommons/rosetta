// Fstring: Fixed-Length Fortran-Compatible String and Substring
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/Fstring.hh>
#include <ObjexxFCL/char.functions.hh>

// C++ Headers
#include <algorithm>
//#include <iostream>


namespace ObjexxFCL {


// Constants
char const SPACE( ' ' );
char const TAB( '\t' );
char const NULLC( '\000' );
std::string const WHITESPACE( " \t\0", 3 );


// Fstring


	/// @brief Copy Constructor
	Fstring::Fstring( Fstring const & s ) :
		len_( s.len_ ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memcpy( str_, s.str_, len_ );
		}
	}


	/// @brief string Constructor
	Fstring::Fstring( std::string const & s ) :
		len_( s.length() ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		s.copy( str_, len_ );
	}


	/// @brief cstring Constructor
	Fstring::Fstring( c_cstring const s ) :
		len_( std::strlen( s ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memcpy( str_, s, len_ );
		}
	}


	/// @brief Length Constructor
	Fstring::Fstring( short int const len_a ) :
		len_( static_cast< size_type >( std::max( static_cast< int >( len_a ), 0 ) ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}


	/// @brief Length Constructor
	Fstring::Fstring( int const len_a ) :
		len_( static_cast< size_type >( std::max( len_a, 0 ) ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}


	/// @brief Length Constructor
	Fstring::Fstring( long int const len_a ) :
		len_( static_cast< size_type >( std::max( len_a, 0l ) ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}


	/// @brief Length Constructor
	Fstring::Fstring( unsigned short int const len_a ) :
		len_( static_cast< size_type >( len_a ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}


	/// @brief Length Constructor
	Fstring::Fstring( unsigned int const len_a ) :
		len_( static_cast< size_type >( len_a ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}


	/// @brief Length Constructor
	Fstring::Fstring( unsigned long int const len_a ) :
		len_( static_cast< size_type >( len_a ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}

	Fstring::Fstring( unsigned long long const len_a ) :
		len_( static_cast< size_type >( len_a ) ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
	}

	/// @brief Length + Fstring Constructor
	Fstring::Fstring( size_type const len_a, Fstring const & s ) :
		len_( len_a ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > s.len_ ) {
			if ( s.len_ > 0 ) std::memcpy( str_, s.str_, s.len_ );
			std::memset( str_ + s.len_, SPACE, len_ - s.len_ ); // Space pad
		} else if ( len_ > 0 ) {
			std::memcpy( str_, s.str_, len_ );
		}
	}


	/// @brief Length + string Constructor
	Fstring::Fstring( size_type const len_a, std::string const & s ) :
		len_( len_a ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		size_type const s_len( s.length() );
		if ( len_ > s_len ) {
			if ( s_len > 0 ) s.copy( str_, s_len );
			std::memset( str_ + s_len, SPACE, len_ - s_len ); // Space pad
		} else if ( len_ > 0 ) {
			s.copy( str_, len_ );
		}
	}


	/// @brief Length + cstring Constructor
	Fstring::Fstring( size_type const len_a, c_cstring const s ) :
		len_( len_a ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		size_type const s_len( std::strlen( s ) );
		if ( len_ > s_len ) {
			if ( s_len > 0 ) std::memcpy( str_, s, s_len );
			std::memset( str_ + s_len, SPACE, len_ - s_len ); // Space pad
		} else if ( len_ > 0 ) {
			std::memcpy( str_, s, len_ );
		}
	}


	/// @brief Length + char Constructor
	/// @note  Fills with specified char => Use Fstring( len_a, "c" ) for space-padded single character
	Fstring::Fstring( size_type const len_a, char const c ) :
		len_( len_a ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, c, len_ );
		}
	}


	/// @brief Length + Initializer Constructor
	Fstring::Fstring( size_type const len_a, initializer_function init ) :
		len_( len_a ),
		str_( len_ > 0 ? new char[ len_ ] : 0 ),
		c_str_( 0 ),
		sub_( false )
	{
		if ( len_ > 0 ) {
			std::memset( str_, SPACE, len_ );
		}
		init( *this );
	}


	/// @brief Substring Range Constructor
	Fstring::Fstring( Fstring const & s, size_type const i, size_type const j ) :
		len_( i <= std::min( j, s.len_ ) ? std::min( j, s.len_ ) - i + 1u : static_cast< size_type >( 0 ) ),
		str_( i <= s.len_ ? s.str_ + i - 1u : s.str_ ),
		c_str_( 0 ),
		sub_( true )
	{
		assert( i > 0 );
	}


	/// @brief Substring Tail Constructor
	Fstring::Fstring( Fstring const & s, size_type const i ) :
		len_( i <= s.len_ ? s.len_ - i + 1u : static_cast< size_type >( 0 ) ),
		str_( i <= s.len_ ? s.str_ + i - 1u : s.str_ ),
		c_str_( 0 ),
		sub_( true )
	{
		assert( i > 0 );
	}


	/// @brief Copy Assignment
	Fstring &
	Fstring::operator =( Fstring const & s )
	{
		if ( this != &s ) {
			if ( len_ > s.len_ ) {
				if ( s.len_ > 0 ) std::memmove( str_, s.str_, s.len_ );
				std::memset( str_ + s.len_, SPACE, len_ - s.len_ ); // Space pad
			} else if ( len_ > 0 ) {
				std::memmove( str_, s.str_, len_ );
			} else if ( ( ! str_ ) && ( ! sub_ ) ) { // Uninitialized
				len_ = s.len_;
				str_ = ( len_ > 0 ? new char[ len_ ] : 0 );
				if ( len_ > 0 ) {
					std::memcpy( str_, s.str_, len_ );
				}
			}
		}
		return *this;
	}


	/// @brief = string
	Fstring &
	Fstring::operator =( std::string const & s )
	{
		size_type const s_len( s.length() );
		if ( len_ > s_len ) {
			if ( s_len > 0 ) s.copy( str_, s_len );
			std::memset( str_ + s_len, SPACE, len_ - s_len ); // Space pad
		} else if ( len_ > 0 ) {
			s.copy( str_, len_ );
		} else if ( ( ! str_ ) && ( ! sub_ ) ) { // Uninitialized
			len_ = s_len;
			str_ = ( len_ > 0 ? new char[ len_ ] : 0 );
			s.copy( str_, len_ );
		}
		return *this;
	}


	/// @brief = cstring
	Fstring &
	Fstring::operator =( c_cstring const s )
	{
		size_type const s_len( std::strlen( s ) );
		if ( len_ > s_len ) {
			if ( s_len > 0 ) std::memmove( str_, s, s_len );
			std::memset( str_ + s_len, SPACE, len_ - s_len ); // Space pad
		} else if ( len_ > 0 ) {
			std::memmove( str_, s, len_ );
		} else if ( ( ! str_ ) && ( ! sub_ ) ) { // Uninitialized
			len_ = s_len;
			str_ = ( len_ > 0 ? new char[ len_ ] : 0 );
			if ( len_ > 0 ) {
				std::memcpy( str_, s, len_ );
			}
		}
		return *this;
	}


	/// @brief = char
	Fstring &
	Fstring::operator =( char const c )
	{
		if ( len_ > 0 ) {
			str_[ 0 ] = c;
			if ( len_ > 1 ) std::memset( str_ + 1, SPACE, len_ - 1 ); // Space pad
		} else if ( ( ! str_ ) && ( ! sub_ ) ) { // Uninitialized
			len_ = 1;
			str_ = new char[ 1 ];
			str_[ 0 ] = c;
		}
		return *this;
	}


	/// @brief Has an Fstring?
	bool
	Fstring::has( Fstring const & s ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			if ( (*this)( i ) == s ) return true;
		}
		return false; // No matches
	}


	/// @brief Has a string?
	bool
	Fstring::has( std::string const & s ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			if ( (*this)( i ) == s ) return true;
		}
		return false; // No matches
	}


	/// @brief Has a cstring?
	bool
	Fstring::has( c_cstring const s ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			if ( (*this)( i ) == s ) return true;
		}
		return false; // No matches
	}


	/// @brief Has a Character?
	bool
	Fstring::has( char const c ) const
	{
		for ( size_type i = 0; i < len_; ++i ) {
			if ( str_[ i ] == c ) return true;
		}
		return false; // No matches
	}


	/// @brief Has Any Character of an Fstring?
	bool
	Fstring::has_any_of( Fstring const & s ) const
	{
		for ( size_type i = 0; i < len_; ++i ) {
			for ( size_type j = 0, e = s.len_; j < e; ++j ) {
				if ( str_[ i ] == s.str_[ j ] ) return true;
			}
		}
		return false; // No matches
	}


	/// @brief Has Any Character of a string?
	bool
	Fstring::has_any_of( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		for ( size_type i = 0; i < len_; ++i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i ] == s[ j ] ) return true;
			}
		}
		return false; // No matches
	}


	/// @brief Has Any Character of a cstring?
	bool
	Fstring::has_any_of( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		for ( size_type i = 0; i < len_; ++i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i ] == s[ j ] ) return true;
			}
		}
		return false; // No matches
	}


	/// @brief Has a Character?
	bool
	Fstring::has_any_of( char const c ) const
	{
		for ( size_type i = 0; i < len_; ++i ) {
			if ( str_[ i ] == c ) return true;
		}
		return false; // No matches
	}


	/// @brief Has a Prefix Case-Optionally?
	bool
	Fstring::has_prefix( Fstring const & s, bool const exact_case ) const
	{
		if ( s.len_ == 0 ) {
			return false;
		} else if ( len_ < s.len_ ) {
			return false;
		} else if ( exact_case ) {
			return ( (*this)( 1, s.len_ ) == s );
		} else {
			return ( lowercased()( 1, s.len_ ) == s.lowercased() );
		}
	}


	/// @brief Has a Prefix Case-Optionally?
	bool
	Fstring::has_prefix( c_cstring const s, bool const exact_case ) const
	{
		size_type const s_len( std::strlen( s ) );
		if ( s_len == 0 ) {
			return false;
		} else if ( len_ < s_len ) {
			return false;
		} else if ( exact_case ) {
			return ( (*this)( 1, s_len ) == s );
		} else {
			return ( lowercased()( 1, s_len ) == Fstring( s ).lowercase() );
		}
	}


	/// @brief Length Space-Trimmed
	Fstring::size_type
	Fstring::len_trim() const
	{
		for ( size_type i = len_; i > 0; --i ) {
			if ( str_[ i - 1 ] != SPACE ) return i;
		}
		return 0;
	}


	/// @brief Length Whitespace-Trimmed
	Fstring::size_type
	Fstring::len_trim_whitespace() const
	{
		for ( size_type i = len_; i > 0; --i ) {
			char const c( str_[ i - 1 ] );
			if ( ( c != SPACE ) && ( c != TAB ) && ( c != NULLC ) ) return i;
		}
		return 0;
	}


	/// @brief Find First Occurrence of a Whitespace Character
	Fstring::size_type
	Fstring::find_whitespace() const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			char const c( str_[ i - 1 ] );
			if ( ( c == SPACE ) || ( c == TAB ) || ( c == NULLC ) ) return i;
		}
		return 0; // All are non-whitespace
	}


	/// @brief Find First Occurrence of a Non-Whitespace Character
	Fstring::size_type
	Fstring::find_non_whitespace() const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			char const c( str_[ i - 1 ] );
			if ( ( c != SPACE ) && ( c != TAB ) && ( c != NULLC ) ) return i;
		}
		return 0; // All are whitespace
	}


	/// @brief Find Last Occurrence of a Whitespace Character
	Fstring::size_type
	Fstring::find_last_whitespace() const
	{
		for ( size_type i = len_; i > 0; --i ) {
			char const c( str_[ i - 1 ] );
			if ( ( c == SPACE ) || ( c == TAB ) || ( c == NULLC ) ) return i;
		}
		return 0; // All are non-whitespace
	}


	/// @brief Find Last Occurrence of a Non-Whitespace Character
	Fstring::size_type
	Fstring::find_last_non_whitespace() const
	{
		for ( size_type i = len_; i > 0; --i ) {
			char const c( str_[ i - 1 ] );
			if ( ( c != SPACE ) && ( c != TAB ) && ( c != NULLC ) ) return i;
		}
		return 0; // All are whitespace
	}


	/// @brief Get Range of Whitespace-Trimmed Portion and Return its Length
	Fstring::size_type
	Fstring::trimmed_whitespace_range( size_type & b, size_type & e ) const
	{
		b = std::max( find_non_whitespace(), static_cast< size_type >( 1 ) );
		e = len_trim_whitespace();
		return e - b + 1;
	}


	/// @brief Find First Occurrence of an Fstring
	Fstring::size_type
	Fstring::find( Fstring const & s ) const
	{
		size_type const s_len( s.length() );
		if ( ( s_len > 0 ) && ( s_len <= len_ ) ) {
			for ( size_type i = 1, e = len_ - s_len + 1; i <= e; ++i ) {
				if ( (*this)( i, i + s_len - 1 ) == s ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of a string
	Fstring::size_type
	Fstring::find( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		if ( ( s_len > 0 ) && ( s_len <= len_ ) ) {
			for ( size_type i = 1, e = len_ - s_len + 1; i <= e; ++i ) {
				if ( (*this)( i, i + s_len - 1 ) == s ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of a cstring
	Fstring::size_type
	Fstring::find( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		if ( ( s_len > 0 ) && ( s_len <= len_ ) ) {
			for ( size_type i = 1, e = len_ - s_len + 1; i <= e; ++i ) {
				if ( (*this)( i, i + s_len - 1 ) == s ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of a Character
	Fstring::size_type
	Fstring::find( char const c ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			if ( str_[ i - 1 ] == c ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of an Fstring
	Fstring::size_type
	Fstring::find_last( Fstring const & s ) const
	{
		size_type const s_len( s.length() );
		if ( ( s_len > 0 ) && ( s_len <= len_ ) ) {
			for ( size_type i = len_ - s_len + 1; i > 0; --i ) {
				if ( (*this)( i, i + s_len - 1 ) == s ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of a string
	Fstring::size_type
	Fstring::find_last( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		if ( ( s_len > 0 ) && ( s_len <= len_ ) ) {
			for ( size_type i = len_ - s_len + 1; i > 0; --i ) {
				if ( (*this)( i, i + s_len - 1 ) == s ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of a cstring
	Fstring::size_type
	Fstring::find_last( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		if ( ( s_len > 0 ) && ( s_len <= len_ ) ) {
			for ( size_type i = len_ - s_len + 1; i > 0; --i ) {
				if ( (*this)( i, i + s_len - 1 ) == s ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of a Character
	Fstring::size_type
	Fstring::find_last( char const c ) const
	{
		for ( size_type i = len_; i > 0; --i ) {
			if ( str_[ i - 1 ] == c ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of Any Character of an Fstring
	Fstring::size_type
	Fstring::find_first_of( Fstring const & s ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			for ( size_type j = 0, e = s.len_; j < e; ++j ) {
				if ( str_[ i - 1 ] == s.str_[ j ] ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of Any Character of a string
	Fstring::size_type
	Fstring::find_first_of( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		for ( size_type i = 1; i <= len_; ++i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of Any Character of a cstring
	Fstring::size_type
	Fstring::find_first_of( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		for ( size_type i = 1; i <= len_; ++i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of a Character
	Fstring::size_type
	Fstring::find_first_of( char const c ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			if ( str_[ i - 1 ] == c ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of Any Character not of an Fstring
	Fstring::size_type
	Fstring::find_first_not_of( Fstring const & s ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			bool found( false );
			for ( size_type j = 0, e = s.len_; j < e; ++j ) {
				if ( str_[ i - 1 ] == s.str_[ j ] ) {
					found = true;
					break;
				}
			}
			if ( ! found ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of Any Character not of a string
	Fstring::size_type
	Fstring::find_first_not_of( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		for ( size_type i = 1; i <= len_; ++i ) {
			bool found( false );
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) {
					found = true;
					break;
				}
			}
			if ( ! found ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence of Any Character not of a cstring
	Fstring::size_type
	Fstring::find_first_not_of( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		for ( size_type i = 1; i <= len_; ++i ) {
			bool found( false );
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) {
					found = true;
					break;
				}
			}
			if ( ! found ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find First Occurrence not of a Character
	Fstring::size_type
	Fstring::find_first_not_of( char const c ) const
	{
		for ( size_type i = 1; i <= len_; ++i ) {
			if ( str_[ i - 1 ] != c ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of Any Character of an Fstring
	Fstring::size_type
	Fstring::find_last_of( Fstring const & s ) const
	{
		for ( size_type i = len_; i > 0; --i ) {
			for ( size_type j = 0, e = s.len_; j < e; ++j ) {
				if ( str_[ i - 1 ] == s.str_[ j ] ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of Any Character of a string
	Fstring::size_type
	Fstring::find_last_of( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		for ( size_type i = len_; i > 0; --i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of Any Character of a cstring
	Fstring::size_type
	Fstring::find_last_of( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		for ( size_type i = len_; i > 0; --i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) return i;
			}
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of a Character
	Fstring::size_type
	Fstring::find_last_of( char const c ) const
	{
		for ( size_type i = len_; i > 0; --i ) {
			if ( str_[ i - 1 ] == c ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of Any Character not of an Fstring
	Fstring::size_type
	Fstring::find_last_not_of( Fstring const & s ) const
	{
		for ( size_type i = len_; i > 0; --i ) {
			bool found( false );
			for ( size_type j = 0, e = s.len_; j < e; ++j ) {
				if ( str_[ i - 1 ] == s.str_[ j ] ) {
					found = true;
					break;
				}
			}
			if ( ! found ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of Any Character not of a string
	Fstring::size_type
	Fstring::find_last_not_of( std::string const & s ) const
	{
		std::string::size_type const s_len( s.length() );
		for ( size_type i = len_; i > 0; --i ) {
			bool found( false );
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) {
					found = true;
					break;
				}
			}
			if ( ! found ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence of Any Character not of a cstring
	Fstring::size_type
	Fstring::find_last_not_of( c_cstring const s ) const
	{
		size_type const s_len( std::strlen( s ) );
		for ( size_type i = len_; i > 0; --i ) {
			bool found( false );
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i - 1 ] == s[ j ] ) {
					found = true;
					break;
				}
			}
			if ( ! found ) return i;
		}
		return 0; // No matches
	}


	/// @brief Find Last Occurrence not of a Character
	Fstring::size_type
	Fstring::find_last_not_of( char const c ) const
	{
		for ( size_type i = len_; i > 0; --i ) {
			if ( str_[ i - 1 ] != c ) return i;
		}
		return 0; // No matches
	}


	/// @brief Lowercase
	Fstring &
	Fstring::lowercase()
	{
		for ( size_type i = 0; i < len_; ++i ) {
			str_[ i ] = std::tolower( str_[ i ] );
		}
		return *this;
	}


	/// @brief Uppercase
	Fstring &
	Fstring::uppercase()
	{
		for ( size_type i = 0; i < len_; ++i ) {
			str_[ i ] = std::toupper( str_[ i ] );
		}
		return *this;
	}


	/// @brief Left Justify
	Fstring &
	Fstring::left_justify()
	{
		for ( size_type i = 0; i < len_; ++i ) {
			if ( str_[ i ] != SPACE ) {
				if ( i > 0 ) {
					std::memmove( str_, str_ + i, len_ - i );
					std::memset( str_ + len_ - i, SPACE, i );
				}
				return *this;
			}
		}
		return *this;
	}


	/// @brief Right Justify
	Fstring &
	Fstring::right_justify()
	{
		for ( size_type i = len_; i > 0; --i ) {
			if ( str_[ i - 1 ] != SPACE ) {
				if ( i < len_ ) {
					std::memmove( str_ + len_ - i, str_, i );
					std::memset( str_, SPACE, len_ - i );
				}
				return *this;
			}
		}
		return *this;
	}


	/// @brief Center
	Fstring &
	Fstring::center()
	{
		left_justify();
		size_type const len_t( len_trim() );
		size_type const pad( ( len_ - len_t ) / 2 );
		if ( pad > 0 ) {
			std::memmove( str_ + pad, str_, len_t );
			std::memset( str_, SPACE, pad );
		}
		return *this;
	}


	/// @brief Compress Out Whitespace
	Fstring &
	Fstring::compress()
	{
		size_type j( 0 );
		for ( size_type i = 0; i < len_; ++i ) {
			char const c( str_[ i ] );
			if ( ( c != SPACE ) && ( c != TAB ) && ( c != NULLC ) ) str_[ j++ ] = c;
		}
		if ( j < len_ ) std::memset( str_ + j, SPACE, len_ - j );
		return *this;
	}


	/// @brief Trim Trailing Whitespace Replacing it with Space
	Fstring &
	Fstring::trim_whitespace()
	{
		if ( len_ > 0 ) {
			size_type const ie( len_trim_whitespace() );
			if ( ie < len_ ) {
				std::memset( str_ + ie, SPACE, len_ - ie ); // Space pad
			}
		}
		return *this;
	}


	/// @brief Strip Specified Characters from the Tails
	Fstring &
	Fstring::strip( std::string const & chars )
	{
		if ( len_ > 0 ) {
			size_type const ib( find_first_not_of( chars ) );
			if ( ib == 0 ) {
				clear();
			} else {
				size_type const ie( find_last_not_of( chars ) );
				assert( ie >= ib );
				size_type const len_sub( ie - ib + 1 );
				if ( len_sub < len_ ) {
					if ( ib > 1 ) std::memmove( str_, str_ + ib - 1, len_sub );
					std::memset( str_ + len_sub, SPACE, len_ - len_sub ); // Space pad
				}
			}
		}
		return *this;
	}


	/// @brief Strip Specified Characters from the Left Tail
	Fstring &
	Fstring::lstrip( std::string const & chars )
	{
		if ( len_ > 0 ) {
			size_type const ib( find_first_not_of( chars ) );
			if ( ib == 0 ) {
				clear();
			} else {
				size_type const len_sub( len_ - ib + 1 );
				if ( len_sub < len_ ) {
					if ( ib > 1 ) std::memmove( str_, str_ + ib - 1, len_sub );
					std::memset( str_ + len_sub, SPACE, len_ - len_sub ); // Space pad
				}
			}
		}
		return *this;
	}


	/// @brief Strip Specified Characters from the Right Tail
	Fstring &
	Fstring::rstrip( std::string const & chars )
	{
		if ( len_ > 0 ) {
			size_type const ie( find_last_not_of( chars ) );
			if ( ie == 0 ) {
				clear();
			} else if ( ie < len_ ) {
				std::memset( str_ + ie, SPACE, len_ - ie ); // Space pad
			}
		}
		return *this;
	}


	/// @brief Strip Space from the Tails
	Fstring &
	Fstring::strip()
	{
		if ( len_ > 0 ) {
			size_type const ib( find_first_not_of( SPACE ) );
			if ( ib == 0 ) {
				clear();
			} else {
				size_type const ie( find_last_not_of( SPACE ) );
				assert( ie >= ib );
				size_type const len_sub( ie - ib + 1 );
				if ( len_sub < len_ ) {
					if ( ib > 1 ) std::memmove( str_, str_ + ib - 1, len_sub );
					std::memset( str_ + len_sub, SPACE, len_ - len_sub ); // Space pad
				}
			}
		}
		return *this;
	}


	/// @brief Strip Space from the Left Tail
	Fstring &
	Fstring::lstrip()
	{
		if ( len_ > 0 ) {
			size_type const ib( find_first_not_of( SPACE ) );
			if ( ib == 0 ) {
				clear();
			} else {
				size_type const len_sub( len_ - ib + 1 );
				if ( len_sub < len_ ) {
					if ( ib > 1 ) std::memmove( str_, str_ + ib - 1, len_sub );
					std::memset( str_ + len_sub, SPACE, len_ - len_sub ); // Space pad
				}
			}
		}
		return *this;
	}


	/// @brief Strip Space from the Right Tail
	Fstring &
	Fstring::rstrip()
	{
		if ( len_ > 0 ) {
			size_type const ie( find_last_not_of( SPACE ) );
			if ( ie == 0 ) {
				clear();
			} else if ( ie < len_ ) {
				std::memset( str_ + ie, SPACE, len_ - ie ); // Space pad
			}
		}
		return *this;
	}


	/// @brief Strip Whitespace from the Tails
	Fstring &
	Fstring::strip_whitespace()
	{
		if ( len_ > 0 ) {
			size_type const ib( find_first_not_of( WHITESPACE ) );
			if ( ib == 0 ) {
				clear();
			} else {
				size_type const ie( find_last_not_of( WHITESPACE ) );
				assert( ie >= ib );
				size_type const len_sub( ie - ib + 1 );
				if ( len_sub < len_ ) {
					if ( ib > 1 ) std::memmove( str_, str_ + ib - 1, len_sub );
					std::memset( str_ + len_sub, SPACE, len_ - len_sub ); // Space pad
				}
			}
		}
		return *this;
	}


	/// @brief Strip Whitespace from the Left Tail
	Fstring &
	Fstring::lstrip_whitespace()
	{
		if ( len_ > 0 ) {
			size_type const ib( find_first_not_of( WHITESPACE ) );
			if ( ib == 0 ) {
				clear();
			} else {
				size_type const len_sub( len_ - ib + 1 );
				if ( len_sub < len_ ) {
					if ( ib > 1 ) std::memmove( str_, str_ + ib - 1, len_sub );
					std::memset( str_ + len_sub, SPACE, len_ - len_sub ); // Space pad
				}
			}
		}
		return *this;
	}


	/// @brief Strip Whitespace from the Right Tail
	Fstring &
	Fstring::rstrip_whitespace()
	{
		if ( len_ > 0 ) {
			size_type const ie( find_last_not_of( WHITESPACE ) );
			if ( ie == 0 ) {
				clear();
			} else if ( ie < len_ ) {
				std::memset( str_ + ie, SPACE, len_ - ie ); // Space pad
			}
		}
		return *this;
	}


	/// @brief Overlay an Fstring
	Fstring &
	Fstring::overlay( Fstring const & s, size_type const pos )
	{
		(*this)( pos, std::min( ( pos + s.len_ ) - 1, len_ ) ) = s;
		return *this;
	}


	/// @brief Overlay a string
	Fstring &
	Fstring::overlay( std::string const & s, size_type const pos )
	{
		(*this)( pos, std::min<size_type>( ( pos + s.length() ) - 1, len_ ) ) = s;
		return *this;
	}


	/// @brief Overlay a cstring
	Fstring &
	Fstring::overlay( c_cstring const s, size_type const pos )
	{
		(*this)( pos, std::min<size_type>( ( pos + std::strlen( s ) ) - 1, len_ ) ) = s;
		return *this;
	}


	/// @brief Null-Terminated cstring Copy of the Fstring that is Owned by the Fstring
	c_cstring
	Fstring::c_str() const
	{
		delete[] c_str_; c_str_ = new char[ len_ + 1 ];
		if ( len_ > 0 ) std::memmove( c_str_, str_, len_ ); // Copy the string data
		c_str_[ len_ ] = '\0'; // Null-terminate
		return c_str_;
	}


	/// @brief Whitespace-Trimmed Null-Terminated cstring Copy of the Fstring that is Owned by the Fstring
	/// @note This shares data/pointer with c_str()
	c_cstring
	Fstring::t_str() const
	{
		size_type const len_trim_whitespace_( len_trim_whitespace() );
		delete[] c_str_; c_str_ = new char[ len_trim_whitespace_ + 1 ];
		if ( len_trim_whitespace_ > 0 ) std::memmove( c_str_, str_, len_trim_whitespace_ ); // Copy the string data
		c_str_[ len_trim_whitespace_ ] = '\0'; // Null-terminate
		return c_str_;
	}


	/// @brief Copy to a Pre-Allocated String
	Fstring::size_type
	Fstring::copy( cstring str, size_type const len_a, size_type const off ) const
	{
		assert( off <= len_ );
		size_type const len_copied( std::min( len_ - std::min( off, len_ ), len_a ) );
		if ( len_copied > 0 ) std::memmove( str, str_ + off, len_copied );
		return len_copied;
	}


	/// @brief Fstring == Fstring
	bool
	operator ==( Fstring const & s, Fstring const & t )
	{
		Fstring::size_type const min_len( std::min( s.len_, t.len_ ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			if ( s.str_[ i ] != t.str_[ i ] ) return false;
		}
		if ( s.len_ < t.len_ ) {
			for ( Fstring::size_type i = s.len_, e = t.len_; i < e; ++i ) {
				if ( t.str_[ i ] != SPACE ) return false;
			}
		} else if ( s.len_ > t.len_ ) {
			for ( Fstring::size_type i = t.len_, e = s.len_; i < e; ++i ) {
				if ( s.str_[ i ] != SPACE ) return false;
			}
		}
		return true;
	}


	/// @brief Fstring == string
	bool
	operator ==( Fstring const & s, std::string const & t )
	{
		Fstring::size_type const t_len( t.length() );
		Fstring::size_type const min_len( std::min( s.len_, t_len ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			if ( s.str_[ i ] != t[ i ] ) return false;
		}
		if ( s.len_ < t_len ) {
			for ( Fstring::size_type i = s.len_; i < t_len; ++i ) {
				if ( t[ i ] != SPACE ) return false;
			}
		} else if ( s.len_ > t_len ) {
			for ( Fstring::size_type i = t_len, e = s.len_; i < e; ++i ) {
				if ( s.str_[ i ] != SPACE ) return false;
			}
		}
		return true;
	}


	/// @brief Fstring == cstring
	bool
	operator ==( Fstring const & s, c_cstring const t )
	{
		Fstring::size_type const t_len( std::strlen( t ) );
		Fstring::size_type const min_len( std::min( s.len_, t_len ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			if ( s.str_[ i ] != t[ i ] ) return false;
		}
		if ( s.len_ < t_len ) {
			for ( Fstring::size_type i = s.len_; i < t_len; ++i ) {
				if ( t[ i ] != SPACE ) return false;
			}
		} else if ( s.len_ > t_len ) {
			for ( Fstring::size_type i = t_len, e = s.len_; i < e; ++i ) {
				if ( s.str_[ i ] != SPACE ) return false;
			}
		}
		return true;
	}


	/// @brief Fstring == char
	bool
	operator ==( Fstring const & s, char const c )
	{
		if ( s.empty() ) { // Zero-length Fstring
			return false;
		} else if ( s.str_[ 0 ] == c ) { // First character matches
			return ( s( 2 ).is_blank() ); // Rest is blank
		} else { // First character doesn't match
			return false;
		}
	}


	/// @brief Fstring <= Fstring
	bool
	operator <=( Fstring const & s, Fstring const & t )
	{
		Fstring::size_type const min_len( std::min( s.len_, t.len_ ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			unsigned char const s_i( s.str_[ i ] );
			unsigned char const t_i( t.str_[ i ] );
			if ( s_i < t_i ) {
				return true;
			} else if ( s_i > t_i ) {
				return false;
			}
		}
		if ( s.len_ < t.len_ ) {
			for ( Fstring::size_type i = s.len_, e = t.len_; i < e; ++i ) {
				unsigned char const t_i( t.str_[ i ] );
				if ( SPACE < t_i ) {
					return true;
				} else if ( SPACE > t_i ) {
					return false;
				}
			}
		} else if ( s.len_ > t.len_ ) {
			for ( Fstring::size_type i = t.len_, e = s.len_; i < e; ++i ) {
				unsigned char const s_i( s.str_[ i ] );
				if ( s_i < SPACE ) {
					return true;
				} else if ( s_i > SPACE ) {
					return false;
				}
			}
		}
		return true; // Equal
	}


	/// @brief Fstring < Fstring
	bool
	operator <( Fstring const & s, Fstring const & t )
	{
		Fstring::size_type const min_len( std::min( s.len_, t.len_ ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			unsigned char const s_i( s.str_[ i ] );
			unsigned char const t_i( t.str_[ i ] );
			if ( s_i < t_i ) {
				return true;
			} else if ( s_i > t_i ) {
				return false;
			}
		}
		if ( s.len_ < t.len_ ) {
			for ( Fstring::size_type i = s.len_, e = t.len_; i < e; ++i ) {
				unsigned char const t_i( t.str_[ i ] );
				if ( SPACE < t_i ) {
					return true;
				} else if ( SPACE > t_i ) {
					return false;
				}
			}
		} else if ( s.len_ > t.len_ ) {
			for ( Fstring::size_type i = t.len_, e = s.len_; i < e; ++i ) {
				unsigned char const s_i( s.str_[ i ] );
				if ( s_i < SPACE ) {
					return true;
				} else if ( s_i > SPACE ) {
					return false;
				}
			}
		}
		return false; // Equal
	}


	/// @brief Fstring <= string
	bool
	operator <=( Fstring const & s, std::string const & t )
	{
		Fstring::size_type const t_len( t.length() );
		Fstring::size_type const min_len( std::min( s.len_, t_len ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			unsigned char const s_i( s.str_[ i ] );
			unsigned char const t_i( t[ i ] );
			if ( s_i < t_i ) {
				return true;
			} else if ( s_i > t_i ) {
				return false;
			}
		}
		if ( s.len_ < t_len ) {
			for ( Fstring::size_type i = s.len_; i < t_len; ++i ) {
				unsigned char const t_i( t[ i ] );
				if ( SPACE < t_i ) {
					return true;
				} else if ( SPACE > t_i ) {
					return false;
				}
			}
		} else if ( s.len_ > t_len ) {
			for ( Fstring::size_type i = t_len, e = s.len_; i < e; ++i ) {
				unsigned char const s_i( s.str_[ i ] );
				if ( s_i < SPACE ) {
					return true;
				} else if ( s_i > SPACE ) {
					return false;
				}
			}
		}
		return true; // Equal
	}


	/// @brief Fstring < string
	bool
	operator <( Fstring const & s, std::string const & t )
	{
		Fstring::size_type const t_len( t.length() );
		Fstring::size_type const min_len( std::min( s.len_, t_len ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			unsigned char const s_i( s.str_[ i ] );
			unsigned char const t_i( t[ i ] );
			if ( s_i < t_i ) {
				return true;
			} else if ( s_i > t_i ) {
				return false;
			}
		}
		if ( s.len_ < t_len ) {
			for ( Fstring::size_type i = s.len_; i < t_len; ++i ) {
				unsigned char const t_i( t[ i ] );
				if ( SPACE < t_i ) {
					return true;
				} else if ( SPACE > t_i ) {
					return false;
				}
			}
		} else if ( s.len_ > t_len ) {
			for ( Fstring::size_type i = t_len, e = s.len_; i < e; ++i ) {
				unsigned char const s_i( s.str_[ i ] );
				if ( s_i < SPACE ) {
					return true;
				} else if ( s_i > SPACE ) {
					return false;
				}
			}
		}
		return false; // Equal
	}


	/// @brief Fstring <= cstring
	bool
	operator <=( Fstring const & s, c_cstring const t )
	{
		Fstring::size_type const t_len( std::strlen( t ) );
		Fstring::size_type const min_len( std::min( s.len_, t_len ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			unsigned char const s_i( s.str_[ i ] );
			unsigned char const t_i( t[ i ] );
			if ( s_i < t_i ) {
				return true;
			} else if ( s_i > t_i ) {
				return false;
			}
		}
		if ( s.len_ < t_len ) {
			for ( Fstring::size_type i = s.len_; i < t_len; ++i ) {
				unsigned char const t_i( t[ i ] );
				if ( SPACE < t_i ) {
					return true;
				} else if ( SPACE > t_i ) {
					return false;
				}
			}
		} else if ( s.len_ > t_len ) {
			for ( Fstring::size_type i = t_len, e = s.len_; i < e; ++i ) {
				unsigned char const s_i( s.str_[ i ] );
				if ( s_i < SPACE ) {
					return true;
				} else if ( s_i > SPACE ) {
					return false;
				}
			}
		}
		return true; // Equal
	}


	/// @brief Fstring < cstring
	bool
	operator <( Fstring const & s, c_cstring const t )
	{
		Fstring::size_type const t_len( std::strlen( t ) );
		Fstring::size_type const min_len( std::min( s.len_, t_len ) );
		for ( Fstring::size_type i = 0; i < min_len; ++i ) {
			unsigned char const s_i( s.str_[ i ] );
			unsigned char const t_i( t[ i ] );
			if ( s_i < t_i ) {
				return true;
			} else if ( s_i > t_i ) {
				return false;
			}
		}
		if ( s.len_ < t_len ) {
			for ( Fstring::size_type i = s.len_; i < t_len; ++i ) {
				unsigned char const t_i( t[ i ] );
				if ( SPACE < t_i ) {
					return true;
				} else if ( SPACE > t_i ) {
					return false;
				}
			}
		} else if ( s.len_ > t_len ) {
			for ( Fstring::size_type i = t_len, e = s.len_; i < e; ++i ) {
				unsigned char const s_i( s.str_[ i ] );
				if ( s_i < SPACE ) {
					return true;
				} else if ( s_i > SPACE ) {
					return false;
				}
			}
		}
		return false; // Equal
	}


	/// @brief Constant Substring: s( i, j )
	Fsubstring const
	Fstring::operator ()( size_type const i, size_type const j ) const
	{
		return Fsubstring( *this, i, j );
	}


	/// @brief Substring: s( i, j )
	Fsubstring
	Fstring::operator ()( size_type const i, size_type const j )
	{
		return Fsubstring( *this, i, j );
	}


	/// @brief Constant Tail Substring: s( i )
	Fstring const
	Fstring::operator ()( size_type const i ) const
	{
		return Fsubstring( *this, i );
	}


	/// @brief Tail Substring: s( i )
	Fsubstring
	Fstring::operator ()( size_type const i )
	{
		return Fsubstring( *this, i );
	}


	/// @brief Space-Free Head Constant Substring
	Fsubstring const
	Fstring::head() const
	{
		size_type const ie( find( SPACE ) );
		if ( ie == 0 ) {
			return Fsubstring( *this, 1, len_ );
		} else {
			return Fsubstring( *this, 1, ie - 1 );
		}
	}


	/// @brief Space-Free Head Substring
	Fsubstring
	Fstring::head()
	{
		size_type const ie( find( SPACE ) );
		if ( ie == 0 ) {
			return Fsubstring( *this, 1, len_ );
		} else {
			return Fsubstring( *this, 1, ie - 1 );
		}
	}


	/// @brief Space Tail Constant Substring
	Fsubstring const
	Fstring::tail() const
	{
		return Fsubstring( *this, len_trim() + 1 );
	}


	/// @brief Space Tail Substring
	Fsubstring
	Fstring::tail()
	{
		return Fsubstring( *this, len_trim() + 1 );
	}


	/// @brief Stream Input
	std::istream &
	operator >>( std::istream & stream, Fstring & s )
	{
		std::string ss;
		stream >> std::setw( s.len_ ) >> ss;
		s = ss;
		return stream;
	}


	/// @brief Get from Stream
	std::istream &
	get( std::istream & stream, Fstring & s )
	{
		if ( s.len_ > 0 ) {
			char * const buff( new char[ s.len_ + 1 ] );
			stream.get( buff, s.len_ + 1 ); // get adds null-terminator
			std::size_t const lb( std::strlen( buff ) );
			std::memcpy( s.str_, buff, lb );
			std::memset( s.str_ + lb, SPACE, s.len_ - lb );
			delete[] buff;
		}
		return stream;
	}


	/// @brief Get Line from Stream
	std::istream &
	getline( std::istream & stream, Fstring & s )
	{
		std::string ss;
		stream.width( s.len_ );
		std::getline( stream, ss );
		s = ss;
		return stream;
	}


	/// @brief Read from Stream
	std::istream &
	read( std::istream & stream, Fstring & s )
	{
		stream.read( s.str_, s.len_ );
		return stream;
	}


	/// @brief Read Available Characters from Stream
	std::istream &
	readsome( std::istream & stream, Fstring & s )
	{
		stream.readsome( s.str_, s.len_ );
		return stream;
	}


	/// @brief Stream Output
	std::ostream &
	operator <<( std::ostream & stream, Fstring const & s )
	{
		stream.write( s.str_, s.len_ );
		return stream;
	}


// Fstring


// Fsubstring


	/// @brief Copy Assignment
	Fsubstring &
	Fsubstring::operator =( Fsubstring const & s )
	{
		if ( this != &s ) {
			if ( len_ > s.len_ ) {
				if ( s.len_ > 0 ) std::memmove( str_, s.str_, s.len_ );
				std::memset( str_ + s.len_, SPACE, len_ - s.len_ ); // Space pad
			} else if ( len_ > 0 ) {
				std::memmove( str_, s.str_, len_ );
			}
		}
		return *this;
	}


	/// @brief = Fstring
	Fsubstring &
	Fsubstring::operator =( Fstring const & s )
	{
		if ( this != &s ) {
			if ( len_ > s.len_ ) {
				if ( s.len_ > 0 ) std::memmove( str_, s.str_, s.len_ );
				std::memset( str_ + s.len_, SPACE, len_ - s.len_ ); // Space pad
			} else if ( len_ > 0 ) {
				std::memmove( str_, s.str_, len_ );
			}
		}
		return *this;
	}


	/// @brief = string
	Fsubstring &
	Fsubstring::operator =( std::string const & s )
	{
		size_type const s_len( s.length() );
		if ( len_ > s_len ) {
			if ( s_len > 0 ) s.copy( str_, s_len );
			std::memset( str_ + s_len, SPACE, len_ - s_len ); // Space pad
		} else if ( len_ > 0 ) {
			s.copy( str_, len_ );
		}
		return *this;
	}


	/// @brief = cstring
	Fsubstring &
	Fsubstring::operator =( c_cstring const s )
	{
		size_type const s_len( std::strlen( s ) );
		if ( len_ > s_len ) {
			if ( s_len > 0 ) std::memmove( str_, s, s_len );
			std::memset( str_ + s_len, SPACE, len_ - s_len ); // Space pad
		} else if ( len_ > 0 ) {
			std::memmove( str_, s, len_ );
		}
		return *this;
	}


	/// @brief = char
	Fsubstring &
	Fsubstring::operator =( char const c )
	{
		if ( len_ > 0 ) {
			str_[ 0 ] = c;
			if ( len_ > 1 ) std::memset( str_ + 1, SPACE, len_ - 1 ); // Space pad
		}
		return *this;
	}


// Fsubstring


} // namespace ObjexxFCL
