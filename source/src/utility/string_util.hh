// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/string_util.hh
///
/// @brief  Some std::string helper functions.
/// @author Sergey Lyskov
#ifndef INCLUDED_utility_string_util_hh
#define INCLUDED_utility_string_util_hh

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/stream_util.hh>

// Boost headers
#include <boost/algorithm/string/erase.hpp>

// C++ headers
#include <list>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>

namespace utility {

//These are useful string utilities from the ObjexxFCL namespace - transclude them here so you have a one-stop shop for string functions.

using ObjexxFCL::uppercase;
using ObjexxFCL::lowercase;
using ObjexxFCL::strip_whitespace;
using ObjexxFCL::uppercased;
using ObjexxFCL::lowercased;
using ObjexxFCL::stripped_whitespace;


/// @brief Reads the contents of <filename> into <contents>, preserving newline
/// characters. Aborts if an error is encoutered.
void ReadFromFileOrDie(const std::string& filename, std::string* contents);

/// @brief split given std::string using ' ' symbol.
utility::vector1< std::string > split(const std::string &s);

/// @brief split given std::string using whitespace as a separator.
/// Unlike string_split_multi_delim(), any group of mixed whitespace counts only as a single seperator.
utility::vector1< std::string > split_whitespace(const std::string &s);

/// @brief combine strings with anything
std::string join(utility::vector1<std::string> const & s, std::string const & connector);

/// @brief combine vector with anything
template<class T>
std::string join(utility::vector1<T> const & vector, std::string const & connector)
{
	std::ostringstream os;
	typename utility::vector1<T>::const_iterator begin= vector.begin();
	os << *begin++;
	for ( ; begin != vector.end(); ++begin ) {
		os<< connector<< *begin;
	}
	return os.str();
}

/// @brief combine strings with anything
std::string join(std::vector<std::string> const & s, std::string const & connector);

/// @brief Join vector of strings in to single string
template< platform::SSize L>
std::string join(vectorL<L, std::string> const & s, std::string const & connector){
	std::ostringstream os;
	typename utility::vectorL<L, std::string>::const_iterator begin= s.begin();
	os << *begin++;
	for ( ; begin != s.end(); ++begin ) {
		os<< connector<< *begin;
	}
	return os.str();
}


/// @brief replace space separations in a string with a connector such as '_'
std::string replace_spaces(std::string const & string_w_spaces, std::string const & replacement);

/// @brief split given std::string using ' ' symbol.
std::list< std::string > split_to_list(const std::string &s);

/// @brief split given std::string to a set using ' ' symbol.
std::set< std::string > split_to_set(std::string const & s);

/// @details split to vector1< std::string > using arbitrary split character
utility::vector1< std::string >
string_split( std::string const & in, char splitchar = ' ' );

/// @brief split to vector1< std::string > using arbitrary split character, but no empty strings (closer to python string::split)
utility::vector1< std::string >
string_split_simple( std::string const & in, char splitchar = ' ' );

/// @details split to vector< std::string > using any of arbitrary split characters
utility::vector1< std::string >
string_split_multi_delim( std::string const & in, std::string splitchars = " \t" );

/// @brief convert a string to a float, returns -1 on failure
float string2float( std::string st );

/// @brief convert a string to an int, returns -1 on failure
int string2int( std::string st );

/// @brief convert a string to a Size, returns numeric::get_undefined_size() on failure
platform::Size string2Size( std::string st );

/// @brief convert a string to a Real, returns numeric::get_undefined_real() on failure
platform::Real string2Real( std::string st );

// @brief Reads an unsigned int from string <x>, writing the result
// to output parameter <y>, which must be non-NULL. If the read was not
// successful, this function call has no effect on the value of <y> that
// was present prior to invokation.
void string2uint(const std::string& x, unsigned int *y);

/// @brief True iff haystack starts with needle
bool startswith(std::string const & haystack, std::string const & needle);

/// @brief True iff haystack ends with needle
bool endswith(std::string const & haystack, std::string const & needle);

void slurp(std::istream & in, std::string & out);

/// @brief Remove any charachters in drop from the front and back of the string
/// Use strip() for the value return version
void trim( std::string & s, const std::string & drop = " " );

/// @brief Return a copy of the string with leading and trailing characters removed
std::string strip(std::string const & source, char c=' ');

/// @brief Return a copy of the string with leading and trailing characters removed
/// Any charachters in drop will be removed
/// For the in place version, see trim()
std::string strip(std::string const & source, std::string const & drop);

/// @brief Ambiguious with the trim( std::string & s ) -- Deprecated:
/// use strip() instead for return-value trimming
inline
std::string
trim( std::string const & s, std::string const & drop = " " ) {
	return strip( s, drop );
}

/// @brief compares two strings ignoring leading and trailing spaces
bool trimmed_compare( std::string const & s1, std::string const & s2 );

/// @brief adds spaces to a left aligned string until a given length is reached
void add_spaces_left_align( std::string & st, std::size_t const newlen );

/// @brief adds spaces to a right aligned string until a given length is reached
void add_spaces_right_align( std::string & st, std::size_t const newlen );

// @brief return true of the string has only [0-9], ,'+','-','.' or '[Ee]'
bool is_string_numeric(std::string const & input);

/// @brief Read the entire contents of a file into a string.  All end-of-line characters are replaced
/// by "\n".  Throws a utility::excn::EXCN_msg_exception if the file cannot be opened.
std::string file_contents( std::string const & file_name );

std::string file_basename( std::string const & full_path );

// "/foo/bar/baz" => "baz"
// "/foo/bar/baz.cc" => "baz.cc"
std::string filename(const std::string& path);

// "/foo/bar/baz" => "/foo/bar/"
std::string pathname(const std::string& path);


/// @brief find all environment variables with the form ${VARIABLE}
/// and replace with the contents of that environment variable.
/// if the environment variable does not exist, return string::npos
std::string replace_environment_variables(std::string input);

/// @brief Compares two strings, ignoring spaces.  Useful for comparing atom
/// name strings which have pdb-alignment built into them.  Slightly dangerous
/// if you consider the fact that atom names in the PDB are different for
/// different indentation rules: ' CA ' is c-alpha.  'CA  ' is calcium.
inline
bool same_ignoring_spaces( std::string const & s1, std::string const & s2 ) {
	std::string t1 = boost::algorithm::erase_all_copy(s1, " ");
	std::string t2 = boost::algorithm::erase_all_copy(s2, " ");
	return t1 == t2;
}

//@brief compute the sha1 hash of a string and return it as a string in hexadecimal form
std::string string_to_sha1(std::string const & input_string);

inline
void replace_in( std::string & s, const char from, const char *to )
{
	// fix string
	for ( unsigned int c = 0; c < s.length(); ++c ) {
		if ( s[c] == from ) s.replace(c,1,to);
	}
}

/// @brief Generate new string from 'source' by replacing all occurrences of 'from' to 'to' string.
std::string replace_in( std::string const & source, std::string const & from, std::string const & to );


template <class T>
inline std::string to_string (const T & t)
{
	std::ostringstream ss;
	ss << t;
	return ss.str();
}

template <class T>
inline T const from_string (std::string const & s, T )
{
	T t;
	std::istringstream ss(s);
	ss >> t;
	if ( ss.fail() ) {
		const char* type = typeid(T).name();
		utility_exit_with_message("cannot convert string "+s+" to type "+type);
	}

	return t;
}

template <class T>
inline utility::vector1<T> const string_split (std::string const &in,char splitchar,T)
{
	utility::vector1<T> parts;
	if ( in.size()==0 ) {
		return parts;
	}

	size_t i(0), j(0);
	while ( j != std::string::npos ) {
		j = in.find( splitchar, i );
		std::string item = in.substr(i,j-i);
		T t;
		std::istringstream ss(item);
		ss >> t;
		if ( ss.fail() ) {
			const char* type = typeid(T).name();
			utility_exit_with_message("cannot convert string '"+item+"' to type "+type);
		}

		parts.push_back( t );
		i = j+1;
	}
	return parts;
}

/// @brief String accepted as a true value?
bool inline
is_true_string( std::string const & value_str )
{
	return (
		( value_str == "true" ) ||
		( value_str == "True" ) ||
		( value_str == "TRUE" ) ||
		( value_str == "t" ) ||
		( value_str == "T" ) ||
		( value_str == "1" ) ||
		( value_str == "on" ) ||
		( value_str == "On" ) ||
		( value_str == "ON" ) ||
		( value_str == "y" ) ||
		( value_str == "Y" ) ||
		( value_str == "yes" ) ||
		( value_str == "Yes" ) ||
		( value_str == "YES" ) );
}

/// @brief String accepted as a false value?
bool inline
is_false_string( std::string const & value_str )
{
	return (
		( value_str == "false" ) ||
		( value_str == "False" ) ||
		( value_str == "FALSE" ) ||
		( value_str == "f" ) ||
		( value_str == "F" ) ||
		( value_str == "0" ) ||
		( value_str == "off" ) ||
		( value_str == "Off" ) ||
		( value_str == "OFF" ) ||
		( value_str == "n" ) ||
		( value_str == "N" ) ||
		( value_str == "no" ) ||
		( value_str == "No" ) ||
		( value_str == "NO" ) );
}

/// @brief Compactifies vectors of ints:  1 2 3 9 10 11 to "1-3 9-11"
std::string
make_tag_with_dashes( utility::vector1< int > res_vector,
	char const delimiter = ' ' );

// Compactifies vectors of ints and chars (resnum and chain):  1A 2A 3A 9B 10B 11B to "A:1-3 B:9-11"
std::string
make_tag_with_dashes( utility::vector1< int > res_vector,
	utility::vector1< char > chain_vector,
	char const delimiter = ' ' );

std::string
make_segtag_with_dashes( utility::vector1< int > res_vector,
	utility::vector1< std::string > segid_vector,
	char const delimiter = ' ');

std::string
make_tag( utility::vector1< int > res_vector );

/// @brief  converts string like "1-3 20-22" or "A:1-5 B:20-22" to vectors containing resnums and chains.
std::pair< std::vector< int >, std::vector< char > >
get_resnum_and_chain( std::string const & s, bool & string_is_ok );

/// @brief  converts string like "1-3 20-22" or "A:1-5 B:20-22" to vectors containing resnums and chains.
std::pair< std::vector< int >, std::vector< char > >
get_resnum_and_chain( std::string const & s );

/// @brief helper function for get_resnum_and_chain
bool
get_resnum_and_chain_from_one_tag( std::string const & tag,
	std::vector< int > & resnum,
	std::vector< char > & chains );

/// @brief  converts string like "1-3 20-22" or "A:1-5 B:20-22" to vectors containing resnums and chains.
std::pair< std::vector< int >, std::vector< std::string > >
get_resnum_and_segid( std::string const & s, bool & string_is_ok );

/// @brief helper function for get_resnum_and_chain
bool
get_resnum_and_segid_from_one_tag( std::string const & tag,
	std::vector< int > & resnum,
	std::vector< std::string > & chains );

platform::Size
get_num_digits( platform::Size value);

}  // namespace utility

#endif  // INCLUDED_utility_string_util_HH
