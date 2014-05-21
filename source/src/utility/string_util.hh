// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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
#include <utility/stream_util.hh>
#include <utility/exit.hh>

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

/// @brief Reads the contents of <filename> into <contents>, preserving newline
/// characters. Aborts if an error is encoutered.
void ReadFromFileOrDie(const std::string& filename, std::string* contents);

/// @brief split given std::string using ' ' symbol.
utility::vector1< std::string > split(const std::string &s);

/// @brief split given std::string using whitespace as a separator.
/// Unlike string_split_multi_delim(), any group of mixed whitespace counts only as a single seperator.
utility::vector1< std::string > split_whitespace(const std::string &s);

///@brief combine strings with anything
std::string join(utility::vector1<std::string> const & s, std::string const & connector);

///@brief combine vector with anything
template<class T>
std::string join(utility::vector1<T> const & vector, std::string const & connector)
{
	std::ostringstream os;
	typename utility::vector1<T>::const_iterator begin= vector.begin();
	os << *begin++;
	for(; begin != vector.end(); ++begin){
		os<< connector<< *begin;
	}
	return os.str();
}

///@brief combine strings with anything
std::string join(std::vector<std::string> const & s, std::string const & connector);

/// @brief join space separations in a string with a connector such as '_'
std::string join(std::string const & string_w_spaces, std::string const & connector);

/// @brief split given std::string using ' ' symbol.
std::list< std::string > split_to_list(const std::string &s);

/// @brief split given std::string to a set using ' ' symbol.
std::set< std::string > split_to_set(std::string const & s);

/// @details split to vector1< std::string > using arbitrary split character
utility::vector1< std::string >
string_split( std::string const & in, char splitchar = ' ' );

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

void trim( std::string & s, const std::string & drop = " " );

/// @brief create a new string that drops all the unwanted substrings of
/// the original string.
std::string
trim( std::string const & s, std::string const & drop = " " );

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


///@brief find all environment variables with the form ${VARIABLE}
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
		if( s[c] == from ) s.replace(c,1,to);
	}
}

/// @brief find/replace strings within input string.
inline
std::string
replace_in( std::string const name_in, std::string const find_string, std::string const replace_string ){
	std::string name = name_in;
	platform::Size pos = name.find( find_string );
	while ( pos != std::string::npos ){
		name = name.replace( pos, find_string.size(), replace_string );
		pos = name.find( find_string );
	}
	return name;
}


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
	if (in.size()==0) {
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
make_tag_with_dashes( utility::vector1< int > res_vector );

// Compactifies vectors of ints and chars (resnum and chain):  1A 2A 3A 9B 10B 11B to "A:1-3 B:9-11"
std::string
make_tag_with_dashes( utility::vector1< int > res_vector,
											utility::vector1< char > chain_vector );

std::string
make_tag( utility::vector1< int > res_vector );

/// @brief  converts string like "1-3 20-22" or "A:1-5 B:20-22" to vectors containing resnums and chains.
std::pair< std::vector< int >, std::vector< char > >
get_resnum_and_chain( std::string const & s, bool & string_is_ok );

/// @brief helper function for get_resnum_and_chain
bool
get_resnum_and_chain_from_one_tag( std::string const & tag,
																	 std::vector< int > & resnum,
																	 std::vector< char > & chains );


platform::Size
get_num_digits( platform::Size value);

}  // namespace utility

#endif  // INCLUDED_utility_string_util_HH
