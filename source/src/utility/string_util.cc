// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/string_util.cc
///
/// @brief  Some std::string helper functions.
/// @author Sergey Lyskov
#include <platform/types.hh>

#include <utility/numbers.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <locale>

// C/C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/uuid/sha1.hpp>

namespace utility {

void ReadFromFileOrDie(const std::string& filename, std::string* contents) {
	using std::ifstream;
	using std::string;
	using std::stringstream;

	assert(contents);

	ifstream in(filename.c_str());
	if ( !in ) {
		stringstream ss;
		ss << "The specified file " << filename
			<< "does not exist or lacks sufficient permissions";
		utility_exit_with_message(ss.str());
	}

	string line;
	while ( in.good() ) {
		getline(in, line);
		(*contents) +=line + "\n";
	}
	in.close();
}

utility::vector1< std::string > split(const std::string &s)
{
	utility::vector1<std::string> r;
	unsigned int start=0, i=0;
	while ( start < s.size() ) {
		if ( s[i] == ' ' /*|| i==s.size()-1 */ ) {
			std::string add(s.begin()+start, s.begin()+i);
			if ( add.size() != 0 ) r.push_back( add );
			start = i+1;
		}
		i++;
		if ( i == s.size() ) {
			std::string add(s.begin()+start, s.begin()+i);
			if ( add.size() != 0 ) r.push_back( add );
			break;
		}
	}
	return r;
}

utility::vector1< std::string > split_whitespace(const std::string &s){
	std::istringstream ss( s );
	utility::vector1<std::string> r;
	while ( !ss.fail() ) {
		std::string e;
		ss >> e;
		if ( ss.fail() ) break;
		r.push_back(e);
	}
	return r;
}

std::string join(utility::vector1<std::string> const & s, std::string const & connector){
	std::ostringstream os;
	utility::vector1<std::string>::const_iterator begin= s.begin();
	os << *begin++;
	for ( ; begin != s.end(); ++begin ) {
		os<< connector<< *begin;
	}
	return os.str();
}

std::string join(std::vector<std::string> const & s, std::string const & connector){
	std::ostringstream os;
	utility::vector1<std::string>::const_iterator begin= s.begin();
	os << *begin++;
	for ( ; begin != s.end(); ++begin ) {
		os<< connector<< *begin;
	}
	return os.str();
}

std::string replace_spaces(std::string const & string_w_spaces, std::string const & replacement){
	//std::string trimmed= trim(string_w_spaces);
	utility::vector1<std::string> pieces= split(string_w_spaces);
	return join(pieces, replacement);
}

/// @details split given std::string using ' ' symbol.
std::list< std::string > split_to_list(const std::string &s) {
	std::list<std::string> r;
	unsigned int start=0, i=0;
	while ( start < s.size() ) {
		if ( s[i] == ' ' /*|| i==s.size()-1 */ ) {
			std::string add(s.begin()+start, s.begin()+i);
			if ( add.size() != 0 ) r.push_back( add );
			start = i+1;
		}
		i++;
		if ( i == s.size() ) {
			std::string add(s.begin()+start, s.begin()+i);
			if ( add.size() != 0 ) r.push_back( add );
			break;
		}
	}
	return r;
}

/// @details split given std::string using ' ' symbol.
std::set< std::string > split_to_set(const std::string &s) {
	std::set<std::string> r;
	unsigned int start=0, i=0;
	while ( start < s.size() ) {
		if ( s[i] == ' ' /*|| i==s.size()-1 */ ) {
			std::string add(s.begin()+start, s.begin()+i);
			if ( add.size() != 0 ) r.insert( add );
			start = i+1;
		}
		i++;
		if ( i == s.size() ) {
			std::string add(s.begin()+start, s.begin()+i);
			if ( add.size() != 0 ) r.insert( add );
			break;
		}
	}
	return r;
}

utility::vector1< std::string >
string_split( std::string const & in, char splitchar /* = ' ' */ )
{
	utility::vector1< std::string > parts;
	size_t i(0), j(0);
	while ( j != std::string::npos ) {
		j = in.find( splitchar, i );
		std::string const part = in.substr(i,j-i);
		parts.push_back( part );
		i = j+1;
	}
	return parts;
}

/// @details split to vector1< std::string > using arbitrary split character, but no empty strings (closer to python string::split)
utility::vector1< std::string >
string_split_simple( std::string const & in, char splitchar /* = ' ' */ )
{
	utility::vector1< std::string > parts;
	size_t i(0), j(0);
	while ( j != std::string::npos ) {
		j = in.find( splitchar, i );
		std::string const part = in.substr(i,j-i);
		if ( part.size() > 0 ) {
			parts.push_back( part );
		}
		i = j+1;
	}
	return parts;
}


//overloaded to split on any of an array of chars, useful to split on any whitespace
utility::vector1< std::string >
string_split_multi_delim( std::string const & in, std::string splitchars )
{
	utility::vector1< std::string > parts;
	size_t i(0), j(0);
	while ( j != std::string::npos ) {
		//find first instance of any of the splitchars
		j = in.find_first_of( splitchars, i );
		parts.push_back( in.substr(i,j-i) );
		i = j+1;
	}
	return parts;
}


/// @details convert a string to a float
float string2float( std::string st ){
	float i;
	std::stringstream ss( st );
	ss >> i;
	if ( !ss ) {
		return -1;
	}
	return i;
}

/// @details convert a string to an int
int string2int( std::string st ){
	int i;
	std::stringstream ss( st );
	ss >> i;
	if ( !ss ) {
		return -1;
	}
	return i;
}

/// @brief convert a string to a Size, returns numeric::get_undefined_size() on failure
platform::Size string2Size( std::string st ){
	platform::Size i;
	std::stringstream ss( st );
	ss >> i;
	if ( !ss ) {
		return get_undefined_size();
	}
	return i;
}

/// @brief convert a string to a Real, returns numeric::get_undefined_real() on failure
platform::Real string2Real( std::string st ){
	platform::Real i;
	std::stringstream ss( st );
	ss >> i;
	if ( !ss ) {
		return get_undefined_real();
	}
	return i;
}

// @brief Reads an unsigned int from string <x>, writing the result
// to output parameter <y>, which must be non-NULL. The result is
// undefined if the input string is malformed.
void string2uint(const std::string& x, unsigned int* y) {
	assert(y != NULL);
	std::stringstream ss(x);
	ss >> *y;
}

/// @details compares two strings ignoring leading and trailing spaces
bool trimmed_compare( std::string const & s1, std::string const & s2 )
{
	std::string const space( " " );

	std::size_t const s1_start( s1.find_first_not_of( space ) );
	std::size_t const s1_end( s1.find_last_not_of( space ) );
	std::size_t const s2_start( s2.find_first_not_of( space ) );
	std::size_t const s2_end( s2.find_last_not_of( space ) );

	std::size_t const s1_len( s1_end - s1_start + 1 );
	std::size_t const s2_len( s2_end - s2_start + 1 );

	return ( ( s1_len == s2_len ) && ( s1.compare( s1_start, s1_len, s2, s2_start, s2_len ) == 0 ) );
}

bool startswith(std::string const & haystack, std::string const & needle)
{
	if ( haystack.length() < needle.length() ) return false;
	else return ( haystack.compare(0, needle.length(), needle) == 0 );
}

bool endswith(std::string const & haystack, std::string const & needle)
{
	if ( haystack.length() < needle.length() ) return false;
	else return ( haystack.compare(haystack.size()-needle.size(),needle.size(),needle) == 0 );
}

void slurp(std::istream & in, std::string & out)
{
	std::string line;
	std::ostringstream os;
	while ( std::getline(in,line) ) {
		os << line << std::endl;
	}
	out.append( os.str());
}

void trim( std::string & s, const std::string & drop)
{
	std::string r = s.erase( s.find_last_not_of(drop)+1 );
	r.erase( 0, r.find_first_not_of(drop) );
	s = r;
}

std::string
trim( std::string const & s, std::string const & drop )
{
	std::string copystr( s );
	trim( copystr, drop );
	return copystr;
}

void add_spaces_left_align( std::string & st, std::size_t const newlen )
{
	std::size_t const to_add = newlen - st.length();
	if ( to_add > 0 ) {
		std::string st_to_add("");
		st_to_add.append(to_add,' ');
		st = st + st_to_add;
	}
}

void add_spaces_right_align( std::string & st, std::size_t const newlen )
{
	std::size_t const to_add = newlen - st.length();
	if ( to_add > 0 ) {
		std::string st_to_add("");
		st_to_add.append(to_add,' ');
		st = st_to_add + st;
	}
}

bool is_string_numeric(std::string const & input)
{
	std::locale loc;
	for ( platform::Size i = 0 ; i < input.size(); ++i ) {
		char current = input[i];
		if ( std::isdigit(current,loc) || current == '-' || current == '+' || current =='E' ||current=='e' ) {
			continue;
		} else {
			return false;
		}
	}
	return true;
}

std::string
file_contents( std::string const & file_name )
{
	vector1< std::string > text;
	std::string line;
	io::izstream textstream( file_name );
	if ( ! textstream ) {
		throw excn::EXCN_Msg_Exception( "Could not open file " + file_name  );
	}
	int strsize( 0 );
	while ( getline(textstream, line) ) {
		text.push_back(line + "\n");
		strsize += line.size() + 1;
	}
	textstream.close();

	std::string alltext;
	alltext.reserve( strsize );
	for ( unsigned int ii = 1; ii <= text.size(); ++ ii ) {
		alltext += text[ii];
	}
	return alltext;
}

std::string file_basename(const std::string& full_path) {
	return filename(full_path);
}

std::string filename(const std::string& path) {
	utility::file::FileName f(path);
	return f.base() + f.extension();
}

std::string pathname(const std::string& path) {
	return utility::file::FileName(path).path();
}

std::string replace_environment_variables(std::string input)
{
	const std::string start("${");
	const std::string end("}");

	platform::Size start_position = 0;
	while ( true )
			{
		start_position = input.find(start);
		if ( start_position != std::string::npos ) {
			platform::Size end_position = input.find(end,start_position);
			if ( start_position == std::string::npos ) {
				utility_exit_with_message("opening ${ but no closing } around an environment variable, check your options file");
			}

			platform::Size env_length = end_position-start_position;

			std::string env_name = input.substr(start_position+2,env_length-2);
			char * env_value = getenv(env_name.c_str());
			if ( !env_value ) {
				utility_exit_with_message("environment variable "+env_name+" does not exist");
			}

			input.replace(start_position, env_length+1,env_value);

		} else {
			return input;
		}
	}
}

std::string string_to_sha1(std::string const & input_string)
{
	//Based on https://gist.github.com/990731
	unsigned int digest[5];
	boost::uuids::detail::sha1 hasher;
	hasher.process_bytes(input_string.c_str(),input_string.size());


	char hash[20];
	std::stringstream output_hash;

	hasher.get_digest(digest);

	for ( int i = 0; i < 5; ++i ) {
		const char* tmp = reinterpret_cast<char*>(digest);
		hash[i*4] = tmp[i*4+3];
		hash[i*4+1] = tmp[i*4+2];
		hash[i*4+2] = tmp[i*4+1];
		hash[i*4+3] = tmp[i*4];
	}

	output_hash << std::hex;

	for ( int i = 0; i < 20 ; ++i ) {
		output_hash << ((hash[i] & 0x000000F0) >> 4) <<  (hash[i] & 0x0000000F);
	}

	return output_hash.str();
}

////////////////////////////////////////////////////////////////////////////
// Compactifies vectors of ints:  1 2 3 9 10 11 to "1-3 9-11"
// The function to go the other way (from string to vector) is available in
// ObjexxFCL::get_ints(), I think.
//
std::string
make_tag_with_dashes( utility::vector1< int > res_vector,
	char const delimiter /* = ' ' */){
	utility::vector1< char > chains;
	for ( platform::Size n = 0; n < res_vector.size(); n++ ) chains.push_back( ' ' );
	return make_tag_with_dashes( res_vector, chains, delimiter );
}

////////////////////////////////////////////////////////////////////////////
// Compactifies vectors of ints and chars (resnum and chain):  1A 2A 3A 9B 10B 11B to "A:1-3 B:9-11"
// The function to go the other way (from string to two vectors) is available below in get_resnum_and_chain()
//
std::string
make_tag_with_dashes( utility::vector1< int > res_vector,
	utility::vector1< char > chain_vector,
	char const delimiter /* = ' ' */){

	using namespace ObjexxFCL;
	std::string tag = "";

	if ( res_vector.size() == 0 ) return tag;
	runtime_assert( res_vector.size() == chain_vector.size() );

	int start_segment = res_vector[1];
	int last_res = res_vector[1];
	int last_chain = chain_vector[1];
	utility::vector1< std::pair<int,int> > res_vector_segments;
	utility::vector1< char > chains_for_segments;
	for ( platform::Size n = 2; n<= res_vector.size(); n++ ) {
		if ( res_vector[n] != last_res+1  || chain_vector[n] != last_chain ) {
			res_vector_segments.push_back( std::make_pair( start_segment, last_res ) );
			chains_for_segments.push_back( last_chain );
			start_segment = res_vector[n];
		}
		last_res = res_vector[n];
		last_chain = chain_vector[n];
	}
	res_vector_segments.push_back( std::make_pair( start_segment, last_res ) );
	chains_for_segments.push_back( last_chain );

	for ( platform::Size n = 1; n <= res_vector_segments.size(); n++ ) {
		if ( n > 1 ) tag += delimiter;
		std::pair< int, int > const & segment = res_vector_segments[n];
		if ( chains_for_segments[n] != '\0' &&
				chains_for_segments[n] != ' '  &&
				chains_for_segments[n] != '_' ) tag += std::string(1,chains_for_segments[n]) + ":";
		if ( segment.first == segment.second ) {
			tag += string_of( segment.first );
		} else {
			tag += string_of( segment.first )+"-"+string_of(segment.second);
		}
	}

	return tag;
}


std::string
make_segtag_with_dashes( utility::vector1< int > res_vector,
	utility::vector1< std::string > segid_vector,
	char const delimiter /* = ' ' */){

	using namespace ObjexxFCL;
	std::string tag = "";

	if ( res_vector.size() == 0 ) return tag;
	runtime_assert( res_vector.size() == segid_vector.size() );

	int start_segment = res_vector[1];
	int last_res = res_vector[1];
	std::string last_segid = segid_vector[1];
	utility::vector1< std::pair<int,int> > res_vector_segments;
	utility::vector1< std::string > segids_for_segments;
	for ( platform::Size n = 2; n<= res_vector.size(); n++ ) {
		if ( res_vector[n] != last_res+1  || segid_vector[n] != last_segid ) {
			res_vector_segments.push_back( std::make_pair( start_segment, last_res ) );
			segids_for_segments.push_back( last_segid );
			start_segment = res_vector[n];
		}
		last_res = res_vector[n];
		last_segid = segid_vector[n];
	}
	res_vector_segments.push_back( std::make_pair( start_segment, last_res ) );
	segids_for_segments.push_back( last_segid );

	for ( platform::Size n = 1; n <= res_vector_segments.size(); n++ ) {
		if ( n > 1 ) tag += delimiter;
		std::pair< int, int > const & segment = res_vector_segments[n];
		if ( segids_for_segments[n] != "" &&
				segids_for_segments[n] != "    " ) {
			tag += segids_for_segments[n] + ":";
		}
		if ( segment.first == segment.second ) {
			tag += string_of( segment.first );
		} else {
			tag += string_of( segment.first )+"-"+string_of(segment.second);
		}
	}

	return tag;
}

/////////////////////////////////////////////////////////////////////////////////
std::string
make_tag( utility::vector1< int > res_vector ){

	using namespace ObjexxFCL;
	std::string tag = "";

	for ( platform::Size n = 1; n <= res_vector.size(); ++n ) {
		if ( n > 1 ) tag += " ";
		tag += string_of( res_vector[n] );
	}

	return tag;
}

/// @brief  converts string like "1-3 20-22" or "A:1-5 B:20-22" to vectors containing resnums and chains.
/// @author rhiju
//
//  #detailed  several kinds of tags are OK, including "A:1-5 B:20-22" and "A1-5 B20,21,22".
//
std::pair< std::vector< int >, std::vector< char > >
get_resnum_and_chain( std::string const & s, bool & string_is_ok ){

	string_is_ok = true;
	std::vector< int  > resnum;
	std::vector< char > chain;

	std::string s_nocommas = replace_in( s, ",", " " ); // order of operations issue?
	utility::vector1< std::string > const tags = utility::string_split( s_nocommas );
	for ( platform::Size n = 1; n <= tags.size(); n++ ) {
		string_is_ok = get_resnum_and_chain_from_one_tag( tags[n], resnum, chain );
		if ( !string_is_ok ) break;
	}
	return std::make_pair( resnum,chain );
}

std::pair< std::vector< int >, std::vector< std::string > >
get_resnum_and_segid( std::string const & s, bool & string_is_ok ){

	string_is_ok = true;
	std::vector< int  > resnum;
	std::vector< std::string > segid;

	std::string s_nocommas = replace_in( s, ",", " " ); // order of operations issue?
	// Each tag looks like
	// abcd:n-m, bcd could be spaces...
	// Strategy: add substring from start to colon (skips spaces!)
	// then from the space on forward.
	// utility::string_split( s_nocommas );
	utility::vector1< std::string > tags;

	// populate tags
	while ( true ) {
		std::string tag = s_nocommas.substr( 0, s_nocommas.find_first_of( ':' )+1 );
		s_nocommas = s_nocommas.substr( s_nocommas.find_first_of( ':' ) +1 );

		if ( s_nocommas.find_first_of( ' ' ) == std::string::npos ) {
			tag += s_nocommas;
			tags.push_back( tag );
			break;
		}

		tag += s_nocommas.substr( s_nocommas.find_first_of( ' ' ) );
		tags.push_back( tag );
		s_nocommas = s_nocommas.substr( s_nocommas.find_first_of( ' ' ) +1 );
	}

	for ( platform::Size n = 1; n <= tags.size(); n++ ) {
		string_is_ok = get_resnum_and_segid_from_one_tag( tags[n], resnum, segid );
		if ( !string_is_ok ) break;
	}
	return std::make_pair( resnum, segid );
}

/// @brief helper function for get_resnum_and_chain
bool
get_resnum_and_chain_from_one_tag( std::string const & tag,
	std::vector< int > & resnum,
	std::vector< char > & chains ){
	bool string_is_ok( false );
	std::vector< int > resnum_from_tag;
	char chain( ' ' );
	std::string const numerical("-0123456789");

	size_t found_colon = tag.find( ":" );
	if ( found_colon == std::string::npos ) {
		if ( numerical.find( tag[0] ) == std::string::npos ) { // looks like a chain character at beginning
			chain= tag[0];
			resnum_from_tag = ObjexxFCL::ints_of( tag.substr(1), string_is_ok );
		} else {
			resnum_from_tag = ObjexxFCL::ints_of( tag, string_is_ok );
		}
	} else {
		if ( found_colon != 1 ) return false;
		chain = tag[0];
		resnum_from_tag = ObjexxFCL::ints_of( tag.substr(2), string_is_ok );
	}

	for ( platform::Size n = 0; n < resnum_from_tag.size(); n++ ) {
		resnum.push_back( resnum_from_tag[ n ] );
		chains.push_back( chain );
	}

	return string_is_ok;
}

bool
get_resnum_and_segid_from_one_tag( std::string const & tag,
	std::vector< int > & resnum,
	std::vector< std::string > & segids ){
	bool string_is_ok( false );
	std::vector< int > resnum_from_tag;
	std::string segid( "" );

	size_t found_colon = tag.find( ":" );

	if ( found_colon != 1 ) return false;

	segid = tag.substr(0,4);
	resnum_from_tag = ObjexxFCL::ints_of( tag.substr(4), string_is_ok );

	for ( platform::Size n = 0; n < resnum_from_tag.size(); n++ ) {
		resnum.push_back( resnum_from_tag[ n ] );
		segids.push_back( segid );
	}

	return string_is_ok;
}

platform::Size
get_num_digits( platform::Size value){
	return (platform::Size)ceil( std::log10(value + 1.) );
}

std::string
replace_in( std::string const & name_in, std::string const & find_string, std::string const & replace_string ){
	std::string name = name_in;
	// WARNING WARNING WARNING: Do not change pos back to platform::Size. In general platform::Size type
	//                          is not the same as size_t and algorithm below is sensitive to this
	// platform::Size pos = name.find( find_string );
	size_t pos = name.find( find_string );
	while ( pos != std::string::npos ) {
		name = name.replace( pos, find_string.size(), replace_string );

		//pos = name.find( find_string );

		// No, need to look for LATER instances of the replacement string.
		pos = name.find( find_string, pos + replace_string.size() );
	}
	return name;
}


// @brief Return a copy of the string with leading and trailing characters removed
std::string strip(std::string const & source, char c)
{
	std::string::const_iterator begin = source.begin();
	std::string::const_iterator end = source.end();
	for ( std::string::const_iterator p = source.begin(); p!=source.end(); ++p ) {
		if ( *p == c ) begin = p+1;
		else break;
	}

	for ( std::string::const_iterator p = source.end(); p!=begin; --p ) {
		if ( *(p-1) == c ) end = p-1;
		else break;
	}

	return std::string(begin, end);
}


} // namespace utility
