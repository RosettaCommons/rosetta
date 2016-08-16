// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/JobOptions.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <basic/resource_manager/JobOptions.hh>

// Platform Headers
#include <platform/types.hh>
#include <utility/vector1.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>


//C++ Headers
#include <string>
#include <map>
#include <iomanip>

namespace basic {
namespace resource_manager {

JobOptions::JobOptions() : track_insertion_order_( false ) {}

/// @details Auto-generated virtual destructor
JobOptions::~JobOptions() {}

using platform::Real;
using utility::options::BooleanOptionKey;
using utility::options::BooleanVectorOptionKey;
using utility::options::FileOptionKey;
using utility::options::FileVectorOptionKey;
using utility::options::IntegerOptionKey;
using utility::options::IntegerVectorOptionKey;
using utility::options::PathOptionKey;
using utility::options::PathVectorOptionKey;
using utility::options::RealOptionKey;
using utility::options::RealVectorOptionKey;
using utility::options::StringOptionKey;
using utility::options::StringVectorOptionKey;
using utility::excn::EXCN_Msg_Exception;
using utility::file::FileName;
using utility::file::PathName;
using utility::vector1;
using std::map;
using std::string;
using std::endl;
using std::setw;

using namespace utility::options;

/// @brief simple enough function that calls the bit-shift operator in the input stream; why bother?
/// For the sake of writing out boolean values as "true" or "false"
template< class S >
void
write_type( std::ostream & out, S val )
{
	out << val;
}

/// @brief Template specialization to write boolean values as "true" or "false"
template <>
void
write_type< bool > ( std::ostream & out, bool val )
{
	out << ( val ? "true" : "false" );
}

template< class T, class S >
void
show_option_map(
	std::map< T, S > const & option_map,
	std::ostream & out,
	char const * option_name
)
{
	if ( option_map.empty() ) return;
	for (
			typename map< T, S >::const_iterator
			o = option_map.begin(), oe = option_map.end(); o != oe; ++o ) {
		out
			<< std::setiosflags(std::ios::left) << setw(16) << option_name
			<< std::setiosflags(std::ios::left) << setw(16) << o->first.id();
		write_type( out, o->second );
		out << endl;
	}
}

template< class T, class S >
void
show_option_vector_map(
	std::map< T, vector1< S > > const & option_map,
	std::ostream & out,
	char const * option_name )
{
	if ( option_map.empty() ) return;
	for (
			typename map< T, vector1< S > >::const_iterator
			o = option_map.begin(), oe = option_map.end();
			o != oe; ++o ) {
		out
			<< std::setiosflags(std::ios::left) << setw(16) << option_name
			<< std::setiosflags(std::ios::left) << setw(16) << o->first.id();
		bool first(true);
		for (
				typename vector1< S >::const_iterator
				k = o->second.begin(), ke = o->second.end(); k != ke; ++k ) {
			if ( !first ) {
				out << " ";
			} else {
				first = false;
			}
			write_type( out, *k );
		}
		out << endl;
	}
}

void
JobOptions::show(
	std::ostream & out
) const {
	using namespace utility::options;
	out
		<< std::setiosflags(std::ios::left) << setw(16) << "OptionKeyType"
		<< std::setiosflags(std::ios::left) << setw(16) << "OptionKey"
		<< "OptionValue" << endl;

	show_option_map( boolean_options_, out, "Boolean" );
	show_option_vector_map( boolean_vector_options_, out, "BooleanVector" );
	show_option_map( file_options_, out, "File" );
	show_option_vector_map( file_vector_options_, out, "FileVector" );
	show_option_map( integer_options_, out, "Integer" );
	show_option_vector_map( integer_vector_options_, out, "IntegerVector" );
	show_option_map( path_options_, out, "Path" );
	show_option_vector_map( path_vector_options_, out, "PathVector" );
	show_option_map( real_options_, out, "Real" );
	show_option_vector_map( real_vector_options_, out, "RealVector" );
	show_option_map( string_options_, out, "String" );
	show_option_vector_map( string_vector_options_, out, "StringVector" );
}

std::ostream &
operator << (
	std::ostream & out,
	const JobOptions & job_options
) {
	job_options.show(out);
	return out;
}

template < class T, class S >
void
add_option_to_map( std::map< T, S > & option_map, T const & key, S const & val )
{
	option_map[ key ] = val;
}

template < class T, class S >
void
remove_option_from_map( std::map< T, S > & option_map, T const & key )
{
	typename std::map< T, S >::iterator iter = option_map.find( key );
	if ( iter != option_map.end() ) {
		option_map.erase( iter );
	}
}

template < class T, class S >
bool
option_map_contains_key( std::map< T, S > const & option_map, T const & key )
{
	return option_map.find( key ) != option_map.end();
}

template < class T, class S >
S const &
get_option_from_map( std::map< T, S > const & option_map, T const & key, char const * option_class_name )
{
	typename std::map< T, S >::const_iterator iter = option_map.find( key );
	if ( iter == option_map.end() ) {
		throw EXCN_Msg_Exception( std::string( option_class_name ) + "'" + key.identifier() + "' not found in JobOptions");
	}
	return iter->second;
}

void
JobOptions::add_option(
	BooleanOptionKey key,
	bool val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( BOOLEAN_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	BooleanOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	BooleanOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

bool
JobOptions::get_option(
	BooleanOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "BooleanOptionKey" );
}


void
JobOptions::add_option(
	BooleanVectorOptionKey key,
	vector1< bool > const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( BOOLEAN_VECTOR_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	BooleanVectorOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	BooleanVectorOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

vector1< bool > const &
JobOptions::get_option(
	BooleanVectorOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "BooleanVectorOptionKey" );
}


void
JobOptions::add_option(
	FileOptionKey key,
	FileName const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( FILE_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	FileOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	FileOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

FileName const &
JobOptions::get_option(
	FileOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "FileOptionKey" );
}


void
JobOptions::add_option(
	FileVectorOptionKey key,
	vector1< FileName > const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( FILE_VECTOR_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	FileVectorOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}


bool
JobOptions::has_option(
	FileVectorOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}


vector1< FileName > const &
JobOptions::get_option(
	FileVectorOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "FileVectorOptionKey" );
}


void
JobOptions::add_option(
	IntegerOptionKey key,
	int val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( INTEGER_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	IntegerOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}


bool
JobOptions::has_option(
	IntegerOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

int
JobOptions::get_option(
	IntegerOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "IntegerOptionKey" );
}


void
JobOptions::add_option(
	IntegerVectorOptionKey key,
	vector1< int > const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( INTEGER_VECTOR_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	IntegerVectorOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	IntegerVectorOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

vector1< int > const &
JobOptions::get_option(
	IntegerVectorOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "IntegerVectorOptionKey" );
}


void
JobOptions::add_option(
	PathOptionKey key,
	PathName const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( PATH_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	PathOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	PathOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

PathName const &
JobOptions::get_option(
	PathOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "PathOptionKey" );
}


void
JobOptions::add_option(
	PathVectorOptionKey key,
	vector1< PathName > const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( PATH_VECTOR_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	PathVectorOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	PathVectorOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

vector1< PathName > const &
JobOptions::get_option(
	PathVectorOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "PathVectorOptionKey" );
}


void
JobOptions::add_option(
	RealOptionKey key,
	Real val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( REAL_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	RealOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	RealOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

Real
JobOptions::get_option(
	RealOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "RealOptionKey" );
}


void
JobOptions::add_option(
	RealVectorOptionKey key,
	vector1< Real > const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( REAL_VECTOR_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	RealVectorOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}


bool
JobOptions::has_option(
	RealVectorOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

vector1< Real > const &
JobOptions::get_option(
	RealVectorOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "RealVectorOptionKey" );
}


void
JobOptions::add_option(
	StringOptionKey key,
	string val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( STRING_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	StringOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}


bool
JobOptions::has_option(
	StringOptionKey key
) const {
	return option_map_contains_key( map_for_key( key ), key );
}

string const &
JobOptions::get_option(
	StringOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "StringOptionKey" );
}


void
JobOptions::add_option(
	StringVectorOptionKey key,
	vector1< string > const & val
) {
	if ( track_insertion_order_ ) insertion_order_.push_back( std::make_pair( STRING_VECTOR_OPTION, key.id() ));
	add_option_to_map( map_for_key( key ), key, val );
}

void
JobOptions::remove_option(
	StringVectorOptionKey key
) {
	remove_option_from_map( map_for_key( key ), key );
}

bool
JobOptions::has_option(
	StringVectorOptionKey key
) const {
	return option_map_contains_key( string_vector_options_, key );
}

vector1< string > const &
JobOptions::get_option(
	StringVectorOptionKey key
) const {
	return get_option_from_map( map_for_key( key ), key, "StringVectorOptionKey" );
}

bool
JobOptions::operator == ( JobOptions const & rhs ) const {
	return
		boolean_options_ == rhs.boolean_options_ &&
		boolean_vector_options_ == rhs.boolean_vector_options_ &&
		file_options_ == rhs.file_options_ &&
		file_vector_options_ == rhs.file_vector_options_ &&
		integer_options_ == rhs.integer_options_ &&
		integer_vector_options_ == rhs.integer_vector_options_ &&
		path_options_ == rhs.path_options_ &&
		path_vector_options_ == rhs.path_vector_options_ &&
		real_options_ == rhs.real_options_ &&
		real_vector_options_ == rhs.real_vector_options_ &&
		string_options_ == rhs.string_options_ &&
		string_vector_options_ == rhs.string_vector_options_;
}

void JobOptions::track_insertion_order( bool setting )
{
	track_insertion_order_ = setting;
}

std::list< std::pair< utility::options::OptionTypes, std::string > > const &
JobOptions::insertion_order() const
{
	return insertion_order_;
}


// Type resolution functions to make the above 12 functions much easier to write / maintain
// These functions basically map an input type to an output type for the sake of passing
// type information to the templated functions above -- all of this type resolution
// happens at compile time and should just be optimized out.

std::map< utility::options::BooleanOptionKey, bool > const &
JobOptions::map_for_key( utility::options::BooleanOptionKey const & ) const {
	return boolean_options_;
}

std::map< utility::options::BooleanVectorOptionKey, utility::vector1< bool > > const &
JobOptions::map_for_key( utility::options::BooleanVectorOptionKey const & ) const {
	return boolean_vector_options_;
}

std::map< utility::options::FileOptionKey, utility::file::FileName > const &
JobOptions::map_for_key( utility::options::FileOptionKey const & ) const {
	return file_options_;
}

std::map< utility::options::FileVectorOptionKey, utility::vector1< utility::file::FileName > > const &
JobOptions::map_for_key( utility::options::FileVectorOptionKey const & ) const {
	return file_vector_options_;
}

std::map< utility::options::IntegerOptionKey, int > const &
JobOptions::map_for_key( utility::options::IntegerOptionKey const & ) const {
	return integer_options_;
}

std::map< utility::options::IntegerVectorOptionKey, utility::vector1< int > > const &
JobOptions::map_for_key( utility::options::IntegerVectorOptionKey const & ) const {
	return integer_vector_options_;
}

std::map< utility::options::PathOptionKey, utility::file::PathName > const &
JobOptions::map_for_key( utility::options::PathOptionKey const & ) const {
	return path_options_;
}

std::map< utility::options::PathVectorOptionKey, utility::vector1< utility::file::PathName > > const &
JobOptions::map_for_key( utility::options::PathVectorOptionKey const & ) const {
	return path_vector_options_;
}

std::map< utility::options::RealOptionKey, platform::Real > const &
JobOptions::map_for_key( utility::options::RealOptionKey const & ) const {
	return real_options_;
}

std::map< utility::options::RealVectorOptionKey, utility::vector1< platform::Real > > const &
JobOptions::map_for_key( utility::options::RealVectorOptionKey const & ) const {
	return real_vector_options_;
}

std::map< utility::options::StringOptionKey, std::string > const &
JobOptions::map_for_key( utility::options::StringOptionKey const & ) const {
	return string_options_;
}

std::map< utility::options::StringVectorOptionKey, utility::vector1< std::string > > const &
JobOptions::map_for_key( utility::options::StringVectorOptionKey const & ) const {
	return string_vector_options_;
}

// Non-const versions of all the above functions

std::map< utility::options::BooleanOptionKey, bool > &
JobOptions::map_for_key( utility::options::BooleanOptionKey const & ) {
	return boolean_options_;
}

std::map< utility::options::BooleanVectorOptionKey, utility::vector1< bool > > &
JobOptions::map_for_key( utility::options::BooleanVectorOptionKey const & ) {
	return boolean_vector_options_;
}

std::map< utility::options::FileOptionKey, utility::file::FileName > &
JobOptions::map_for_key( utility::options::FileOptionKey const & ) {
	return file_options_;
}

std::map< utility::options::FileVectorOptionKey, utility::vector1< utility::file::FileName > > &
JobOptions::map_for_key( utility::options::FileVectorOptionKey const & ) {
	return file_vector_options_;
}

std::map< utility::options::IntegerOptionKey, int > &
JobOptions::map_for_key( utility::options::IntegerOptionKey const & ) {
	return integer_options_;
}

std::map< utility::options::IntegerVectorOptionKey, utility::vector1< int > > &
JobOptions::map_for_key( utility::options::IntegerVectorOptionKey const & ) {
	return integer_vector_options_;
}

std::map< utility::options::PathOptionKey, utility::file::PathName > &
JobOptions::map_for_key( utility::options::PathOptionKey const & ) {
	return path_options_;
}

std::map< utility::options::PathVectorOptionKey, utility::vector1< utility::file::PathName > > &
JobOptions::map_for_key( utility::options::PathVectorOptionKey const & ) {
	return path_vector_options_;
}

std::map< utility::options::RealOptionKey, platform::Real > &
JobOptions::map_for_key( utility::options::RealOptionKey const & ) {
	return real_options_;
}

std::map< utility::options::RealVectorOptionKey, utility::vector1< platform::Real > > &
JobOptions::map_for_key( utility::options::RealVectorOptionKey const & ) {
	return real_vector_options_;
}

std::map< utility::options::StringOptionKey, std::string > &
JobOptions::map_for_key( utility::options::StringOptionKey const & ) {
	return string_options_;
}

std::map< utility::options::StringVectorOptionKey, utility::vector1< std::string > > &
JobOptions::map_for_key( utility::options::StringVectorOptionKey const & ) {
	return string_vector_options_;
}



} // namespace
} // namespace
