// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/options/PyRosettaFunctions.hh
/// @brief  Additional functions to set/get Options from PyRosetta
/// @author Sergey Lyskov


#include <basic/options/PyRosettaFunctions.hh>
#include <basic/options/option.hh>                         // for option
#include <iosfwd>                                          // for string
#include <string>                                          // for allocator
#include <utility/exit.hh>                                 // for utility_ex...
#include <utility/file/FileName.hh>                        // for operator<
#include <utility/keys/AutoKey.hh>                         // for operator<
#include <utility/keys/KeyLookup.hh>                       // for has, key
#include <utility/options/keys/BooleanOptionKey.hh>        // for BooleanOpt...
#include <utility/options/keys/BooleanVectorOptionKey.hh>  // for BooleanVec...
#include <utility/options/keys/FileOptionKey.hh>           // for FileOptionKey
#include <utility/options/keys/FileVectorOptionKey.hh>     // for FileVector...
#include <utility/options/keys/IntegerOptionKey.hh>        // for IntegerOpt...
#include <utility/options/keys/IntegerVectorOptionKey.hh>  // for IntegerVec...
#include <utility/options/keys/OptionKeys.hh>              // for KeyType
#include <utility/options/keys/RealOptionKey.hh>           // for RealOptionKey
#include <utility/options/keys/RealVectorOptionKey.hh>     // for RealVector...
#include <utility/options/keys/StringOptionKey.hh>         // for StringOpti...
#include <utility/options/keys/StringVectorOptionKey.hh>   // for StringVect...
#include <utility/vector1.hh>                              // for vector1
#include <utility/vector1_bool.hh>                         // for vector1
#include <basic/options/keys/OptionKeys.hh>


namespace basic {
namespace options {


template< typename T, typename K >
T get_option(std::string const & id)
{
	if ( !utility::options::OptionKeys::has( id ) ) utility_exit_with_message( "get_option: OptionKey with id " + id + " not found!" );

	return basic::options::option[ dynamic_cast<K const &>( utility::options::OptionKeys::key( id ) ) ].value();
}

template< typename T, typename K >
void set_option(std::string const & id, T const & value)
{
	if ( !utility::options::OptionKeys::has( id ) )  utility_exit_with_message( "set_option: OptionKey with id " + id + " not found!" );

	basic::options::option[ dynamic_cast<K const &>( utility::options::OptionKeys::key( id ) ) ].value( value );
}


bool        get_boolean_option(std::string const & id) { return get_option<        bool, utility::options::BooleanOptionKey> (id); }
int         get_integer_option(std::string const & id) { return get_option<         int, utility::options::IntegerOptionKey> (id); }
double      get_real_option(std::string const & id)    { return get_option<      double, utility::options::RealOptionKey>    (id); }
std::string get_string_option(std::string const & id)  { return get_option< std::string, utility::options::StringOptionKey>  (id); }
std::string get_file_option(std::string const & id)    { return get_option< std::string, utility::options::FileOptionKey>    (id); }

void set_boolean_option(std::string const & id, bool v)                { return set_option<        bool, utility::options::BooleanOptionKey> (id, v); }
void set_integer_option(std::string const & id, int v)                 { return set_option<         int, utility::options::IntegerOptionKey> (id, v); }
void set_real_option(std::string const & id, double v)                 { return set_option<      double, utility::options::RealOptionKey>    (id, v); }
void set_string_option(std::string const & id, std::string const & v)  { return set_option< std::string, utility::options::StringOptionKey>  (id, v); }
void set_file_option(std::string const & id, std::string const & v)      { return set_option< std::string, utility::options::FileOptionKey>    (id, v); }


utility::vector1<bool>        get_boolean_vector_option(std::string const & id) { return get_option< utility::vector1<bool>,        utility::options::BooleanVectorOptionKey> (id); }
utility::vector1<int>         get_integer_vector_option(std::string const & id) { return get_option< utility::vector1<int>,         utility::options::IntegerVectorOptionKey> (id); }
utility::vector1<double>      get_real_vector_option(std::string const & id)    { return get_option< utility::vector1<double>,      utility::options::RealVectorOptionKey>    (id); }
utility::vector1<std::string> get_string_vector_option(std::string const & id)  { return get_option< utility::vector1<std::string>, utility::options::StringVectorOptionKey>  (id); }
utility::vector1<std::string> get_file_vector_option(std::string const & id)    { return get_option< utility::vector1<std::string>, utility::options::FileVectorOptionKey>    (id); }


void set_boolean_vector_option(std::string const & id, utility::vector1<bool> const & v)       { return set_option< utility::vector1<bool>,        utility::options::BooleanVectorOptionKey> (id, v); }
void set_integer_vector_option(std::string const & id, utility::vector1<int> const & v)        { return set_option< utility::vector1<int>,         utility::options::IntegerVectorOptionKey> (id, v); }
void set_real_vector_option(std::string const & id, utility::vector1<double> const & v)        { return set_option< utility::vector1<double>,      utility::options::RealVectorOptionKey>    (id, v); }
void set_string_vector_option(std::string const & id, utility::vector1<std::string> const & v) { return set_option< utility::vector1<std::string>, utility::options::StringVectorOptionKey>  (id, v); }
void set_file_vector_option(std::string const & id, utility::vector1<std::string> const & v)   { return set_option< utility::vector1<std::string>, utility::options::FileVectorOptionKey>    (id, v); }


} // namespace options
} // namespace basic
