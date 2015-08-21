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

#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/keys/KeyLookup.hh>
#include <basic/options/option.hh>

#include <utility/exit.hh>

//Auto Headers
#include <platform/types.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
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
