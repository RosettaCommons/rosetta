// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/after_opts.cc
/// @brief  Option lookup functions emulating Rosetta++ equivalents for transitional use
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com) (port to basic::options)
/// @author Original author(s) unknown


// Unit headers
#include <basic/options/after_opts.hh>

// Package headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

// Utility headers
#include <utility/exit.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>

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
#include <utility/keys/KeyLookup.hh>
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
#include <utility/options/keys/BooleanOptionKey.hh>
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
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>



namespace basic {
namespace options {


/// @brief Get key for an option name
OptionKey const &
key( std::string const & str )
{
	return OptionKeys::key( OptionCollection::find_key_cl( str, std::string(), true ) );
}


bool
truefalseoption( std::string const & str )
{
	return truefalseoption( key( str ) );
}


bool
truefalseoption( OptionKey const & key )
{
	if ( option[ key ].user() ) {
		std::cout << "[T/F  OPT]New TRUE value for [-" << key.id() << ']' << std::endl;
		return true;
	} else {
		std::cout << "[T/F  OPT]Default FALSE value for [-" << key.id() << ']' << std::endl;
		return false;
	}
}


double
realafteroption(
	std::string const & str
)
{
	if ( ! truefalseoption( str ) ) {
		std::cerr << "realafteroption: missing a required argument: " << str << std::endl;
		utility_exit();
	}
	double val;
	realafteroption( str, 100, val );
	return val;
}


double
realafteroption(
	std::string const & str,
	double const opt_default
)
{
	double value;
	realafteroption( key( str ), opt_default, value );
	return value;
}


void
realafteroption(
	std::string const & str,
	double const opt_default,
	double & rnum
)
{
	realafteroption( key( str ), opt_default, rnum );
}


void
realafteroption(
	OptionKey const & key,
	double const opt_default,
	double & rnum
)
{
	using namespace ObjexxFCL::format;
	if ( option[ key ].user() ) {
		rnum = option.option< RealOption >( key );
		std::cout << "[REAL OPT]New value for [-" << key.id() << "] " << SS( rnum ) << std::endl;
	} else {
		rnum = opt_default;
		std::cout << "[REAL OPT]Default value for [-" << key.id() << "] " << SS( rnum ) << std::endl;
	}
}


void
real2afteroption(
	std::string const & str,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2
)
{
	real2afteroption( key( str ), default1, rnum1, default2, rnum2 );
}


void
real2afteroption(
	OptionKey const & key,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2
)
{
	using namespace ObjexxFCL::format;
	if ( option[ key ].user() ) {
		RealVectorOption const & opt( option.option< RealVectorOption >( key ) );
		if ( opt.size() < 2 ) {
			std::cerr << "real2afteroption: option specified with < 2 values: -" << key.id() << std::endl;
			utility_exit();
		}
		rnum1 = opt[ 1 ];
		rnum2 = opt[ 2 ];
		std::cout << "[REAL OPT]New values for [-" << key.id() << "] " << SS( rnum1 ) << SS( rnum2 ) << std::endl;
	} else {
		rnum1 = default1;
		rnum2 = default2;
		std::cout << "[REAL OPT]Default value for [-" << key.id() << "] " << SS( rnum1 ) << SS( rnum2 ) << std::endl;
	}
}


void
real3afteroption(
	std::string const & str,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2,
	double const default3,
	double & rnum3
)
{
	real3afteroption( key( str ), default1, rnum1, default2, rnum2, default3, rnum3 );
}


void
real3afteroption(
	OptionKey const & key,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2,
	double const default3,
	double & rnum3
)
{
	using namespace ObjexxFCL::format;
	if ( option[ key ].user() ) {
		RealVectorOption const & opt( option.option< RealVectorOption >( key ) );
		if ( opt.size() < 3 ) {
			std::cerr << "real3afteroption: option specified with < 3 values: -" << key.id() << std::endl;
			utility_exit();
		}
		rnum1 = opt[ 1 ];
		rnum2 = opt[ 2 ];
		rnum3 = opt[ 3 ];
		std::cout << "[REAL OPT]New values for [-" << key.id() << "] " << SS( rnum1 ) << SS( rnum2 ) << SS( rnum3 ) << std::endl;
	} else {
		rnum1 = default1;
		rnum2 = default2;
		rnum3 = default3;
		std::cout << "[REAL OPT]Default value for [-" << key.id() << "] " << SS( rnum1 ) << SS( rnum2 ) << SS( rnum3 ) << std::endl;
	}
}


int
intafteroption(
	std::string const & str
)
{
	if ( !truefalseoption( str ) ) {
		std::cout << "intafteroption: missing required arg: " << str << std::endl;
		utility_exit();
	}
	int val;
	intafteroption( key( str ), 100, val );
	return val;
}


int
intafteroption(
	std::string const & str,
	int const opt_default
)
{
	int val;
	intafteroption( key( str ), opt_default, val );
	return val;
}


void
intafteroption(
	std::string const & str,
	int const opt_default,
	int & inum
)
{
	intafteroption( key( str ), opt_default, inum );
}


void
intafteroption(
	OptionKey const & key,
	int const opt_default,
	int & inum
)
{
	using namespace ObjexxFCL::format;
	if ( option[ key ].user() ) {
		inum = option.option< IntegerOption >( key );
		std::cout << "[INT OPT]New value for [-" << key.id() << "] " << SS( inum ) << std::endl;
	} else {
		inum = opt_default;
		std::cout << "[INT OPT]Default value for [-" << key.id() << "] " << SS( inum ) << std::endl;
	}
}


void
optional_positive_intafteroption(
	std::string const & str,
	int const opt_default,
	int & inum
)
{
	optional_positive_intafteroption( key( str ), opt_default, inum );
}


void
optional_positive_intafteroption(
	OptionKey const & key,
	int const opt_default,
	int & inum
)
{
	using namespace ObjexxFCL::format;
	if ( option[ key ].user() ) {
		IntegerVectorOption const & opt( option.option< IntegerVectorOption >( key ) );
		if ( opt.size() >= 1 ) {
			inum = opt[ 1 ];
			std::cout << "[INT OPT]New value for [-" << key.id() << "] " << SS( inum ) << std::endl;
		} else {
			inum = opt_default;
			std::cout << "[INT OPT]Default value for [-" << key.id() << "] " << SS( inum ) << std::endl;
		}
	} else {
		inum = opt_default;
		std::cout << "[INT OPT]Default value for [-" << key.id() << "] " << SS( inum ) << std::endl;
	}
}


void
int2afteroption(
	std::string const & str,
	int const opt_default,
	int & inum,
	int const opt_default2,
	int & inum2
)
{
	int2afteroption( key( str ), opt_default, inum, opt_default2, inum2 );
}


void
int2afteroption(
	OptionKey const & key,
	int const opt_default,
	int & inum,
	int const opt_default2,
	int & inum2
)
{
	using namespace ObjexxFCL::format;
	if ( option[ key ].user() ) {
		IntegerVectorOption const & opt( option.option< IntegerVectorOption >( key ) );
		if ( opt.size() < 2 ) {
			std::cerr << "int2afteroption: option specified with < 2 values: -" << key.id() << std::endl;
			utility_exit();
		}
		inum  = opt[ 1 ];
		inum2 = opt[ 2 ];
		std::cout << "[INT OPT]New values for [-" << key.id() << "] " << SS( inum ) << SS( inum2 ) << std::endl;
	} else {
		inum  = opt_default;
		inum2 = opt_default2;
		std::cout << "[INT OPT]Default values for [-" << key.id() << "] " << SS( inum ) << SS( inum2 ) << std::endl;
	}
}


// require presence of the option

std::string
stringafteroption(
	std::string const & str
)
{
	if ( ! truefalseoption( str ) ) {
		std::cout << "STOP:: stringafteroption: missing a required argument: " << str << std::endl;
		utility_exit();
	}
	std::string val;
	stringafteroption( key( str ), "dummy", val );
	return val;
}


std::string
stringafteroption(
	std::string const & str,
	std::string const & opt_default
)
{
	std::string val;
	stringafteroption( key( str ), opt_default, val );
	return val;
}


void
stringafteroption(
	std::string const & str,
	std::string const & opt_default,
	std::string & cval
)
{
	stringafteroption( key( str ), opt_default, cval );
}


void
stringafteroption(
	OptionKey const & key,
	std::string const & opt_default,
	std::string & cval
)
{
	if ( option[ key ].user() ) {
		cval = option.option< StringOption >( key );
		std::cout << "[STR  OPT]New value for [-" << key.id() << "] " << cval << '.' << std::endl;
	} else {
		cval = opt_default;
		std::cout << "[STR  OPT]Default value for [-" << key.id() << "] " << cval << '.' << std::endl;
	}
}


void
stringafteroption(
	std::string const & str,
	char const opt_default,
	char & cval
)
{
	stringafteroption( key( str ), opt_default, cval );
}


void
stringafteroption(
	OptionKey const & key,
	char const opt_default,
	char & cval
)
{
	if ( option[ key ].user() ) {
		StringOption const & opt( option.option< StringOption >( key ) );
		cval = ( opt().length() > 0 ? opt()[ 0 ] : ' ' );
		std::cout << "[STR  OPT]New value for [-" << key.id() << "] " << cval << '.' << std::endl;
	} else {
		cval = opt_default;
		std::cout << "[STR  OPT]Default value for [-" << key.id() << "] " << cval << '.' << std::endl;
	}
}


} // namespace options
} // namespace basic
