// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/after_opts.cc
/// @brief  Option lookup functions emulating Rosetta++ equivalents for transitional use
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com) (port to basic::options)
/// @author Original author(s) unknown


// Unit headers
#include <basic/options/after_opts.hh>
#include <ObjexxFCL/format.hh>                 // for SS
#include <basic/options/keys/OptionKeys.hh>    // for OptionKey, KeyType
#include <basic/options/option.hh>             // for OptionCollection, Real...
#include <iostream>                            // for operator<<, basic_ostream
#include <string>                              // for char_traits, allocator
#include <utility/exit.hh>                     // for utility_exit
#include <utility/keys/AutoKey.hh>             // for operator<, operator==
#include <utility/keys/KeyLookup.hh>           // for key
#include <utility/options/IntegerOption.hh>    // for IntegerOption
#include <utility/options/Option.hh>           // for Option
#include <utility/options/RealOption.hh>       // for RealOption
#include <utility/options/ScalarOption_T_.hh>  // for ScalarOption_T_, Scala...
#include <utility/options/StringOption.hh>     // for StringOption
#include <utility/options/VectorOption_T_.hh>  // for VectorOption_T_

#include <basic/options/keys/OptionKeys.hh>


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
		std::cerr << "intafteroption: missing required arg: " << str << std::endl;
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
