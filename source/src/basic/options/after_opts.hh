// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/after_opts.hh
/// @brief  Option lookup functions emulating Rosetta++ equivalents for transitional use
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com) (port to basic::options)
/// @author Original author(s) unknown


#ifndef INCLUDED_basic_options_after_opts_hh
#define INCLUDED_basic_options_after_opts_hh


// Utility headers
#include <utility/options/keys/OptionKey.fwd.hh>

// C++ headers
#ifdef WIN32
#include <string>
#else
#include <iosfwd>
#endif


namespace basic {
namespace options {


typedef  utility::options::OptionKey  OptionKey;


bool
truefalseoption( std::string const & str );


bool
truefalseoption( OptionKey const & key );


double
realafteroption(
	std::string const & str
);


double
realafteroption(
	std::string const & str,
	double const opt_default
);


void
realafteroption(
	std::string const & str,
	double const opt_default,
	double & rnum
);


void
realafteroption(
	OptionKey const & key,
	double const opt_default,
	double & rnum
);


void
real2afteroption(
	std::string const & str,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2
);


void
real2afteroption(
	OptionKey const & key,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2
);


void
real3afteroption(
	std::string const & str,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2,
	double const default3,
	double & rnum3
);


void
real3afteroption(
	OptionKey const & key,
	double const default1,
	double & rnum1,
	double const default2,
	double & rnum2,
	double const default3,
	double & rnum3
);


int
intafteroption(
	std::string const & str
);


int
intafteroption(
	std::string const & str,
	int const opt_default
);


void
intafteroption(
	std::string const & str,
	int const opt_default,
	int & inum
);


void
intafteroption(
	OptionKey const & key,
	int const opt_default,
	int & inum
);


void
int2afteroption(
	std::string const & str,
	int const opt_default,
	int & inum,
	int const opt_default2,
	int & inum2
);


void
int2afteroption(
	OptionKey const & key,
	int const opt_default,
	int & inum,
	int const opt_default2,
	int & inum2
);


void
optional_positive_intafteroption(
	std::string const & str,
	int const opt_default,
	int & inum
);


void
optional_positive_intafteroption(
	OptionKey const & key,
	int const opt_default,
	int & inum
);


std::string
stringafteroption(
	std::string const & str
);


std::string
stringafteroption(
	std::string const & str,
	std::string const & opt_default
);


void
stringafteroption(
	std::string const & str,
	std::string const & opt_default,
	std::string & cval
);


void
stringafteroption(
	OptionKey const & key,
	std::string const & opt_default,
	std::string & cval
);


void
stringafteroption(
	std::string const & str,
	char const opt_default,
	char & cval
);


void
stringafteroption(
	OptionKey const & key,
	char const opt_default,
	char & cval
);


} // namespace options
} // namespace basic


#endif // INCLUDED_basic_options_after_opts_HH
