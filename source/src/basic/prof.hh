// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_basic_prof_hh
#define INCLUDED_basic_prof_hh

#include <basic/prof.fwd.hh> // For some enum definitions

#include <basic/Tracer.fwd.hh>
#include <time.h>                                    // for clock
#include <basic/options/keys/run.OptionKeys.gen.hh>  // for profile
#include <basic/options/option.hh>                   // for OptionCollection
#include <string>                                    // for allocator
#include <utility/options/BooleanOption.hh>          // for BooleanOption
#include <utility/vector1.hh>                        // for vector1
#include <utility/vectorL.hh>                        // for vectorL


namespace basic {
/**
not intended for profiling inside tight loops
the clock() routine currently being used has fairly crappy
resolution and it will introduce some small overhead with the
function calls and the if-check even if not using -profile on
the command line

you can wrap it around large-ish chunks of code, like fullatom_energy
or rotamer_trials...

A simple setup for timing code fragments. Probably not optimal timing
functions -- I'm open to suggestions.

looks like (see eg fullatom_energy or scorefxn)

PROF_START( prof::TAG );

<function call>

PROF_STOP( prof::TAG );

where TAG is in the enum "Prof_tag" below (feel free to add new ones)
also add to tag2string if you want friendly output.

PROF_STOP checks the time and increments the total time assigned to TAG


2. later on, in your simulation code you can do:

prof_reset();

-- miscellaneous simulation --

prof_show();

The final call to prof::show() will display the time usage measured
by all the PROF_* calls between reset() and show()
**/

#ifdef NO_PROF
#define PROF_START(expr)
#define PROF_STOP(expr)
#else
#define PROF_START(expr) ( prof_start_function_body( expr) )
#define PROF_STOP(expr) ( prof_stop_function_body( expr ) )
#endif

extern utility::vector1< std::string > tag2string;
extern utility::vector1< clock_t > start_clock;
extern utility::vector1< double > total_clock;
extern utility::vector1< int > calls;
extern utility::vector1< int > bad_calls;

extern double const clock_factor;
extern clock_t const SHRINK_FACTOR; // prevent overflow

inline
void
prof_start_function_body( ProfTag const tag )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	start_clock[ tag ] = clock() / SHRINK_FACTOR;
	++calls[ tag ];
}


inline
void
prof_stop_function_body( ProfTag const tag )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	clock_t const current( clock() / SHRINK_FACTOR );
	clock_t const start( start_clock[tag] );

	if ( current >= start ) {
		total_clock[ tag ] += clock_factor * ( current - start );
	} else {
		--calls[ tag ];
		++bad_calls[ tag ];
	}
}

void prof_reset();

void prof_show();

class ProfileThis {
public:
	ProfileThis( ProfTag tag ) : tag_( tag ) {
		PROF_START( tag_ );
	}
	~ProfileThis() {
		PROF_STOP( tag_ );
	}
private:
	ProfTag tag_;
};

class DynamicProfileThis {
public:
	DynamicProfileThis( std::string const& prof_tag );
	~DynamicProfileThis();
private:
	clock_t start_clock_;
	std::string tag_;
};

/// @brief print "TIME_STAMP: Www Mmm dd hh:mm:ss yyyy msg" on tr.Error and on std::cerr (if boolean is true)
extern bool show_time_on_cerr;
void show_time( basic::Tracer& tr, std::string const& msg );

} // basic

#endif
