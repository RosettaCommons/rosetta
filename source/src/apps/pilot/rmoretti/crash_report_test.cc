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

#include <core/types.hh>
#include <basic/options/option.hh>

#include <utility/exit.hh>
#include <utility/crash_report.hh>
#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <string>
#include <list>

int return_zero(int x) {
	if ( x >= 0 ) { return 0; }
	else { return 1; }
}

// HIGHLY dangerous function which *will* exhaust your memory.
// (Use something like `ulimit -m 2048000; ulimit -H -v 2048000` in a new terminal shell to keep things from bringing down your machine entirely.)
template< class Container >
void do_oom() {
	Container values;
	core::Size counter = 0;
	while ( true ) {
		values.push_back( ++counter );
		if ( counter % 10000000 == 0 ) {
			std::cout << "Causing OOM: " << values.size() << std::endl;
		}
	}
}


void do_two(std::string const & type) {
	if ( type == "manual" ) {
		do_crash("Die Hard");
	} else if ( type == "utility_exit_with_message" ) {
		utility_exit_with_message("Crash Exit");
	} else if ( type == "utility_exit" ) {
		utility_exit();
	} else if ( type == "user_exit" ) {
		user_fixable_issue_exit("Try fixing issue number ", 1, " and re-running.");
	} else if ( type == "user_assert" ) {
		user_fixable_issue_assert( 2 + 2 == 5, "Try running with slighly larger values of ", 2 );
	} else if ( type == "runtime_assert" ) {
		runtime_assert(false == true);
	} else if ( type == "debug_assert" ) {
		debug_assert(1 == 2);
	} else if ( type == "throw" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Bad Dates, Indy!");
	} else if ( type == "throw_user" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "User made a boo-boo.");
	} else if ( type == "throw_std" ) {
		throw std::logic_error("Manual throw of non-Rosetta Error");
	} else if ( type == "throw_map" ) {
		std::map<int,int> const mapp = { {1,2}, {3,4} };
		std::cout << mapp.at(2) << std::endl; // Should throw std::out_of_range
	} else if ( type == "throw_string" ) {
		throw std::string("This is a literal string I am throwing");
	} else if ( type == "terminate" ) {
		std::terminate(); // Direct call
	} else if ( type == "zero" ) {
		int f = 5 / return_zero(1);
		int g = 5 / ( return_zero(-1) - 1);
		std::cout << "Computed: " << f << " " << g << std::endl;
	} else if ( type == "segfault" ) {
		int * p = nullptr;
		std::cout << "I am " << return_zero( *p ) << std::endl;
	} else if ( type == "OOM" ) {
		do_oom< std::list<core::Size> >(); // List will grab memory in small bites
	} else if ( type == "OOM_vect" ) {
		do_oom< utility::vector1<core::Size> >(); // The exponential resize of vector should leave enough memory for reporter.
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Crash type "+type+" is unknown.");
	}
	// Other things to check for (edge cases in pre-initialization loading.)
	// * -HCF
	// * extra 'free' option on the commandline
	// * miss-spelled options
	// * missing/malformed option file
}

void do_one(std::string const & type) {
	do_two(type);
}

int do_main(std::string const & type) {
	do_one(type);
	return 0;
}

int
main( int argc, char * argv [] )
{

	try {

		devel::init(argc, argv);

		if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].size() == 0 ) {
			std::cout << "NEED TO PASS -s" << std::endl;
			return -1;
		}

		return do_main(basic::options::option[ basic::options::OptionKeys::in::file::s ][1] );

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}

