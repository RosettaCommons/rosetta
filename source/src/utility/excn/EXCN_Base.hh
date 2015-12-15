// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/excn/EXCN_Base.hh
/// @brief  base class for Exception system
/// @author Oliver Lange


#ifndef INCLUDED_utility_excn_EXCN_Base_HH
#define INCLUDED_utility_excn_EXCN_Base_HH

// Unit Headers
#include <utility/excn/EXCN_Base.fwd.hh>

// Package Headers
#include <utility/assert.hh>
#include <cstdlib>
#include <ostream>
#include <sstream>

namespace utility {
namespace excn {

/* *********************************************************************************************************
************************************************************************************************************
*********************                                                                  *********************
*********************                        W  A  R  N  I  N  G                       *********************
*********************                                                                  *********************
************************************************************************************************************
************************************************************************************************************
don't use this Exception definition anywhere apart from
in the jd2::JobDistributor
( or in your main if you don't use the JobDistributor )
instead derive your Exceptions from EXCN_Exception
and only catch EXCN_Exception or child-classes
( except in main and jdist, where you should catch EXCN_Base )

Oliver <olange@u.washington.edu>
Matthew O'Meara <mattjomeara@gmail.com>
************************************************************************************************************
************************************************************************************************************
*/
class EXCN_Base {
protected:
	EXCN_Base() {
		//would like to add an option run:no_exceptions
		// * if ( option[ run::no_exceptions ] ) {
		//exit - the hard way:
		//damn can't call virtual from constructor....
		//show( std::cerr );


		//* THE ASSERT MACRO is here that one can find the origin of the EXCEPTION in gdb */
		// debug_assert( false );
		// a better method for this is to issue the command
		// catch throw
		// in the gdb command line... now gdb will stop execution when an Exception
		// http://www.delorie.com/gnu/docs/gdb/gdb_31.html
		/* IN RELEASE MODE THIS HAS CONSTRUCTOR MUST NOT FAIL! --- otherwise the ERROR Msg get's lost! */
	};
public:
	virtual ~EXCN_Base() {};
	virtual void show( std::ostream& ) const = 0;
	virtual std::string const msg() const {
		std::string msg;
		std::ostringstream os;
		show( os );
		msg = os.str();
		return msg;
	};
};

inline std::ostream& operator << ( std::ostream& os, EXCN_Base const & excn ) {
	excn.show( os );
	return os;
}

}
}
#endif
