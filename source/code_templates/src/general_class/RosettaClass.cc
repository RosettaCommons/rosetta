// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#include <--path--/--class--.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );


--namespace--

--class--::--class--():
 utility::pointer::ReferenceCount()
{

}

--class--::~--class--(){}

--class--::--class--( --class-- const & ) {

}



--class--OP
--class--::clone() const {
	return --class--OP( new --class--( *this ) );
}


--end_namespace--






