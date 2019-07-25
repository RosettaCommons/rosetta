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

// Project headers:
#include <--path--/--class--.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "--namespace_dot--.--class--" );


--namespace--

/// @brief Default constructor.
--class--::--class--():
 utility::pointer::ReferenceCount()
{}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
--class--::--class--( --class-- const & )=default;

/// @brief Destructor.
--class--::~--class--(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
--class--OP
--class--::clone() const {
	return utility::pointer::make_shared< --class-- >( *this );
}

--end_namespace--