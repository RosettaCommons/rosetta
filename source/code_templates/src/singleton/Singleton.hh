// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--
/// @brief --brief--
/// @author --name-- (--email--)


#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

// Unit headers

// Utility header
#include <utility/SingletonBase.hh>

// C++ header

--namespace--

/// @brief --brief--
class --class-- : public utility::SingletonBase< --class-- > {
	friend class utility::SingletonBase< --class-- >;

private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	--class--();
	--class--(--class-- const & ) = delete;
	--class-- operator=(--class-- const & ) = delete;

private:  // Private data /////////////////////////////////////////////////////

};

--end_namespace--

#endif //INCLUDED_--path--_--class--_fwd_hh



