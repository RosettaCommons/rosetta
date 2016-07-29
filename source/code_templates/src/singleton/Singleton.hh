// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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

	// Singleton-creation function for use with utility::thread::threadsafe_singleton
	static --class-- * create_singleton_instance();


private:  // Private data /////////////////////////////////////////////////////

};

--end_namespace--


#endif //INCLUDED_--path--_--class--_fwd_hh



