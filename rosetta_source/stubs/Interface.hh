// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   project/subproject/Interface.hh
///
/// @brief
/// @author
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.


#ifndef INCLUDED_project_subproject_Interface_HH
#define INCLUDED_project_subproject_Interface_HH


// Project forward headers
#include <project/subproject/Interface.fwd.hh>


// Other project headers


// External library headers


// C++ headers


// Operating system headers


// Forward declarations




namespace project {
namespace subproject {


	/// @brief
class Interface
{
  // Friends


public: // Types


protected: // Types




public: // Constants


protected: // Constants




public: // Creation


	/// @brief Destructor
	virtual
	~Interface() {};


protected: // Creation


	/// @brief Prevent direct instantiation: No other constructors allowed.
	Interface() {};


public: // Methods
	// Further subsections of methods allowed


protected: // Methods
	// Further subsections of methods allowed




public: // Properties


protected: // Properties


	// NO FIELDS ALLOWED


}; // Interface


} // namespace subproject
} // namespace project


#endif // INCLUDED_project_subproject_Interface_HH
