// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   project/subproject/FinalClass.hh
///
/// @brief
/// @author


#ifndef INCLUDED_project_subproject_FinalClass_HH
#define INCLUDED_project_subproject_FinalClass_HH


// Project forward headers
#include <project/subproject/FinalClass.fwd.hh>


// Project headers
#include <project/subproject/BaseClass.hh>


// External library headers


// C++ headers


// Operating system headers


// Forward declarations




namespace project {
namespace subproject {


	/// @brief
class FinalClass : public BaseClass
{

  // Friends


public: // Types


	typedef BaseClass Super;


private: // Types




public: // Constants


private: // Constants




public: // Creation


	/// @brief Constructor
	FinalClass();


	/// @brief Destructor
	virtual ~FinalClass();


	/// @brief Copy constructor
	FinalClass( FinalClass const & );


private: // Creation




public: // Methods: assignment


	/// @brief operator=
	FinalClass&
	operator=( FinalClass const & );


public: // Methods: comparison




public: // Methods
	// Further subsections of methods allowed


private: // Methods
	// Further subsections of methods allowed



public: // Properties




private: // Fields


}; // FinalClass


} // namespace subproject
} // namespace project


#endif // INCLUDED_project_subproject_FinalClass_HH
