// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   project/subproject/BaseClass.hh
///
/// @brief
/// @author


#ifndef INCLUDED_project_subproject_BaseClass_HH
#define INCLUDED_project_subproject_BaseClass_HH


// Project forward headers
#include <project/subproject/BaseClass.fwd.hh>


// Project headers


// External library headers


// C++ headers


// Operating system headers


// Forward declarations




namespace project {
namespace subproject {


	/// @brief
class BaseClass
{

  // Friends


public: // Types


protected: // Types


private: // Types




public: // Constants


protected: // Constants


private: // Constants




public: // Creation


	/// @brief Destructor.
	virtual
	~BaseClass();


protected: // Creation


	/// @brief Constructor.
	BaseClass();


private: // Creation


	/// @brief Noncopyable.  Undefined
	BaseClass( BaseClass const & );




private: // Methods: assignment


	/// @brief Noncopyable.  Undefined.
	BaseClass&
	operator=( BaseClass const & );


public: // Methods: comparison




public: // Methods
	// Further subsections of methods allowed


protected: // Methods
	// Further subsections of methods allowed


private: // Methods
	// Further subsections of methods allowed




public: // Properties


protected: // Properties




private: // Fields


}; // BaseClass


} // namespace subproject
} // namespace project


#endif // INCLUDED_project_subproject_BaseClass_HH
