// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   project/subproject/DerivedClass.hh
///
/// @brief
/// @author


#ifndef INCLUDED_project_subproject_DerivedClass_HH
#define INCLUDED_project_subproject_DerivedClass_HH


// Project forward headers
#include <project/subproject/DerivedClass.fwd.hh>


// Project headers
#include <project/subproject/BaseClass.hh>


// External library headers


// C++ headers


// Operating system headers


// Forward declarations




namespace project {
namespace subproject {


	/// @brief
class DerivedClass : public BaseClass
{

  // Friends


public: // Types


	typedef BaseClass Super;


protected: // Types


private: // Types




public: // Constants


protected: // Constants


private: // Constants




public: // Creation


	/// @brief Constructor
	DerivedClass();


	/// @brief Destructor
	virtual
	~DerivedClass();


protected: // Creation


private: // Creation


	/// @brief Noncopyable.  Undefined.
	DerivedClass( DerivedClass const & );




private: // Methods: assignment


	/// @brief Noncopyable.  Undefined.
	DerivedClass&
	operator=( DerivedClass const & );


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


}; // DerivedClass


} // namespace subproject
} // namespace project


#endif // INCLUDED_project_subproject_DerivedClass_HH
