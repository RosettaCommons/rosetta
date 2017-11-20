// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Exceptions.hh
/// @brief  Fold tree class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_Exceptions_hh
#define INCLUDED_core_kinematics_Exceptions_hh


// Unit headers
#include <core/kinematics/FoldTree.hh>

// Package Headers

// utility headers
#include <utility/excn/Exceptions.hh>

// ObjexxFCL Headers

// // C++ Headers
#include <utility/assert.hh>
#include <vector>


namespace core {
namespace kinematics {

class EXCN_InvalidFoldTree : public utility::excn::Exception {
	typedef utility::excn::Exception Parent;
public:
	EXCN_InvalidFoldTree(char const *file, int line, std::string const& msg, FoldTree f )
	: Exception(file, line, msg ), bad_tree_( f ) {};

	using utility::excn::Exception::show;

	virtual void show( std::ostream& os ) {
		os << msg() << "\nInvalid FoldTree: "<< bad_tree() << std::endl;
	}

	FoldTree const& bad_tree() { return bad_tree_; };
private:
	FoldTree bad_tree_;
};

}
}

#endif
