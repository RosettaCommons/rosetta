// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/min.hh
/// @brief  Kinematics
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_AtomNode_hh
#define INCLUDED_core_optimization_AtomNode_hh


// Package headers

// Project headers
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID_Map.hh>

// Rosetta headers
// #include <util_basic.hh>
// #include <jump_classes.hh>
// #include <core/kinematics/Stub.hh>
// #include <id.hh>

// // Numeric headers
// #include <numeric/all.fwd.hh>
// #include <numeric/conversions.hh>
// #include <numeric/xyzMatrix.hh>
// #include <numeric/xyzVector.hh>

// // ObjexxFCL headers
// #include <ObjexxFCL/FArray1D.hh>

// // Utility headers
// #include <utility/io/all.fwd.hh>

// // C++ headers
// #include <algorithm>
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// //#include <iosfwd>
// #include <utility/assert.hh>
// #include <vector>
// #include <string>
// #include <map>
// #include <list>


namespace core {
namespace optimization {


class AtomNode
{
public:
	typedef id::AtomID AtomID;


public:
	AtomNode( AtomID const & id_in ):id( id_in ) {}

	inline
	int
	rsd() const { return id.rsd(); }

	inline
	int
	atomno() const { return id.atomno(); }

	inline
	std::vector< AtomNode* >::iterator
	nbr_list_begin() { return nbrs.begin(); }

	inline
	std::vector< AtomNode* >::iterator
	nbr_list_end() { return nbrs.end(); }

	inline
	void
	add_nbr( AtomNode* nbr ) { nbrs.push_back( nbr ); }

private:
	AtomID id;
	std::vector< AtomNode* > nbrs;

}; // AtomNode


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_min_HH
