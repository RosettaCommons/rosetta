// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/StrandAssemblyCommon.hh
/// @brief Has most common/shared inclusions and namespace usings
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_StrandAssemblyCommon_hh
#define INCLUDED_protocols_features_strand_assembly_StrandAssemblyCommon_hh

//Basic rosetta

//C library

// c++
//#include <stdio.h>     //for remove( ) and rename( )
#include <stdlib.h> // for std::abs()

//Core

//DSSP

//External

// for get_sw_can_by_sh_id, get_central_residues_in_each_of_two_edge_strands

// for parse_my_tag

// for string return

//for vector
#include <numeric/xyzVector.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>

//Others

//Protocols

// Utility

// Utility: exception handling


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


template <typename T, size_t N> const T* mybegin(const T (&a)[N]) { return a; }
template <typename T, size_t N> const T* myend  (const T (&a)[N]) { return a+N; }
// reference: http://stackoverflow.com/questions/9874802/how-can-i-get-the-max-or-min-value-in-a-vector-c


namespace protocols {
namespace features {
namespace strand_assembly {

using core::id::NamedAtomID;
using numeric::xyzVector;


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
