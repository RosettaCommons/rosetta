// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_rna_fragments_classes_fwd_hh
#define INCLUDED_protocols_rna_rna_fragments_classes_fwd_hh

// ObjexxFCL Headers
#include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/StaticIndexRange.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ Headers
#include <string>
#include <map>
#include <vector>


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Goal: to make a fragment object that can choose fragments
// "on the fly" for RNA ab inito folding.
//
// After reading in a set of torsions from, e.g., the ribosome crystal structure,
//  should be able to generate fragments of size 1, 2, or 3, with
//  exact sequence matches, partial Y/R matches, or ignoring sequence.
//
namespace protocols {
namespace farna {

	class TorsionSet;
	class FragmentLibrary;
	class RNA_Fragments;

	typedef utility::pointer::owning_ptr< FragmentLibrary > FragmentLibraryOP;

} //farna
} //protocols

#endif
