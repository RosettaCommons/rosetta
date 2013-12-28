// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
//  vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_RNA_FragmentsClasses_FWD_HH
#define INCLUDED_protocols_rna_RNA_FragmentsClasses_FWD_HH

// ObjexxFCL Headers


// C++ Headers

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>



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
	typedef utility::pointer::owning_ptr< RNA_Fragments > RNA_FragmentsOP;

} //farna
} //protocols

#endif
