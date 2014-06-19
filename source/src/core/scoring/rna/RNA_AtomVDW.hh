// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_AtomVDW.hh
/// @brief  Information on which atoms to use for computing clashes
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_AtomVDW_hh
#define INCLUDED_core_scoring_rna_RNA_AtomVDW_hh

// Unit Headers
#include <core/scoring/rna/RNA_AtomVDW.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/AA.hh>

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray4D.hh>

#include <map>

#include <utility/vector1.fwd.hh>



namespace core {
namespace scoring {
namespace rna {


class RNA_AtomVDW : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor, reads data file
	RNA_AtomVDW();

	///
	utility::vector1 < std::string > const
	vdw_atom_list( char const which_nucleotide ) const;

	Real
	bump_parameter( Size const atom1, Size const atom2,
									char const which_nucleotide1, char const which_nucleotide2 ) const;

private: //data

	//Which atoms to loop over during VDW check?
	typedef std::map< char, utility::vector1< std::string > >  AtomList;
	AtomList rna_vdw_atom_;

	//Radii for VDW calculations...
	// I originally had this as a crazy map,
	// but it was really slow... now its an FArray,
	// with a helper function to convert a,c,g,u to indices 1,2,3,4.
	ObjexxFCL::FArray4D < Real > rna_vdw_parameter_;


	};

} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
