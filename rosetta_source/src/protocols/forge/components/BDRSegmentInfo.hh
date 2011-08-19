// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/components/BDRSegmentInfo.hh
/// @brief
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_components_BDRSegmentInfo_hh
#define INCLUDED_protocols_forge_components_BDRSegmentInfo_hh

// unit headers
#include <protocols/forge/components/BDRSegmentInfo.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/forge/build/Interval.hh>

// project headers


// utility headers
#include <utility/exit.hh>

// C++ headers
#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>



namespace protocols {
namespace forge {
namespace components {


struct BDRSegmentInfo {


	typedef protocols::forge::build::Interval Interval;

	typedef std::string String;


	/// @brief default constructor
	inline
	BDRSegmentInfo() {}


	/// @brief value constructor
	inline
	BDRSegmentInfo(
		Interval const & ival,
		String const & secstruct,
		String const & aa_build = String(),
		String const & aa_dr = String()
	) :
		interval( ival ),
		ss( secstruct ),
		aa_during_build( aa_build ),
		aa_during_design_refine( aa_dr )
	{
		runtime_assert(
			aa_during_build.empty() ||
			ss.length() == core::pose::annotated_to_oneletter_sequence( aa_during_build ).length()
		);

		runtime_assert(
			aa_during_design_refine.empty() ||
			ss.length() == core::pose::annotated_to_oneletter_sequence( aa_during_design_refine ).length()
		);
	}


	// @brief default destructor
	inline
	~BDRSegmentInfo() {}


	/// @brief interval specifying the residues to remodel
	Interval interval;


	/// @brief secondary structure string, also specifies the length of the
	///  new section
	String ss;


	/// @brief annotated amino acid string, specifies the amino acid sequence
	///  of the new section during centroid build; can be empty
	/// @remarks if defined, the *one-letter* aa string must be equal in
	///  length to the secondary structure string
	String aa_during_build;


	/// @brief annotated amino acid string, specifies the amino acid sequence
	///  of the new section during design-refine; can be empty
	/// @remarks if defined, the *one-letter* aa string must be equal in
	///  length to the secondary structure string
	String aa_during_design_refine;


};


} // namespace components
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_BDR_HH */
