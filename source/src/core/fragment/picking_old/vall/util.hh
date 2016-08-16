// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/util.hh
/// @brief  utility functions for interacting with VallLibrary
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_util_hh
#define INCLUDED_core_fragment_picking_old_vall_util_hh

// project headers
#include <core/fragment/FragData.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <utility/vector1.hh>


// C++ headers


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief convert a fragment extent to a FragData instance
/// @param[in] begin Iterator pointing the beginning of the VallResidue extent.
/// @param[in] end Iterator pointing just past the end of the VallResidue extent.
/// @param[in] srfd_type The type of BBTorsionSRFD to use.
template< typename VallResidueIterator >
AnnotatedFragDataOP extent_to_fragdata(
	VallResidueIterator begin,
	VallResidueIterator end,
	BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
) {
	AnnotatedFragDataOP fragdata( new AnnotatedFragData( begin->id(), begin->resi() ) );

	for ( VallResidueIterator r = begin; r != end; ++r ) {
		fragdata->add_residue( r->bbtorsion_srfd( srfd_type ) );
	}

	if ( fragdata->size() > 0 ) {
		fragdata->set_valid();
	}

	return fragdata;
}

/// @brief adds score to fragdata object
//
template< typename VallResidueIterator >
AnnotatedFragDataOP extent_to_fragdata(
	VallResidueIterator begin,
	VallResidueIterator end,
	core::Real score,
	BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
) {
	AnnotatedFragDataOP fragdata = extent_to_fragdata(begin,end,srfd_type);
	fragdata->set_score(score);
	return fragdata;
}

/// @brief extract amino acid sequence from a fragment extent into a string
template< typename VallResidueIterator >
std::string extent_aa_str( VallResidueIterator begin, VallResidueIterator end ) {
	std::string aa_str;

	for ( VallResidueIterator r = begin; r != end; ++r ) {
		aa_str.push_back( r->aa() );
	}

	return aa_str;
}


/// @brief extract secondary structure from a fragment extent into a string
template< typename VallResidueIterator >
std::string extent_ss_str( VallResidueIterator begin, VallResidueIterator end ) {
	std::string ss_str;

	for ( VallResidueIterator r = begin; r != end; ++r ) {
		ss_str.push_back( r->ss() );
	}

	return ss_str;
}


/// @brief pick fragments by default sec.struct IdentityScore
/// @param[in] ss secondary structure string of desired frag length
/// @param[in] top_n return the top 'n' fragments, default 200
/// @param[in] randomize add random noise within [0, 0.001) to
/// @param[in] srfd_type The BBTorsionSRFD type to use.
///  score to break up equivalent scores
core::fragment::FragDataOPs
pick_fragments_by_ss(
	std::string const & ss,
	core::Size const top_n,
	bool const randomize = true,
	BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
);


/// @brief pick fragments by default sec.struct IdentityScore
/// @param[in] ss secondary structure string of desired frag length
/// @param[in] aa amino acid string of same length as ss string
/// @param[in] top_n return the top 'n' fragments, default 200
/// @param[in] randomize add random noise within [0, 0.001) to
/// @param[in] srfd_type The BBTorsionSRFD type to use.
///  score to break up equivalent scores
core::fragment::FragDataOPs
pick_fragments_by_ss_plus_aa(
	std::string const & ss,
	std::string const & aa,
	core::Size const top_n,
	bool const randomize = true,
	BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
);


/// @brief pick fragments by default sec.struct IdentityScore
/// @param[in] ss secondary structure string of desired frag length
/// @param[in] aa amino acid string of same length as ss string
/// @param[in] abego vector of string of desired frag length
/// @param[in] top_n return the top 'n' fragments, default 200
/// @param[in] randomize add random noise within [0, 0.001) to
/// @param[in] srfd_type The BBTorsionSRFD type to use.
///  score to break up equivalent scores
core::fragment::FragDataOPs
pick_fragments(
	std::string const & ss,
	std::string const & aa,
	utility::vector1< std::string > const & abego,
	core::Size const top_n = 200,
	bool const randomize = true,
	BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
);


} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_util_HH */
