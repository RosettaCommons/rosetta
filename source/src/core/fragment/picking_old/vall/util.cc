// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/util.cc
/// @brief  utility functions for interacting with VallLibrary
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/util.hh>

// package headers
#include <core/fragment/picking_old/FragmentLibraryManager.hh>
#include <core/fragment/picking_old/vall/VallLibrarian.hh>
#include <core/fragment/picking_old/vall/VallLibrary.hh>
#include <core/fragment/picking_old/vall/eval/ABEGOEval.hh>
#include <core/fragment/picking_old/vall/eval/IdentityEval.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.fwd.hh>
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.fwd.hh>
#include <core/fragment/picking_old/vall/gen/LengthGen.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief pick fragments by default sec.struct IdentityScore
/// @param[in] ss secondary structure string of desired frag length
/// @param[in] top_n return the top 'n' fragments, default 200
/// @param[in] randomize add random noise within [0, 0.001) to
///  score to break up equivalent scores
/// @param[in] srfd_type The BBTorsionSRFD type to use.
core::fragment::FragDataOPs
pick_fragments_by_ss(
	std::string const & ss,
	core::Size const top_n,
	bool const randomize,
	BBTorsionSRFD const & srfd_type
)
{
	using eval::IdentityEval;
	using gen::LengthGen;
	using eval::VallFragmentEvalOP;
	using eval::VallFragmentEvalCOP;
	using gen::VallFragmentGenOP;
	using gen::VallFragmentGenCOP;

	VallLibrary const & library = FragmentLibraryManager::get_instance()->get_Vall();

	VallLibrarian librarian;
	librarian.add_fragment_gen( VallFragmentGenCOP( VallFragmentGenOP( new LengthGen( ss.length() ) ) ) );
	librarian.add_fragment_eval( VallFragmentEvalCOP( VallFragmentEvalOP( new IdentityEval( ss, 1.0, randomize ) ) ) );

	// catalog fragments
	librarian.catalog( library );

	// grab top fragments
	return librarian.top_fragments( top_n, srfd_type );
}


/// @brief pick fragments by default sec.struct IdentityScore
/// @param[in] ss secondary structure string of desired frag length
/// @param[in] aa amino acid string of same length as ss string
/// @param[in] top_n return the top 'n' fragments, default 200
/// @param[in] randomize add random noise within [0, 0.001) to
///  score to break up equivalent scores
/// @param[in] srfd_type The BBTorsionSRFD type to use.
core::fragment::FragDataOPs
pick_fragments_by_ss_plus_aa(
	std::string const & ss,
	std::string const & aa,
	core::Size const top_n,
	bool const randomize,
	BBTorsionSRFD const & srfd_type
)
{
	using eval::IdentityEval;
	using gen::LengthGen;
	using eval::VallFragmentEvalOP;
	using eval::VallFragmentEvalCOP;
	using gen::VallFragmentGenOP;
	using gen::VallFragmentGenCOP;

	VallLibrary const & library = FragmentLibraryManager::get_instance()->get_Vall();

	VallLibrarian librarian;
	librarian.add_fragment_gen( VallFragmentGenCOP( VallFragmentGenOP( new LengthGen( ss.length() ) ) ) );
	librarian.add_fragment_eval( VallFragmentEvalCOP( VallFragmentEvalOP( new IdentityEval( ss, aa, 1.0, 1.0, randomize ) ) ) );

	// catalog fragments
	librarian.catalog( library );

	// grab top fragments
	return librarian.top_fragments( top_n, srfd_type );
}


/// @brief pick fragments by default sec.struct IdentityScore
/// @param[in] ss secondary structure string of desired frag length
/// @param[in] aa amino acid string of same length as ss string
/// @param[in] abego vector of string of desired frag length
/// @param[in] top_n return the top 'n' fragments, default 200
/// @param[in] randomize add random noise within [0, 0.001) to
///  score to break up equivalent scores
/// @param[in] srfd_type The BBTorsionSRFD type to use.
core::fragment::FragDataOPs
pick_fragments(
	std::string const & ss,
	std::string const & aa,
	utility::vector1< std::string > const & abego,
	core::Size const top_n,
	bool const randomize,
	BBTorsionSRFD const & srfd_type
)
{
	using eval::IdentityEval;
	using eval::ABEGOEval;
	using gen::LengthGen;
	using eval::VallFragmentEvalOP;
	using eval::VallFragmentEvalCOP;
	using gen::VallFragmentGenOP;
	using gen::VallFragmentGenCOP;

	VallLibrary const & library = FragmentLibraryManager::get_instance()->get_Vall();

	VallLibrarian librarian;
	librarian.add_fragment_gen( VallFragmentGenCOP( VallFragmentGenOP( new LengthGen( ss.length() ) ) ) );

	if ( !aa.empty() ) {
		librarian.add_fragment_eval( VallFragmentEvalCOP( VallFragmentEvalOP( new IdentityEval( ss, aa, 1.0, 1.0, randomize ) ) ) );
	} else {
		librarian.add_fragment_eval( VallFragmentEvalCOP( VallFragmentEvalOP( new IdentityEval( ss, 1.0, randomize ) ) ) );
	}

	if ( abego.size() > 0 ) {
		librarian.add_fragment_eval( VallFragmentEvalCOP( VallFragmentEvalOP( new ABEGOEval( abego, 1.0, randomize ) ) ) );
	}

	// catalog fragments
	librarian.catalog( library );

	// grab top fragments
	return librarian.top_fragments( top_n, srfd_type );
}


} // vall
} // picking_old
} // fragment
} // core
