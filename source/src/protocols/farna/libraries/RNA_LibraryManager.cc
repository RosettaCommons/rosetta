// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/libraries/RNA_LibraryManager.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/libraries/RNA_LibraryManager.hh>
#include <protocols/farna/fragments/FullAtomRNA_Fragments.hh>
#include <protocols/farna/libraries/RNA_JumpLibrary.hh>
#include <protocols/farna/libraries/BasePairStepLibrary.hh>
#include <protocols/farna/options/RNA_FragmentMonteCarloOptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.libraries.RNA_LibraryManager" );

using namespace protocols::farna::options;
using namespace protocols::farna::fragments;

namespace protocols {
namespace farna {
namespace libraries {

RNA_Fragments const &
RNA_LibraryManager::rna_fragment_library( std::string const & tag ) {
	if ( rna_fragment_libraries_.find( tag ) == rna_fragment_libraries_.end() ) {
		rna_fragment_libraries_[ tag ] = RNA_FragmentsCOP( new FullAtomRNA_Fragments( tag ) );
	}
	return *rna_fragment_libraries_[ tag ];
}

RNA_JumpLibrary const &
RNA_LibraryManager::rna_jump_library( std::string const & tag ) {
	return *( rna_jump_library_cop( tag ) );
}

RNA_JumpLibraryCOP const &
RNA_LibraryManager::rna_jump_library_cop( std::string const & tag ) {
	if ( rna_jump_libraries_.find( tag ) == rna_jump_libraries_.end() ) {
		rna_jump_libraries_[ tag ] = RNA_JumpLibraryCOP( new RNA_JumpLibrary( tag ) );
	}
	return rna_jump_libraries_[ tag ];
}

RNA_JumpLibraryCOP const &
RNA_LibraryManager::rna_jump_library_cop() {
	RNA_FragmentMonteCarloOptions options; // default value is stored here.
	return rna_jump_library_cop( options.jump_library_file() );
}

BasePairStepLibrary const &
RNA_LibraryManager::canonical_base_pair_step_library() {
	if ( canonical_base_pair_step_library_ == 0 )  canonical_base_pair_step_library_ = BasePairStepLibraryCOP( new BasePairStepLibrary( true ) );
	return *canonical_base_pair_step_library_;
}

BasePairStepLibrary const &
RNA_LibraryManager::general_base_pair_step_library() {
	if ( general_base_pair_step_library_ == 0 )  general_base_pair_step_library_ = BasePairStepLibraryCOP( new BasePairStepLibrary( false ) );
	return *general_base_pair_step_library_;
}

}
} //farna
} //protocols
