// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ChemicalReaction.cc
/// @brief A reaction object that keeps track of its fragment sets.
/// @author Tracy Tang (yidan.tang@vanderbilt.edu)

#include <protocols/drug_design/ChemicalReaction.hh>

#include <numeric/random/random.hh>

#include <utility/numbers.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ChemicalReaction");

//Constructor
ChemicalReaction::ChemicalReaction(
		std::string const & dir,
		std::string const & name,
		std::string const & smirks
) :
dir_( dir ),
name_( name )
{
	rxn_ = ::RDKit::ChemicalReactionOP( ::RDKit::RxnSmartsToChemicalReaction( smirks ) );
	if( ! rxn_ ) {
		TR.Error << "ERROR: Cannot parse as chemical reaction: " << smirks << std::endl;
		return;
	}
	rxn_->initReactantMatchers();
	TR << "Reaction " << name << " as parsed: " << ::RDKit::ChemicalReactionToRxnSmarts( *rxn_ ) << std::endl;
	if ( rxn_->getNumReactantTemplates() == 0 ) {
		TR.Error << "ERROR: Reaction has zero reactants." << std::endl;
		rxn_ = nullptr;
		return;
	}
	if ( rxn_->getNumProductTemplates() > 1 ) {
		TR.Error << "ERROR: Reaction has more than one product." << std::endl;
		rxn_ = nullptr;
		return;
	}
	if ( rxn_->getNumProductTemplates() == 0 ) {
		TR.Error << "ERROR: Reaction has zero products." << std::endl;
		rxn_ = nullptr;
		return;
	}
}

ChemicalReaction::ChemicalReaction(
		std::string const & name,
		std::string const & smirks
) : ChemicalReaction( "", name, smirks )
{}

ChemicalReaction::~ChemicalReaction() = default;

bool
ChemicalReaction::reaction_valid() const {
	return rxn_ != nullptr;
}

bool
ChemicalReaction::is_reaction_usable() {
	if ( rxn_ == nullptr ) { return false; }
	if ( reagents_bad_ ) { return false; } // We tried and failed to load the reagents.
	load_reagents();
	return nreagents() == reagent_smiles_.size();
}

core::Size
ChemicalReaction::nreagents() const {
	return rxn_->getNumReactantTemplates();
}

core::Size
ChemicalReaction::n_availible_reagents(core::Size reag_no) {
	load_reagents();
	if ( reag_no >= reagent_smiles_.size() ) {
		return 0; // No reagents loaded
	} else {
		return reagent_smiles_[ reag_no ].size();
	}
}

::RDKit::ROMolOP
 ChemicalReaction::reagent(core::Size reag_no, core::Size reag_index) {
	load_reagents();
	//ASSERT_ALWAYS( reag_no < reagent_smiles_.size() );
	//ASSERT_ALWAYS( 1 <= reag_index && reag_index <= reagent_smiles_[reag_no].size() );
	if ( reag_no >= reagents_.size() ) {
		reagents_.resize( reag_no + 1 ); // 0-indexed
	}
	if ( reag_index > reagents_[ reag_no ].size() ) {
		reagents_[ reag_no ].resize( reag_index );
	}
	::RDKit::RWMolOP reag = reagents_[reag_no][reag_index];
	if ( reag != nullptr ) {
		return reag;
	}

	try {
		reag = ::RDKit::RWMOL_SPTR( ::RDKit::SmilesToMol( reagent_smiles_[reag_no][reag_index] ) );
	} catch ( ::RDKit::MolSanitizeException const& ) {
		TR.Error << "PUZZLE SETUP ISSUE: issue reading SMILES string in reaction design panel " << reagent_smiles_[reag_no][reag_index] << std::endl;
		return nullptr;
	}

	if ( reag->getNumConformers() == 0 ) {
		::RDKit::DGeomHelpers::EmbedMolecule(*reag);
	}

	reagents_[reag_no][reag_index] = reag; // Cache for next time.
	return reag;
}

void
ChemicalReaction::load_reagents() {
	if ( !reagent_smiles_.empty() ) {
		return; // We already loaded them.
	}
	if ( reagents_bad_ ) {
		return; // We tried and failed already.
	}
	if ( rxn_ == nullptr || name_.empty() || dir_.empty() ) {
		return; // Bad reaction -- don't load reagents
	}

	reagent_smiles_.resize( nreagents() );

	core::Size nloaded = 0;

	for ( core::Size reag_no(0); reag_no < nreagents(); ++reag_no ) {
		std::string reagent_file = name_ + "_" + std::to_string(reag_no+1) + ".smi";

		for ( std::string const & line: load_file( dir_, reagent_file ) ) {
			utility::vector1< std::string > parts( utility::split_whitespace(line) );
			if( parts.size() == 0 ) { continue; }
			if( parts[1].size() < 1 || parts[1][0] == '#' ) { continue; }
			//if ( parts.size() > 3 ) {
			//	std::cout << "ERROR: Reaction line malformed: " << line << std::endl;
			//	reagents_bad_ = true;
			//	continue;
			//}

			++nloaded;
			reagent_smiles_[ reag_no ].push_back( parts[1] );
		}

		if ( reagent_smiles_[ reag_no ].size() == 0 ) {
			TR.Error << "ERROR: No valid reagents found for reaction " << name_ << " reagent " << reag_no << " -- skipping reaction." << std::endl;
			reagent_smiles_.clear();
			reagents_bad_ = true;
			return;
		}

	} // Iterate over reagents.

	TR << "Read " << nloaded << " reagents for reaction " << name_ << " in " << dir_ << std::endl;
}

::RDKit::RWMolOP
ChemicalReaction::reaction_result( ::RDKit::MOL_SPTR_VECT const & reagents ) const {

	std::vector< ::RDKit::MOL_SPTR_VECT > products( rxn_->runReactants( reagents ) );
	if ( products.empty() ) {
		TR.Error << "ERROR: Product list empty? Should not happen." << std::endl;
		return nullptr;
	}
//	std::cout << "Created " << products.size() << " products from reaction." << std::endl;
	core::Size prod_num( numeric::random::random_range(0,products.size()-1) ); // Perhaps should shuffle and then iterate?

	if ( products[ prod_num ].size() != 1 ) {
		TR.Error << "ERROR: Product is not a single molecule" << std::endl;
		return nullptr;
	}

	::RDKit::RWMolOP prod( new ::RDKit::RWMol( *products[ prod_num ][0] ) );

	if ( cleanup_product( *prod ) ) {
		TR.Error << "Product " << ::RDKit::MolToSmiles( *prod ) << " failed cleanup." << std::endl;
		return nullptr;
	}

	return prod;
}

bool ChemicalReaction::cleanup_product( ::RDKit::RWMol & prod ) const {
	try {
		::RDKit::MolOps::sanitizeMol(prod);
	} catch (::RDKit::MolSanitizeException const&se) {
		return true; // Cleanup failed.
	}

	// Make a basic conformation - don't put too much effort in it, as we're probably just using it for ghost atoms.
	// This long call is needed such that we can put "useExpTorsionAnglePrefs=true" and "useBasicKnowledge=true" so embedding doesn't mess up rings, etc.
	// We also use a constant seed for embedding, to make the generated conformations consistent.
	int conf_num = ::RDKit::DGeomHelpers::EmbedMolecule(prod, /*maxIterations*/ 0, /*seed*/ 111111, /*clearConfs*/ true,
			/*useRandomCoords*/ false, /*boxSizeMult*/ 2.0, /*randNegEig*/ true, /*numZeroFail*/ 1, /*coordMap*/ 0, /*optimizerForceTol*/ 1e-3,
			/*ignoreSmoothingFailures*/ false, /*enforceChirality*/ true, /*useExpTorsionAnglePrefs*/ true, /*useBasicKnowledge*/ true);
	if ( conf_num == -1 ) {
		return true; // Can't make 3D for this product.
	}

//	std::cout << "Successfully cleaned up the product." << std::endl;
	return false; // All good.
}

utility::vector1< std::string >
ChemicalReaction::load_file( std::string const & reaction_dir, std::string const & filename ) {
	utility::vector1< std::string > lines;

	std::string actual_filename( reaction_dir + "/" + filename );
	utility::io::izstream file( actual_filename );
	if ( ! file.good() ) {
		actual_filename = basic::database::full_name("protocol_data/drug_design/"+filename);
		utility::io::izstream data( actual_filename );
		if ( ! data.good() ) {
			TR.Error << "Can't find reaction file, either (./)" << filename << " or " << actual_filename << std::endl;
			utility_exit_with_message("ERROR: cannot find reaction file.");
		}
	}

	std::string line;
	while( getline( file, line ) ) {
		lines.push_back( line );
	}
	return lines;
}

::RDKit::ChemicalReactionOP const &
 ChemicalReaction::get_reaction() const {
	return rxn_;
}


}
}

