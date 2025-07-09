// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionBasedAnalogSampler.hh
/// @brief Reaction-based and similarity guided analoging in a given chemical library
/// @author Yidan Tang (yidan.tang@vanderbilt.edu)

#include <protocols/drug_design/ReactionBasedAnalogSampler.hh>
#include <protocols/drug_design/ReactionBasedAnalogSamplerCreator.hh>
#include <protocols/drug_design/ReactionChemistry.hh>

#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/chemistries/util.hh>

//  The following are needed for ResidueType::operator=
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/Orbital.hh>
#include <core/chemical/ResidueConnection.hh>
// end operator= includes

#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>

#include <rdkit/DataStructs/BitOps.h>  // For similarity bit wrappers
#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h> // For MolToSmarts
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/Fingerprints/MorganFingerprints.h>
#include <rdkit/GraphMol/Fingerprints/Fingerprints.h>
#include <rdkit/GraphMol/MolAlign/O3AAlignMolecules.h>
#include <rdkit/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>

#include <list>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <chrono>

// debuggg
#include <iostream>
#include <fstream>
#include <queue>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ReactionBasedAnalogSampler");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
ReactionBasedAnalogSamplerCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new ReactionBasedAnalogSampler );
}

std::string
ReactionBasedAnalogSamplerCreator::keyname() const {
	return ReactionBasedAnalogSampler::class_name();
}

void
ReactionBasedAnalogSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ReactionBasedAnalogSampler::provide_xml_schema( xsd );
}

//------------------------- Chemistry -----------------------------

ReactionBasedAnalogSampler::ReactionBasedAnalogSampler():
	Chemistry(class_name()),
	geo_spl_ratio_( 0.1 ),
	dynamic_sampling_( false ),
	minCandidates_( 20 )
{}

ReactionBasedAnalogSampler::~ReactionBasedAnalogSampler() {}

/// @brief load reaction files from a designated path (deprecated)
void
ReactionBasedAnalogSampler::load_reactions( std::string const & reaction_dir, std::string const & filename ) {

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
	while ( getline( file, line ) ) {
		//load all the reactions that are available
		//set up for file is 1st column = reaction name, 2nd column = smirk string
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if ( parts.size() == 0 ) { continue; }
		if ( parts[1].size() < 1 || parts[1][0] == '#' ) { continue; }
		if ( parts.size() != 2 ) {
			std::cout << "ERROR: Reaction line malformed: " << line << std::endl;
			continue;
		}
		std::cout << "Parsing Reaction '" << parts[1] << ": " << parts[2] << std::endl;
		ChemicalReactionOP rxn( new ChemicalReaction(reaction_dir, parts[1], parts[2] ) );
		if ( rxn->reaction_valid() ) {
			rxns_[ parts[1] ] = rxn;
		}
	}
	std::cout << "Total Reactions loaded: " << rxns_.size() << std::endl;
}

/// @brief load reactions from a single file
void
ReactionBasedAnalogSampler::load_reactions( std::string const & filename ) {
	utility::io::izstream file( filename );
	if ( ! file.good() ) {
		std::string actual_filename = basic::database::full_name("protocol_data/drug_design/"+filename);
		utility::io::izstream data( actual_filename );
		if ( ! data.good() ) {
			TR.Error << "Can't find reaction file, either (./)" << filename << " or " << actual_filename << std::endl;
			utility_exit_with_message("ERROR: cannot find reaction file.");
		}
	}

	std::string line;
	getline( file, line ); //first line is column names
	while ( getline( file, line ) ) {
		//set up for file is 1st column = reaction id, 2nd column = No. of components, 3rd column = reaction smirks
		//2nd column (no. of components) is not used here since RDKit can detect this through reaction SMARTS. This is simply to line up with the format (used in REvoLd)
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if ( parts.size() == 0 ) { continue; }
		if ( parts.size() < 3 ) {
			TR.Error << "ERROR: Reaction line malformed: " << line << std::endl;
			continue;
		}
		std::cout << "Parsing Reaction '" << parts[1] << ": " << parts[3] << std::endl;
		ChemicalReactionOP rxn( new ChemicalReaction( parts[1], parts[3] ) );
		if ( rxn->reaction_valid() ) {
			rxns_[ parts[1] ] = rxn;
		}
	}
	TR << "Total Reactions loaded: " << rxns_.size() << std::endl;
}

/// @brief load reagents through corresponding reactions (deprecated)
void
ReactionBasedAnalogSampler::load_all_reagents() {
	for ( auto it : rxns_ ) {
		ChemicalReactionOP rxn( it.second );
		rxn->load_reagents();

		for ( core::Size rr(1); rr <= rxn->nreagents(); ++rr ) {
			for ( core::Size ri(1); ri <= rxn->n_availible_reagents(rr); ++ri ) {
				struct Reagent reagent;
				// Generate fingerprint for each availible reagent
				::RDKit::ROMolOP reag = rxn->reagent( rr, ri );
				reagent.rdmol_ = reag;
				reagent.fp_ = std::shared_ptr<ExplicitBitVect>( ::RDKit::RDKFingerprintMol(*reag) );
				reagent.rxn_ = rxn->reaction_name();
				reagent.no_ = rr;
				reagents_.push_back( reagent );
			}
		}
	}
	geometric_sampling( reagent_sampler_, reagents_.size() );
	TR << "There are in total " << reagents_.size() << " reagents." << std::endl;
}

/// @brief load reagents from a single file
void
ReactionBasedAnalogSampler::load_all_reagents( std::string const & filename ) {
	utility::io::izstream file( filename );
	std::string line;
	getline( file, line ); //first line is column names
	while ( getline( file, line ) ) {
		//set up for file is 1st column = SMILES, 2nd column = reactant id, 3rd column = component #, 4th column = reaction id
		//again, 2nd column (reactant id) is not used here. This is simply to line up with the format (used in REvoLd)
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if ( parts.size() == 0 ) { continue; }
		if ( parts.size() < 4 ) {
			TR.Error << "ERROR: Reactant line malformed: " << line << std::endl;
			continue;
		}
		struct Reagent reagent;
		// Generate fingerprint for each availible reagent
		::RDKit::ROMolOP reag = ::RDKit::RWMOL_SPTR( ::RDKit::SmilesToMol( parts[1] ) );
		if ( check_reagent_validity( reag ) ) {   // Remove any invalid or problematic structures from the chemical database
			reagent.rdmol_ = reag;
			reagent.fp_ = std::shared_ptr<ExplicitBitVect>( ::RDKit::RDKFingerprintMol(*reag) );
			reagent.no_ = (core::Size)std::stoi(parts[3]);
			reagent.rxn_ = parts[4];
			reagents_.push_back( reagent );
		}
	}
	TR << "Total reagents loaded: " << reagents_.size() << std::endl;
	geometric_sampling( reagent_sampler_, reagents_.size() );
}

/// @brief check if the reagent contains bad structures that would later fail RDKit
bool
ReactionBasedAnalogSampler::check_reagent_validity( ::RDKit::ROMOL_SPTR reag ) const {
	// A collection of SMARTS that will fail RDKit later on
	static const utility::vector1< ::RDKit::ROMOL_SPTR > BAD_SMARTS {
		::RDKit::RWMolOP(::RDKit::SmartsToMol("[*](S(F)(F)(F)(F)F)"))
		};
	::RDKit::MatchVectType matchvect;

	for ( ::RDKit::ROMOL_SPTR const & bad_molOP : BAD_SMARTS ) {
		if ( ::RDKit::SubstructMatch( *reag, *bad_molOP, matchvect ) ) {
			return false;
		}
	}
	return true;
}

void
ReactionBasedAnalogSampler::apply( core::chemical::MutableResidueType & rsdtype )
{
	core::chemical::rdkit::RestypeToRDMol to_converter(rsdtype); // Neutralize and remove hydrogens.
	::RDKit::RWMOL_SPTR rdmol( to_converter.Mol() ); // Convert
	core::chemical::rdkit::label_with_index(*rdmol, "Original_Index"); // (Re)Label with "Original_Index" to use quick mapping behavior.
	// Properties carry through the reaction, except for a "bug" with the explicitly mentioned atoms.

	// Generate fingerprint for the input residue molecule
	std::shared_ptr<ExplicitBitVect> fp_rdmol( ::RDKit::RDKFingerprintMol(*rdmol) );

	if ( rxns_.size() == 0 || reagents_.size() == 0 ) {
		TR.Warning << "No reactions found or reagents are empty. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	// container of pairs: first term is similarity score, second term is original index of the fragment (to keep track after sorting)
	std::list< std::pair< core::Real, core::Size > > score_idx_set;

	// Generate similarity score and reagent index pairs for all reagents
	for ( core::Size ri(1); ri <= reagents_.size(); ++ri ) {

		std::shared_ptr<ExplicitBitVect> fp_reag = reagents_[ri].fp_;

		// Compare to the entire molecule and generate a similarity score as the weight
		core::Real score( TverskySimilarity( *fp_rdmol, *fp_reag, 0.1, 0.9 ) );
		score_idx_set.push_front( std::make_pair( score, ri ) );
	}

	// Sort the pairs so that the most similar fragments to the input are on the top of the list
	score_idx_set.sort(std::greater< std::pair< core::Real, core::Size > >());

	// If use dynamic sampling, update the geometric sampling ratio for current cycle
	if ( dynamic_sampling_ ) {
		// If OFF_after_n step is set a valid MC cycle, set the sampling to base ratio and turn off dynamic sampling
		if ( dynamic_spl_ratios_[4] > 0 && ( dynamic_spl_ratios_[4] <= visited_.size() ) ) {
			geo_spl_ratio_ = dynamic_spl_ratios_[5];
			dynamic_sampling_ = false;
			TR << "Dynamic sampling turned off." << std::endl;
		} else {
			std::string curr_smiles = ::RDKit::MolToSmiles( *rdmol );
			if ( last_smiles_ != "" && last_smiles_ != curr_smiles ) {  // last MC cycle accepted
				reset_spl_ratio();
				TR << "Last MC accepted. Reset sampling ratio to " << dynamic_spl_ratios_[1] << std::endl;
			} else {              // last MC cycle rejected
				geo_spl_ratio_ = std::min( dynamic_spl_ratios_[2], geo_spl_ratio_ + dynamic_spl_ratios_[3]);
				TR << "Last MC rejected. Sampling ratio set to " << geo_spl_ratio_ << std::endl;
			}
			geometric_sampling( reagent_sampler_, reagents_.size() - visited_.size() );
			last_smiles_ = curr_smiles;
		}
	}

	// Sample from the fragments and propose a small list of candidates
	utility::vector1< Product > candidates = sample( score_idx_set );

	if ( candidates.empty() ) { // If no candidate found from last step
		TR.Error << "ERROR: No candidate product yielded from sampling.";
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	Product product = sample_candidate( candidates, rdmol );

	visited_.insert( product );
	::RDKit::RWMolOP prod = product.rdmol_;
	utility::vector1< core::Size > prod_frags = product.frags_;

	// // Search for analog being deciding the final output product (deprecated)
	// utility::vector1< Product > analogs = analog_search( product );
	// Product ana = sample_candidate( analogs, rdmol );
	// visited_.insert( ana );
	// ::RDKit::RWMolOP prod = ana.rdmol_;
	// utility::vector1< core::Size > prod_frags = ana.frags_;
	// TR << "Analog made: " << ana.smiles_ << std::endl;

	// Find a mapping from the pre-reaction molecule to the post-reaction molecule. Possibly empty mapping at this time.
	// Hopefully, the "Original_Index" property carries through, and can be used as a quick shortcut.
	core::chemical::IndexIndexMapping rxn_map( find_O3A_mapping( rdmol, prod ) );

	if ( rxn_map.empty() ) {
		TR.Warning << "Cannot find substructure mapping of reactant to product." << std::endl;
	} else {
		TR << "O3A mapping completed." << std::endl;
	}

	// Now convert the residue into a Rosetta residue type.

	core::chemical::VDIndexMapping restype_prod_map( combine( to_converter.vd_to_index(), rxn_map ) );
	core::chemical::rdkit::RDMolToRestype from_converter(*prod);
	TR << "Product converting to Rosetta residue..." << std::endl;

	// Check if re-neighboring is necessary

	core::chemical::VD rsd_nbr ( rsdtype.nbr_vertex() );

	// First try regenerate atom mapping with more conformers for product
	if ( restype_prod_map[ rsd_nbr ] == restype_prod_map.invalid_entry() ) {
		TR << "Cannot map neighbor atom to product. Generating more conformers for product and redo O3A..." << std::endl;
		::RDKit::DGeomHelpers::EmbedMultipleConfs( *prod, 50, ::RDKit::DGeomHelpers::srETKDGv3 );
		rxn_map = find_O3A_mapping( rdmol, prod );
		restype_prod_map = combine( to_converter.vd_to_index(), rxn_map );
	} else {
		TR << "No need for re-neighboring." << std::endl;
	}

	// If still need reneighboring
	if ( restype_prod_map[ rsd_nbr ] == restype_prod_map.invalid_entry() ) {
		TR << "Still cannot map neighbor atom to product. Re-assigning neighbor atom..." << std::endl; // Previous neighbor atom not in common substructure
		std::queue< core::chemical::VD > nbr_queue;
		nbr_queue.push( rsdtype.nbr_vertex() );

		core::chemical::VD found_nbr = core::chemical::INVALID_VD;
		while ( !nbr_queue.empty() ) {
			for ( core::chemical::VD vd: rsdtype.bonded_heavyatoms( nbr_queue.front() ) ) {
				if ( restype_prod_map[ vd ] != restype_prod_map.invalid_entry() ) {
					found_nbr = vd;
					break;
				} else {
					nbr_queue.push( vd );
				}
			}
			nbr_queue.pop(); // Because there's no combined pop-and-return-value.

			if ( found_nbr != core::chemical::INVALID_VD ) { break; }
		}

		if ( found_nbr == core::chemical::INVALID_VD ) {
			utility_exit_with_message( "Unable to find acceptable neighbor atom." );
		}

		TR << "Assigning neighbor atom to '" << rsdtype.atom_name( found_nbr ) << "'" << std::endl;

		from_converter.set_nbr( restype_prod_map[ found_nbr ] );

	} else {
		from_converter.set_nbr( restype_prod_map[ rsdtype.nbr_vertex() ] );
		TR << "New residue neighbor set." << std::endl;
	}

	core::chemical::MutableResidueTypeOP new_resop( from_converter.generate_restype(rsdtype,restype_prod_map) );
	mapping_ = combine( restype_prod_map, from_converter.index_to_vd() );
	rsdtype = *new_resop;
	TR << "New residue conversion completed." << std::endl;

	// Add each chosen fragment as a string property to the residue.
	for ( core::Size rr(1); rr <= prod_frags.size(); ++rr ) {
		core::Size i = prod_frags[ rr ];
		::RDKit::ROMolOP reag = reagents_[ i ].rdmol_;
		rsdtype.add_string_property( "fragment" + std::to_string( rr ), ::RDKit::MolToSmiles( *reag ) );
		if ( rsdtype.properties().string_properties().count( "fragment" + std::to_string( rr ) ) > 0 ) {
			TR << "Fragment property: " + ::RDKit::MolToSmiles( *reag ) + " added." << std::endl;
		} else {
			TR << "Failed to add fragment property: " << ::RDKit::MolToSmiles( *reag ) << std::endl;
		}
	}
	// Add the product SMILES as a string property to the residue.
	rsdtype.add_string_property( "SMILES", product.smiles_ );

	mapping_ = combine( mapping_, combine( core::chemical::VDStringMapping(*new_resop), core::chemical::StringVDMapping(rsdtype)) );

	// That's it, we're successful
	set_last_status( core::chemical::modifications::SUCCESS );
	return;
}

std::shared_ptr<::RDKit::SparseIntVect<unsigned int>>
ReactionBasedAnalogSampler::getMorganFingerprint( ::RDKit::ROMol const & mol, bool useFeatures /*=false*/ ) const {
	// Fingerprint properties
	int radius = 2;

	if ( useFeatures ) { // FCFP
		std::vector<unsigned int> invars( mol.getNumAtoms() );
		::RDKit::MorganFingerprints::getFeatureInvariants( mol, invars );
		return std::shared_ptr<::RDKit::SparseIntVect<unsigned int>>( ::RDKit::MorganFingerprints::getFingerprint( mol, radius, &invars ) );

	} else { // ECFP
		return std::shared_ptr<::RDKit::SparseIntVect<unsigned int>>( ::RDKit::MorganFingerprints::getFingerprint( mol, radius ) );
	}
}

core::chemical::IndexIndexMapping
ReactionBasedAnalogSampler::find_O3A_mapping( ::RDKit::ROMOL_SPTR from, ::RDKit::ROMOL_SPTR to ) const {
	core::chemical::IndexIndexMapping retval( core::Size(-1), core::Size(-1) );
	// Make sure both molecules have at least 1 conformer
	if ( from->getNumConformers() > 0 ) {
		if ( to->getNumConformers() > 0 ) {

			std::vector< boost::shared_ptr< ::RDKit::MolAlign::O3A >> O3As( to->getNumConformers() ); // holder for O3A objects
			::RDKit::ROMolOP fromCopy( new ::RDKit::ROMol(*from) );
			::RDKit::MMFF::MMFFMolProperties fromMP( *fromCopy );
			::RDKit::ROMolOP toCopy( new ::RDKit::ROMol(*to) );
			::RDKit::MMFF::MMFFMolProperties toMP( *toCopy );

			::RDKit::MolAlign::getO3AForProbeConfs( *to, *from, &toMP, &fromMP, O3As );

			// Find the conformer that aligns best
			double max_score( 0 );
			std::vector< boost::shared_ptr< ::RDKit::MolAlign::O3A >>::size_type max_i( 0 );
			for ( std::vector< boost::shared_ptr< ::RDKit::MolAlign::O3A >>::size_type oi(0); oi < O3As.size(); ++oi ) {
				boost::shared_ptr< ::RDKit::MolAlign::O3A > o3a = O3As[ oi ];
				o3a->align();
				double score = o3a->score();
				if ( score > max_score ) {
					max_score = score;
					max_i = oi;
				}
			}
			// Get best conformer's atom mapping
			retval = core::chemical::rdkit::convert_match_vect_to_index_index_map( *(O3As[ max_i ]->matches()) ).reverse();

		} else {
			TR.Error << "New product has invalid number of conformers.";
		}
	} else {
		TR.Error << "Previous product has invalid number of conformers.";
	}

	return retval;
}

/// @brief Draw samples given a list of fragment similarity to reference input
utility::vector1< ReactionBasedAnalogSampler::Product >
ReactionBasedAnalogSampler::sample( std::list< std::pair< core::Real, core::Size > > & score_idx_set ) const {
	// Setup the sampling set where each reaction id has its own bucket of reagent indices
	std::unordered_map< std::string, utility::vector1< utility::vector1< core::Size > > > sample_sets;
	utility::vector1< Product > candidate_set;  // Setup the candidate set

	int revisited = 0;
	int nLoop = 0;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	while ( candidate_set.size() < minCandidates_ ) {
		nLoop++;
		core::Size choice = reagent_sampler_.random_sample() - 1; //because std::list is 0-index
		std::list< std::pair< core::Real, core::Size > >::iterator itr = score_idx_set.begin();
		std::advance( itr, choice );
		core::Size reag_idx = ( *itr ).second;
		Reagent sample = reagents_[ reag_idx ];
		std::string sample_rxn = sample.rxn_;
		core::Size sample_no = sample.no_;
		// Initialize the reaction bucket if first time sampling the reaction
		if ( sample_sets.find( sample_rxn ) == sample_sets.end() ) {
			core::Size rxn_nComponent = rxns_.at( sample_rxn )->nreagents();
			utility::vector1< utility::vector1< core::Size > > rxn_set( rxn_nComponent, utility::vector1< core::Size >() );
			sample_sets[ sample_rxn ] = rxn_set;
		}
		sample_sets[ sample_rxn ][ sample_no ].push_back( reag_idx );
		// Generate all products for current sample's reaction and add to candidate set
		pair( sample_no, sample_sets[ sample_rxn ], candidate_set, revisited );
		// Remove the element so that sampler don't re-pick the same on the next cycle
		score_idx_set.erase( itr );
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	TR << "Sampling done in " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0  << " sec, sampled " << nLoop << " fragments, and generated " << candidate_set.size() << " candidates." << std::endl;

	return candidate_set;
}

utility::vector1< ReactionBasedAnalogSampler::Product >
ReactionBasedAnalogSampler::sample_fragment( utility::vector1< std::list< std::pair< core::Real, core::Size > > > & score_idx_set ) const {
	// Setup the sampling set where each reaction id has its own bucket of reagent indices
	std::unordered_map< std::string, utility::vector1< utility::vector1< core::Size > > > sample_sets;
	utility::vector1< Product > candidate_set;   // Setup the candidate set
	std::unordered_set< core::Size > sampled_reagents; // Avoid resampling the same reagent

	int revisited = 0;
	int nLoop = 0;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	while ( candidate_set.size() < minCandidates_ ) {
		nLoop++;
		for ( std::list< std::pair< core::Real, core::Size > > & frag_set : score_idx_set ) {
			core::Size choice = reagent_sampler_.random_sample() - 1; // because std::list is 0-index
			std::list< std::pair< core::Real, core::Size > >::iterator itr = frag_set.begin();
			std::advance( itr, choice );
			core::Size reag_idx = ( *itr ).second;
			// Remove the element so that sampler don't re-pick the same on the next cycle
			frag_set.erase( itr );
			if ( sampled_reagents.find( reag_idx ) == sampled_reagents.end() ) {
				Reagent sample = reagents_[ reag_idx ];
				std::string sample_rxn = sample.rxn_;
				core::Size sample_no = sample.no_;
				// Initialize the reaction bucket if first time sampling the reaction
				if ( sample_sets.find( sample_rxn ) == sample_sets.end() ) {
					core::Size rxn_nComponent = rxns_.at( sample_rxn )->nreagents();
					utility::vector1< utility::vector1< core::Size > > rxn_set( rxn_nComponent, utility::vector1< core::Size >() );
					sample_sets[ sample_rxn ] = rxn_set;
				}
				sample_sets[ sample_rxn ][ sample_no ].push_back( reag_idx );
				// Generate all products for current sample's reaction and add to candidate set
				pair( sample_no, sample_sets[ sample_rxn ], candidate_set, revisited );
				sampled_reagents.insert( reag_idx );
			}
		}
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	TR << "Sampling by Fragment done in " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0  << " sec and looped " << nLoop << std::endl;
	return candidate_set;
}

utility::vector1< ReactionBasedAnalogSampler::Product >
ReactionBasedAnalogSampler::analog_search( const ReactionBasedAnalogSampler::Product & prod ) const {
	utility::vector1< core::Size > frags = prod.frags_;
	utility::vector1< std::shared_ptr<ExplicitBitVect> > frags_fp;
	utility::vector1< core::Size > ana_frags;
	utility::vector1< std::list< std::pair< core::Real, core::Size > > > score_idx_analog_set;
	// Initialization
	for ( core::Size const & frag_i : frags ) {
		std::list< std::pair< core::Real, core::Size > > score_idx_frag_set;
		score_idx_analog_set.push_back( score_idx_frag_set );
		frags_fp.push_back( reagents_[ frag_i ].fp_ );
	}

	for ( core::Size ri(1); ri <= reagents_.size(); ++ri ) {
		ReactionBasedAnalogSampler::Reagent reag = reagents_[ ri ];
		// pattern fingerprint seems better suited for analoging than ECFP since it is substructure-based
		std::shared_ptr<ExplicitBitVect> fp_reag = reag.fp_;
		for ( core::Size fi(1); fi <= frags.size(); ++fi ) {
			std::shared_ptr<ExplicitBitVect> fp_sample = frags_fp[ fi ];
			// Compare to the chosen reagent and generate a similarity score
			core::Real score( TanimotoSimilarity( *fp_sample, *fp_reag ) );
			score_idx_analog_set[ fi ].push_front( std::make_pair( score, ri ) );
		}
	}

	for ( core::Size fi(1); fi <= frags.size(); ++fi ) {
		// Sort the score_index set before sampling
		score_idx_analog_set[ fi ].sort(std::greater< std::pair< core::Real, core::Size > >());
	}

	// Sample from the fragments and propose a small list of candidates
	utility::vector1< Product > candidates = sample_fragment( score_idx_analog_set );
	return candidates;
}

/// @brief Draw a product from the given candidate set
ReactionBasedAnalogSampler::Product
ReactionBasedAnalogSampler::sample_candidate( const utility::vector1< Product > & candidate_set, const ::RDKit::RWMOL_SPTR rdmol ) {
	// container of pairs for candidate sampling: first term is similarity score, second term is original index of the candidate (to keep track after sorting)
	utility::vector1< std::pair< core::Real, core::Size > > score_idx_candidate_set;
	numeric::random::WeightedSampler product_sampler;
	std::shared_ptr<ExplicitBitVect> fp_rdmol( ::RDKit::RDKFingerprintMol(*rdmol) );

	// Compare each candidate to the input query molecule and rank the set based on similarity

	core::Real totalSim = 0;
	core::Real minSim = 1;
	core::Real maxSim = 0;
	for ( core::Size i(1); i <= candidate_set.size(); ++i ) {
		core::Real sim = TanimotoSimilarity( *fp_rdmol, *(candidate_set[i].fp_) );
		//  core::Real sim = TanimotoSimilarity( *ecfp_rdmol, *(candidate_set[i].ecfp_) );
		score_idx_candidate_set.push_back( std::make_pair( sim, i ) );
		totalSim += sim;
		if ( sim < minSim ) { minSim = sim; }
		if ( sim > maxSim ) { maxSim = sim; }
	}

	std::sort( score_idx_candidate_set.begin(), score_idx_candidate_set.end(), sortbySim );

	core::Real meanSim = totalSim / ( candidate_set.size() );
	core::Real medianSim = score_idx_candidate_set[ candidate_set.size() / 2 ].first;

	TR << "Candidate set statistics: similarity to current input [ ";
	TR << boost::format("min = %.2f ") % minSim;
	TR << boost::format("max = %.2f ") % maxSim;
	TR << boost::format("mean = %.2f ") % meanSim;
	TR << boost::format("median = %.2f ]") % medianSim << std::endl;

	// Pick a product and update the visited collector
	geometric_sampling( product_sampler, candidate_set.size() );
	core::Size choice = product_sampler.random_sample();
	core::Size cand_idx = score_idx_candidate_set[ choice ].second;
	Product prod = candidate_set[ cand_idx ];
	TR << "Product selected: " << prod.smiles_ << " with " << score_idx_candidate_set[ choice ].first << " similarity to current input." << std::endl;
	return prod;
}

bool
ReactionBasedAnalogSampler::sortbySim( const std::pair< core::Real, core::Size >& a, const std::pair< core::Real, core::Size >& b ) {
	return (a.first > b.first);
}

/// @brief Find all pair of reagents that can undergo the specific reaction from the current sets.
void
ReactionBasedAnalogSampler::pair(
	core::Size r_no,
	utility::vector1< utility::vector1< core::Size > > const & sets,
	utility::vector1< ReactionBasedAnalogSampler::Product > & candidates,
	int & revisits
) const {
	utility::vector1< utility::vector1< core::Size > > pairs;

	// Make sure each reactant bucket is not empty. Return if any bucket is empty
	for ( utility::vector1< core::Size > const & s : sets ) {
		if ( s.empty() ) {
			return;
		}
	}
	// Run a DFS to generate all possible pairs, starting from the first component
	pair( 1, r_no, utility::vector1< core::Size >(), pairs, sets );

	for ( utility::vector1< core::Size > & p: pairs ) {
		Product prod = run_reaction( p );
		// Check if the product is nullptr
		if ( prod.rdmol_ == nullptr ) {
			continue;
		}
		// Check if the product has been visited in previous MC cycles
		if ( visited_.find( prod ) == visited_.end() ) {
			candidates.push_back( prod );
		} else {
			revisits++;
		}
	}
}

/// @brief Helper function for pairing. Use a DFS to generate all pairs.
void
ReactionBasedAnalogSampler::pair(
	core::Size curr_no,            // current layer (component No.) of DFS
	core::Size r_no,             // component No. to locate the current sample
	utility::vector1< core::Size > single_pair,      // current pair to work on
	utility::vector1< utility::vector1< core::Size > > & all_pairs, // collection of all pairs
	utility::vector1< utility::vector1< core::Size > > const & sets // the current reaction bucket
) const {
	if ( curr_no > sets.size() ) {
		all_pairs.push_back( single_pair );
		return;
	}
	utility::vector1< core::Size > reagent_set = sets[ curr_no ];
	if ( curr_no == r_no ) {
		core::Size reag_idx = reagent_set.back();
		single_pair.push_back( reag_idx );
		pair( curr_no + 1, r_no, single_pair, all_pairs, sets );
	} else {
		utility::vector1< core::Size > temp_pair = single_pair;
		for ( core::Size const & reag_idx : reagent_set ) {
			single_pair.push_back( reag_idx );
			pair( curr_no + 1, r_no, single_pair, all_pairs, sets );
			single_pair = temp_pair;
		}
	}
}

ReactionBasedAnalogSampler::Product
ReactionBasedAnalogSampler::run_reaction( utility::vector1< core::Size > const & reagents ) const {
	std::string rxn_name = reagents_[ reagents[1] ].rxn_;
	utility::vector1< Reagent > reagent_selection;
	for ( core::Size i : reagents ) {
		reagent_selection.push_back( reagents_[i] );
	}
	std::sort( reagent_selection.begin(), reagent_selection.end() );
	::RDKit::MOL_SPTR_VECT rdreagents;
	for ( Reagent r : reagent_selection ) {
		::RDKit::ROMolOP reag = r.rdmol_;
		if ( reag == nullptr ) {
			TR.Error << "ERROR: Can't run reaction, points to bad reagent.";
		}
		rdreagents.push_back( reag );
	}
	ChemicalReactionOP rxn = rxns_.at( rxn_name );
	::RDKit::RWMolOP prod = rxn->reaction_result( rdreagents );

	struct Product p;

	if ( prod == nullptr ) {
		p.rdmol_ = nullptr;
		return p;
	}

	try {
		// Attempt to clean up silly charges and perform internal sanitization.
		core::chemical::rdkit::neutralize_rdmol( *prod, false );
	} catch ( ::RDKit::MolSanitizeException &se ) {
		TR.Error << "Cannot Sanitize product with RDKit. Skipping: " << ::RDKit::MolToSmiles( *prod ) << std::endl;
		p.rdmol_ = nullptr;
		return p;
	}

	p.rdmol_ = prod;
	p.smiles_ = ::RDKit::MolToSmiles( *prod );
	p.fp_ = std::shared_ptr<ExplicitBitVect>( ::RDKit::RDKFingerprintMol(*prod) );
	p.frags_ = reagents;
	return p;
}

/// @brief Setup weights for a sampler. Weighted sampling with geometric distribution; common ratio dynamically decided by population size
void
ReactionBasedAnalogSampler::geometric_sampling( numeric::random::WeightedSampler & sampler, core::Size N ) {
	sampler.clear();
	core::Real min_ratio = std::max( geo_spl_ratio_, 5.0/N ); // sample among at least top 5
	core::Real n = N * min_ratio;
	core::Real q = 1 - std::pow( min_ratio, 1/n );   // derived from geometric series summation formula, given sum of first n terms being 1-ratio
	for ( unsigned int rank = 1; rank <= N; ++rank ) {
		sampler.add_weight( q * std::pow( 1 - q, rank - 1 ) ); // An = A1 * r ^ ( n - 1 )
	}
}

core::chemical::VDVDMapping
ReactionBasedAnalogSampler::get_mapping() const {
	return mapping_;
}

std::string
ReactionBasedAnalogSampler::class_name() {
	return "ReactionBasedAnalogSampler";
}

void
ReactionBasedAnalogSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("dir", xs_string,
		"The directory containing the SMARTS-based reactions file and corresponding reagent files.", "N/A" )
		+ XMLSchemaAttribute::required_attribute("reactions", xs_string,
		"The filename of the SMARTS-based reactions file.")
		+ XMLSchemaAttribute::attribute_w_default("reagents", xs_string,
		"The filename of the SMILES-based reacgents file.", "N/A" )
		+ XMLSchemaAttribute::attribute_w_default( "sampling_ratio", xsct_real,
		"Sampling ratio for geometric sampler when picking fragments. Assumption invalid when exceeds 0.25.", "0.001" )
		+ XMLSchemaAttribute::attribute_w_default( "minCandidates", xsct_non_negative_integer,
		"The minimum number of candidates required in the sampling set before the output is drawn from each MC cycle. Increasing this will slightly increase runtime. At minimum 20 candidates.", "20" );

	AttributeList sampling_attributes;
	sampling_attributes
		+ XMLSchemaAttribute::attribute_w_default("min", xsct_real, "min sampling ratio for reset.", "0.0001")
		+ XMLSchemaAttribute::attribute_w_default("max", xsct_real, "max sampling ratio limit.", "0.001")
		+ XMLSchemaAttribute::attribute_w_default("step", xsct_real, "step sampling ratio to increase every MC cycle", "0.00002")
		+ XMLSchemaAttribute::attribute_w_default("OFF_after_n_step", xsct_real, "reset to base ratio and turn off dynamic sampling after n MC cycles", "-1")
		+ XMLSchemaAttribute::attribute_w_default("base", xsct_real, "baseline sampling ratio.", "0.0001");

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "DynamicSampling", sampling_attributes, "Set parameters to perform dynamic sampling with varying MC step size");

	protocols::chemistries::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, class_name(),
		"Reaction-based and similarity guided analoging in a given chemical library",
		attlist, subelements );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
ReactionBasedAnalogSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &) {

	geo_spl_ratio_ = tag->getOption<core::Real>("sampling_ratio", 0.001);
	minCandidates_ = tag->getOption<int>("minCandidates", 20);
	if ( minCandidates_ < 20 ) {
		TR.Warning << "minCandidates is too small. Geometrically weighted sampler may not work properly. Recommend at least 20 candidates." << std::endl;
	}
	std::string dir = tag->getOption<std::string>("dir", "N/A");
	std::string reagents = tag->getOption<std::string>("reagents", "N/A");
	if ( dir != "N/A" && reagents == "N/A" ) {
		load_reactions( dir, tag->getOption<std::string>("reactions") );
		load_all_reagents();
	} else if ( dir == "N/A" && reagents != "N/A" ) {
		load_reactions( tag->getOption<std::string>("reactions") );
		load_all_reagents( reagents );
	} else if ( dir == "N/A" && reagents == "N/A" ) {
		TR.Error << "Error XML input: Must specify either a reaction directory or a reagents filename." << std::endl;
	} else {
		TR.Error << "Error XML input: Cannot specify both a reaction directory and a reagents filename." << std::endl;
	}

	for ( utility::tag::TagCOP const & subtag : tag->getTags() ) {
		if ( subtag->getName() == "DynamicSampling" ) {
			dynamic_sampling_ = true;
			dynamic_spl_ratios_.push_back( subtag->getOption<core::Real>("min", 0.0001) );
			dynamic_spl_ratios_.push_back( subtag->getOption<core::Real>("max", 0.001) );
			dynamic_spl_ratios_.push_back( subtag->getOption<core::Real>("step", 0.00002) );
			dynamic_spl_ratios_.push_back( subtag->getOption<core::Real>("OFF_after_n_step", -1.0) );
			dynamic_spl_ratios_.push_back( subtag->getOption<core::Real>("base", 0.0001) );
		} else {
			utility_exit_with_message( "tag name " + subtag->getName() + " unrecognized." );
		}
	}

	last_smiles_ = "";
}


}
}
