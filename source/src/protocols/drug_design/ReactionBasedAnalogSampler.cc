// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionBasedAnalogSampler.hh
/// @brief apply RDKit's reaction mechanism to react reagents based on given reactions
/// @author Tracy Tang (yidan.tang@vanderbilt.edu)

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

#include <rdkit/DataStructs/BitOps.h>		// For similarity bit wrappers
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
	dynamic_sampling_( false )
{}

ReactionBasedAnalogSampler::~ReactionBasedAnalogSampler() {}

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
	while( getline( file, line ) ) {
		//load all the reactions that are available
		//set up for file is 1st column = reaction name, 2nd column = smirk string
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if( parts.size() == 0 ) { continue; }
		if( parts[1].size() < 1 || parts[1][0] == '#' ) { continue; }
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
	getline( file, line );	//first line is column names
	while( getline( file, line ) ) {
		//set up for file is 1st column = reaction id, 2nd column = No. of components, 3rd column = reaction smirks
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if( parts.size() == 0 ) { continue; }
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

void
ReactionBasedAnalogSampler::load_all_reagents() {
	for ( auto it : rxns_ ) {
		ChemicalReactionOP rxn( it.second );
		rxn->load_reagents();

		for ( core::Size rr(0); rr < rxn->nreagents(); ++rr ) {
			for ( core::Size ri(1); ri <= rxn->n_availible_reagents(rr); ++ri ) {
				struct Reagent reagent;
				// Generate fingerprint for each availible reagent
				::RDKit::ROMolOP reag = rxn->reagent( rr, ri );
				reagent.rdmol_ = reag;
//				reagent.ecfp_ = getMorganFingerprint(*reag);
//				reagent.fp_ = ::RDKit::PatternFingerprintMol(*reag);
				reagent.fp_ = ::RDKit::RDKFingerprintMol(*reag);
				reagent.rxn_ = rxn->reaction_name();
				reagent.no_ = rr;
				reagents_.push_back( reagent );
			}
		}
	}
	geometric_sampling( reagent_sampler, reagents_.size() );
	std::cout << "There are in total " << reagents_.size() << " reagents." << std::endl;
}

void
ReactionBasedAnalogSampler::load_all_reagents( std::string const & filename ) {
	utility::io::izstream file( filename );
	std::string line;
	getline( file, line );	//first line is column names
	while( getline( file, line ) ) {
		//set up for file is 1st column = SMILES, 2nd column = reactant id, 3rd column = component #, 4th column = reaction id
		utility::vector1< std::string > parts( utility::split_whitespace(line) );
		if( parts.size() == 0 ) { continue; }
		if ( parts.size() < 4 ) {
			TR.Error << "ERROR: Reactant line malformed: " << line << std::endl;
			continue;
		}
		struct Reagent reagent;
		// Generate fingerprint for each availible reagent
		::RDKit::ROMolOP reag = ::RDKit::RWMOL_SPTR( ::RDKit::SmilesToMol( parts[1] ) );
		if ( check_reagent_validity( reag ) ) {
			reagent.rdmol_ = reag;
//			reagent.ecfp_ = getMorganFingerprint(*reag);
//			reagent.fp_ = ::RDKit::PatternFingerprintMol(*reag);
			reagent.fp_ = ::RDKit::RDKFingerprintMol(*reag);
			reagent.no_ = (core::Size)std::stoi(parts[3]);
			reagent.rxn_ = parts[4];
			reagents_.push_back( reagent );
		}
	}
	geometric_sampling( reagent_sampler, reagents_.size() );
	TR << "There are in total " << reagents_.size() << " reagents." << std::endl;
}

bool
ReactionBasedAnalogSampler::check_reagent_validity( ::RDKit::ROMOL_SPTR reag ) const {
	// A collection of SMARTS that will fail RDKit later on
	static const utility::vector0< ::RDKit::ROMOL_SPTR > BAD_SMARTS {
		::RDKit::RWMolOP(::RDKit::SmartsToMol("[*](S(F)(F)(F)(F)F)"))
	};
	::RDKit::MatchVectType matchvect;

	for ( core::Size i(0); i < BAD_SMARTS.size(); i++ ) {
		if ( ::RDKit::SubstructMatch( *reag, *(BAD_SMARTS[i]), matchvect ) ) {
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
//	::RDKit::SparseIntVect<unsigned int>* ecfp_rdmol = getMorganFingerprint( *rdmol );
//	ExplicitBitVect* fp_rdmol = ::RDKit::PatternFingerprintMol(*rdmol);
	ExplicitBitVect* fp_rdmol = ::RDKit::RDKFingerprintMol(*rdmol);

	if ( rxns_.size() == 0 ) {
		TR.Warning << "No suitable reactions found. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	// container of pairs: first term is similarity score, second term is original index of the fragment (to keep track after sorting)
	std::list< std::pair< core::Real, core::Size > > score_idx_set;

	// Generate similarity score and reagent index pairs for all reagents
	for ( core::Size ri(1); ri <= reagents_.size(); ++ri ) {

//		::RDKit::SparseIntVect<unsigned int>* ecfp_reag = reagents_[ri].ecfp_;
		ExplicitBitVect* fp_reag = reagents_[ri].fp_;

		// Compare to the entire molecule and generate a similarity score as the weight
		core::Real score( TverskySimilarity( *fp_rdmol, *fp_reag, 0.1, 0.9 ) );
		//core::Real score( ::RDKit::TanimotoSimilarity( *fp_rdmol, *fp_reag ) );
//		core::Real score( TverskySimilarity( *ecfp_rdmol, *ecfp_reag, 0.1, 0.9 ) );
//		core::Real score( ( score1 + score2 ) / 2 );
		score_idx_set.push_front( std::make_pair( score, ri ) );
	}

	// Sort the pairs so that the most similar fragments to the input are on the top of the list
	score_idx_set.sort(std::greater< std::pair< core::Real, core::Size > >());
//	std::sort( score_idx_set.begin(), score_idx_set.end(), sortbySim );

	// If use dynamic sampling, update the geometric sampling ratio for current cycle
	if ( dynamic_sampling_ ) {
		// If OFF_after_n step is set a valid MC cycle, set the sampling to base ratio and turn off dynamic sampling
		if ( dynamic_spl_ratios[3] > 0 && ( dynamic_spl_ratios[3] <= visited_.size() ) ) {
			geo_spl_ratio = dynamic_spl_ratios[4];
			dynamic_sampling_ = false;
		} else {
			std::string curr_smiles = ::RDKit::MolToSmiles( *rdmol );
			if ( last_smiles_ != "" && last_smiles_ != curr_smiles ) {		// last MC cycle accepted
				reset_spl_ratio();
			} else {														// last MC cycle rejected
				geo_spl_ratio = std::min( dynamic_spl_ratios[1], geo_spl_ratio + dynamic_spl_ratios[2]);
			}
			geometric_sampling( reagent_sampler, reagents_.size() - visited_.size() );
			last_smiles_ = curr_smiles;
			// debuggg
			//		std::ofstream f;
			//		f.open("spl_ratio.debug", std::fstream::app);
			//		f << geo_spl_ratio << std::endl;
			//		f.close();
		}
	}

	// Sample from the fragments and propose a small list of candidates
	utility::vector1< Product > candidates = sample( score_idx_set );

	if ( candidates.empty() ) {	// If no candidate found from last step
		TR.Error << "ERROR: No candidate product yielded from sampling.";
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	Product product = sample_candidate( candidates, rdmol );
	TR << "Product made: " << product.smiles_ << std::endl;
//	utility::vector1< Product > analogs = analog_search( product );
//	Product ana = sample_candidate( analogs, rdmol );
//	visited_.insert( ana );
	visited_.insert( product );
	::RDKit::RWMolOP prod = product.rdmol_;
//	::RDKit::RWMolOP prod = ana.rdmol_;
//	utility::vector0< core::Size > prod_frags = ana.frags_;
	utility::vector0< core::Size > prod_frags = product.frags_;
//	TR << "Analog made: " << ana.smiles_ << std::endl;

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
		::RDKit::DGeomHelpers::EmbedMultipleConfs( *prod, 50 );
		rxn_map = find_O3A_mapping( rdmol, prod );
		restype_prod_map = combine( to_converter.vd_to_index(), rxn_map );
	} else {
		TR << "No need for re-neighboring." << std::endl;
	}

	// If still need reneighboring
	if ( restype_prod_map[ rsd_nbr ] == restype_prod_map.invalid_entry() ) {

		TR << "Cannot map neighbor atom to product. Re-assigning neighbor atom..." << std::endl; // Previous neighbor atom not in common substructure

		std::queue< core::chemical::VD > nbr_queue;
		nbr_queue.push( rsdtype.nbr_vertex() );

		core::chemical::VD found_nbr = core::chemical::INVALID_VD;
		while (!nbr_queue.empty()) {
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
	for ( core::Size rr(0); rr < prod_frags.size(); ++rr )
	{
		core::Size i = prod_frags[ rr ];
		::RDKit::ROMolOP reag = reagents_[ i ].rdmol_;
		rsdtype.add_string_property( "fragment" + std::to_string( rr ), ::RDKit::MolToSmiles( *reag ) );
		if ( rsdtype.properties().string_properties().count( "fragment" + std::to_string( rr ) ) > 0 ) {
			TR << "Fragment property: " + ::RDKit::MolToSmiles( *reag ) + " added." << std::endl;
		} else {
			TR << "Failed to add fragment property." << std::endl;
		}
	}

	mapping_ = combine( mapping_, combine( core::chemical::VDStringMapping(*new_resop), core::chemical::StringVDMapping(rsdtype)) );

	// That's it, we're successful
	set_last_status( core::chemical::modifications::SUCCESS );
	return;
}

::RDKit::SparseIntVect<unsigned int>*
ReactionBasedAnalogSampler::getMorganFingerprint( ::RDKit::ROMol const & mol ) const {
	// Fingerprint properties
	int radius = 2;
	bool useFeatures = false;	// false for ECFP / true for FCFP

	if ( useFeatures ) {
		std::vector<unsigned int> invars( mol.getNumAtoms() );
		::RDKit::MorganFingerprints::getFeatureInvariants( mol, invars );
		return ::RDKit::MorganFingerprints::getFingerprint( mol, radius, &invars );

	} else {
		return ::RDKit::MorganFingerprints::getFingerprint( mol, radius );
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

		}else {
			TR.Error << "New product has invalid number of conformers.";
		}
	}else {
		TR.Error << "Previous product has invalid number of conformers.";
	}

	return retval;
}

utility::vector1< ReactionBasedAnalogSampler::Product >
ReactionBasedAnalogSampler::sample( std::list< std::pair< core::Real, core::Size > > & score_idx_set ) const {
	// Setup the sampling set where each reaction id has its own bucket of reagent indices
	std::unordered_map< std::string, utility::vector1< utility::vector0< core::Size > > > sample_sets;
	utility::vector1< Product > candidate_set; 	// Setup the candidate set

	int revisited = 0;
	int nLoop = 0;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	while ( candidate_set.size() < minCandidates_ ) {
		nLoop++;
		core::Size choice = reagent_sampler.random_sample() - 1;
		std::list< std::pair< core::Real, core::Size > >::iterator itr = score_idx_set.begin();
		std::advance( itr, choice );
		core::Size reag_idx = ( *itr ).second;
		Reagent sample = reagents_[ reag_idx ];
		std::string sample_rxn = sample.rxn_;
		core::Size sample_no = sample.no_;
		// Initialize the reaction bucket if first time sampling the reaction
		if ( sample_sets.find( sample_rxn ) == sample_sets.end() ) {
			core::Size rxn_nComponent = rxns_.at( sample_rxn )->nreagents();
			utility::vector1< utility::vector0< core::Size > > rxn_set( rxn_nComponent, utility::vector0< core::Size >() );
			sample_sets[ sample_rxn ] = rxn_set;
		}
		sample_sets[ sample_rxn ][ sample_no ].push_back( reag_idx );
		// Generate all products for current sample's reaction and add to candidate set
		pair( sample_no, sample_sets[ sample_rxn ], candidate_set, revisited );
		// Remove the element so that sampler don't re-pick the same on the next cycle
		score_idx_set.erase( itr );
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	TR << "Sampling done in " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0  << " sec and looped " << nLoop << std::endl;

	// debuggg
//	std::ofstream sample_f;
//	sample_f.open("sample.debug", std::fstream::app);
//	sample_f << "Revisited " << revisited << " ligands." << std::endl;
//	sample_f.close();
	return candidate_set;
}

utility::vector1< ReactionBasedAnalogSampler::Product >
ReactionBasedAnalogSampler::sample_fragment( utility::vector0< std::list< std::pair< core::Real, core::Size > > > & score_idx_set ) const {
	// Setup the sampling set where each reaction id has its own bucket of reagent indices
	std::unordered_map< std::string, utility::vector1< utility::vector0< core::Size > > > sample_sets;
	utility::vector1< Product > candidate_set;			// Setup the candidate set
	std::unordered_set< core::Size > sampled_reagents;	// Avoid resampling the same reagent

	int revisited = 0;
	int nLoop = 0;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	while ( candidate_set.size() < minCandidates_ ) {
		nLoop++;
		for ( core::Size fi(0); fi < score_idx_set.size(); ++fi ) {
			core::Size choice = reagent_sampler.random_sample() - 1;
			std::list< std::pair< core::Real, core::Size > >::iterator itr = score_idx_set[ fi ].begin();
			std::advance( itr, choice );
			core::Size reag_idx = ( *itr ).second;
			// Remove the element so that sampler don't re-pick the same on the next cycle
			score_idx_set[ fi ].erase( itr );
			if ( sampled_reagents.find( reag_idx ) == sampled_reagents.end() ) {
				Reagent sample = reagents_[ reag_idx ];
				std::string sample_rxn = sample.rxn_;
				core::Size sample_no = sample.no_;
				// Initialize the reaction bucket if first time sampling the reaction
				if ( sample_sets.find( sample_rxn ) == sample_sets.end() ) {
					core::Size rxn_nComponent = rxns_.at( sample_rxn )->nreagents();
					utility::vector1< utility::vector0< core::Size > > rxn_set( rxn_nComponent, utility::vector0< core::Size >() );
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
	utility::vector0< core::Size > frags = prod.frags_;
	utility::vector0< ExplicitBitVect* > frags_fp;
	utility::vector0< core::Size > ana_frags;
	utility::vector0< std::list< std::pair< core::Real, core::Size > > > score_idx_analog_set;
	// Initialization
	for ( core::Size fi(0); fi < frags.size(); ++fi ) {
		std::list< std::pair< core::Real, core::Size > > score_idx_frag_set;
		score_idx_analog_set.push_back( score_idx_frag_set );
		core::Size frag_i = frags[ fi ];
		frags_fp.push_back( reagents_[ frag_i ].fp_ );
	}

	for ( core::Size ri(1); ri <= reagents_.size(); ++ri ) {
		ReactionBasedAnalogSampler::Reagent reag = reagents_[ ri ];
		// pattern fingerprint seems better suited for analoging than ECFP since it is substructure-based
		ExplicitBitVect* fp_reag = reag.fp_;
		for ( core::Size fi(0); fi < frags.size(); ++fi ) {
			ExplicitBitVect* fp_sample = frags_fp[ fi ];
			// Compare to the chosen reagent and generate a similarity score
			core::Real score( TanimotoSimilarity( *fp_sample, *fp_reag ) );
			score_idx_analog_set[ fi ].push_front( std::make_pair( score, ri ) );
		}
	}

	for ( core::Size fi(0); fi < frags.size(); ++fi ) {
		// Sort the score_index set before sampling
		score_idx_analog_set[ fi ].sort(std::greater< std::pair< core::Real, core::Size > >());
	}

	// Sample from the fragments and propose a small list of candidates
	utility::vector1< Product > candidates = sample_fragment( score_idx_analog_set );
	return candidates;
}

ReactionBasedAnalogSampler::Product
ReactionBasedAnalogSampler::sample_candidate( const utility::vector1< Product > & candidate_set, const ::RDKit::RWMOL_SPTR rdmol ) {
	// container of pairs for candidate sampling: first term is similarity score, second term is original index of the candidate (to keep track after sorting)
	utility::vector1< std::pair< core::Real, core::Size > > score_idx_candidate_set;
	numeric::random::WeightedSampler product_sampler;
//	::RDKit::SparseIntVect<unsigned int>* ecfp_rdmol = getMorganFingerprint( *rdmol );
//	ExplicitBitVect* fp_rdmol = ::RDKit::PatternFingerprintMol(*rdmol);
	ExplicitBitVect* fp_rdmol = ::RDKit::RDKFingerprintMol(*rdmol);

	// Compare each candidate to the input query molecule and rank the set based on similarity

//	core::Real totalSim = 0;
//	core::Real minSim = 1;
//	core::Real maxSim = 0;
	for ( core::Size i(1); i <= candidate_set.size(); ++i ) {
		core::Real sim = TanimotoSimilarity( *fp_rdmol, *(candidate_set[i].fp_) );
//		core::Real sim = TanimotoSimilarity( *ecfp_rdmol, *(candidate_set[i].ecfp_) );
		score_idx_candidate_set.push_back( std::make_pair( sim, i ) );
//		totalSim += sim;
//		if ( sim < minSim ) { minSim = sim; }
//		if ( sim > maxSim ) { maxSim = sim;	}
	}

	std::sort( score_idx_candidate_set.begin(), score_idx_candidate_set.end(), sortbySim );

	// debuggg
//	std::ofstream sim_f;
//	sim_f.open("sample.sim", std::fstream::app);
//	sim_f << minSim << ',' << maxSim << ',' << totalSim / ( candidates.size() ) << std::endl;
//	sim_f.close();

	// Pick a product and update the visited collector
	geometric_sampling( product_sampler, candidate_set.size() );
	core::Size choice = product_sampler.random_sample();
	core::Size cand_idx = score_idx_candidate_set[ choice ].second;
	return candidate_set[ cand_idx ];
}

bool
ReactionBasedAnalogSampler::sortbySim( const std::pair< core::Real, core::Size >& a, const std::pair< core::Real, core::Size >& b ) {
	return (a.first > b.first);
}

void
ReactionBasedAnalogSampler::pair(
		core::Size r_no,
		utility::vector1< utility::vector0< core::Size > > const & sets,
		utility::vector1< ReactionBasedAnalogSampler::Product > & candidates,
		int & revisits
) const {
	utility::vector0< utility::vector0< core::Size > > pairs;

	// Make sure each reactant bucket is not empty. Return if any bucket is empty
	for ( core::Size i(1); i <= sets.size(); i++ ) {
		if ( sets[i].empty() ) {
			return;
		}
	}
	// Run a DFS to generate all possible pairs, starting from the first reactant
	pair( 1, r_no, utility::vector0< core::Size >(), pairs, sets );

	for ( core::Size i(0); i < pairs.size(); i++ ) {
		// Convert each pair to product
		utility::vector0< core::Size > p = pairs[i];
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

void
ReactionBasedAnalogSampler::pair(
		core::Size curr_no, 											// current layer (reactant No.) of DFS
		core::Size r_no, 												// reactant No. to locate the current sample
		utility::vector0< core::Size > single_pair, 					// current pair to work on
		utility::vector0< utility::vector0< core::Size > > & all_pairs, // collection of all pairs
		utility::vector1< utility::vector0< core::Size > > const & sets // the current reaction bucket
) const {
	if ( curr_no > sets.size() ) {
		all_pairs.push_back( single_pair );
		return;
	}
	utility::vector0< core::Size > reagent_set = sets[ curr_no ];
	if ( curr_no == r_no ) {
		core::Size reag_idx = reagent_set.back();
		single_pair.push_back( reag_idx );
		pair( curr_no + 1, r_no, single_pair, all_pairs, sets );
	} else {
		utility::vector0< core::Size > temp_pair = single_pair;
		for ( core::Size i(0); i < reagent_set.size(); i++ ) {
			core::Size reag_idx = reagent_set[ i ];
			single_pair.push_back( reag_idx );
			pair( curr_no + 1, r_no, single_pair, all_pairs, sets );
			single_pair = temp_pair;
		}
	}
}

ReactionBasedAnalogSampler::Product
ReactionBasedAnalogSampler::run_reaction( utility::vector0< core::Size > const & reagents ) const {
	std::string rxn_name = reagents_[ reagents[0] ].rxn_;
	utility::vector0< Reagent > reagent_selection;
	for ( core::Size i : reagents ) {
		reagent_selection.push_back( reagents_[i] );
	}
	std::sort( reagent_selection.begin(), reagent_selection.end() );
	::RDKit::MOL_SPTR_VECT rdreagents;
	for ( core::Size ii(0); ii < reagent_selection.size(); ++ii ) {
		::RDKit::ROMolOP reag = reagent_selection[ii].rdmol_;
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
//	p.ecfp_ = getMorganFingerprint(*prod);
//	p.fp_ = ::RDKit::PatternFingerprintMol(*prod);
	p.fp_ = ::RDKit::RDKFingerprintMol(*prod);
	p.frags_ = reagents;
	return p;
}

void
ReactionBasedAnalogSampler::geometric_sampling( numeric::random::WeightedSampler & sampler, core::Size N ) {
	sampler.clear();
	core::Real min_ratio = std::max( geo_spl_ratio, 5.0/N );	// sample among at least top 5
	core::Real n = N * min_ratio;
	core::Real q = 1 - std::pow( min_ratio, 1/n );			// derived from geometric series summation formula, given sum of first n terms being 1-ratio
	for ( unsigned int rank = 1; rank <= N; ++rank ) {
		sampler.add_weight( q * std::pow( 1 - q, rank - 1 ) );	// An = A1 * r ^ ( n - 1 )
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
			"Start with a given list of reactions, select the proper reagents to react.",
			attlist, subelements );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
ReactionBasedAnalogSampler::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &) {

	geo_spl_ratio = tag->getOption<core::Real>("sampling_ratio", 0.001);
	minCandidates_ = tag->getOption<int>("minCandidates", 20);
	if ( minCandidates_ < 20 ) {
		TR << "minCandidates is too small. Using 20 minCandidates instead..." << std::endl;
		minCandidates_ = 20;
	}
	std::string dir = tag->getOption<std::string>("dir", "N/A");
	std::string reagents = tag->getOption<std::string>("reagents", "N/A");
	if ( dir != "N/A" && reagents == "N/A" ) {
		load_reactions( dir, tag->getOption<std::string>("reactions") );
		load_all_reagents();
	} else if ( dir == "N/A" && reagents != "N/A" ) {
		load_reactions( tag->getOption<std::string>("reactions") );
		load_all_reagents( reagents );
	} else {
		TR.Error << "Error XML input. Either specify reaction directory or specify reagents filename.";
	}

	for ( utility::tag::TagCOP const & subtag : tag->getTags() ) {
		if ( subtag->getName() == "DynamicSampling" ) {
			dynamic_sampling_ = true;
			dynamic_spl_ratios.push_back( subtag->getOption<core::Real>("min", 0.0001) );
			dynamic_spl_ratios.push_back( subtag->getOption<core::Real>("max", 0.001) );
			dynamic_spl_ratios.push_back( subtag->getOption<core::Real>("step", 0.00002) );
			dynamic_spl_ratios.push_back( subtag->getOption<core::Real>("OFF_after_n_step", -1.0) );
			dynamic_spl_ratios.push_back( subtag->getOption<core::Real>("base", 0.0001) );
		} else {
			utility_exit_with_message( "tag name " + subtag->getName() + " unrecognized." );
		}
	}

	last_smiles_ = "";
}


}
}
