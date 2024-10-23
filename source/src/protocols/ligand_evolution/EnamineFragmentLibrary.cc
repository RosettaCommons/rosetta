// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/EnamineFragmentLibrary.cc
/// @brief  The definition for %EnamineFragmentLibrary class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/EnamineFragmentLibrary.hh>

// project headers
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// RDKit headers
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

// C/C++ headers
#include <algorithm>

static basic::Tracer TR( "protocols.ligand_evolution.EnamineFragmentLibrary" ); // NOLINT(cert-err58-cpp)

namespace protocols {
namespace ligand_evolution {

// #########################################################################
// EnamineFragmentLibrary
// #########################################################################

EnamineFragmentLibrary::EnamineFragmentLibrary() = default;

EnamineFragmentLibrary::~EnamineFragmentLibrary(){
	for ( auto & fp : fingerprints_ ) {
		delete fp.second;
	}
}

void EnamineFragmentLibrary::load_data( std::string const& reaction_file_path, std::string const& reagent_file_path, core::Size rank ) {

	// if rank is not 0, we don't want to load data, but we need to know position counts

	load_reactions( reaction_file_path );

	if ( rank != 0 ) {
		TR.Debug << "Don't load and save reagent data for worker only processes." << std::endl;
		return;
	}

	load_reagents( reagent_file_path );

	core::Size total_possible_mols = 0;

	TR.Debug << "Reaction\ttotal\tpos1\tpos2\tpos.." << std::endl;
	for ( core::Size ii = 1; ii <= reactions_.size(); ++ii ) {
		// seed the random sampler with the actual number of possible molecules per reaction
		core::Size possible_mols = reactions_[ ii ]->calculate_possible_molecules();
		weighted_sampler_.add_weight( core::Real( possible_mols ) );

		// print size distribution of loaded reactions
		TR.Debug << reactions_[ ii ]->name_ << "\t" << possible_mols;
		for ( core::Size jj( 1 ); jj <= reactions_[ ii ]->n_positions(); ++jj ) {
			TR.Debug << "\t" << reactions_[ ii ]->reagents_[ jj ].size();
		}
		TR.Debug << std::endl;

		total_possible_mols += possible_mols;
	}


	if ( reagents_.empty() ) {
		TR.Error << "Failed to load reagents." << std::endl;
		utility_exit_with_message( "Failed to load reagents." );
	}
	if ( reactions_.empty() ) {
		TR.Error << "Failed to load reactions." << std::endl;
		utility_exit_with_message( "Failed to load reactions." );
	}

	TR << "Loaded " << reactions_.size() << " reactions and " << reagents_.size() << " reagents." << std::endl;
	TR << total_possible_mols << " molecules are available in total." << std::endl;

}

void EnamineFragmentLibrary::load_reactions( std::string const& reaction_file_path ) {

	// open file
	TR.Debug << "Search for reaction file '" << reaction_file_path << std::endl;
	utility::io::izstream reaction_file( reaction_file_path );
	if ( !reaction_file ) {
		TR.Error << "Can't find reaction file " << reaction_file_path << std::endl;
		utility_exit_with_message( "Unable to find reaction file." );
	}

	std::string line;
	utility::vector1< std::string > header;
	utility::vector1< std::string > required_fields{ "reaction_id", "components", "Reaction" };
	core::Size id_i( 1 ), comp_i( 1 ), reaction_i( 1 );

	while ( getline( reaction_file, line ) ) {

		utility::vector1< std::string > split_line( utility::split_whitespace( line ) );
		// skip empty lines
		if ( split_line.empty() ) { continue; }

		// set header
		if ( header.empty() ) {
			header = split_line;
			for ( std::string const& field : required_fields ) {
				if ( !header.has_value( field ) ) {
					TR.Error
						<< "Header line is either not present or does not contain the required field \""
						<< field
						<< "\""
						<< std::endl;
					utility_exit_with_message("Reaction file header is different from expectations.");
				}
			}
			id_i = header.index( required_fields.at( 1 ) );
			comp_i = header.index( required_fields.at( 2 ) );
			reaction_i = header.index( required_fields.at( 3 ) );

			continue;
		}

		// allow uncommenting
		if ( split_line[ 1 ] == "#" ) { continue; }
		// skip flawed line
		if ( split_line.size() != header.size() ) {
			TR.Warning << "Discarded flawed reaction line: " << line << std::endl;
			continue;
		}

		std::string reaction_name = split_line[ id_i ];
		std::string reaction_smiles = split_line[ reaction_i ];

		core::Size reagent_count = utility::string2Size( split_line[ comp_i ] );
		if ( reagent_count > max_reagents_ ) max_reagents_ = reagent_count;

		// construct new EnamineReaction object in place at the back of the reactions vector
		TR.Debug << "Create reaction " << reaction_name << " defined as " << reaction_smiles << " with " << reagent_count << " reagents." << std::endl;
		reactions_.emplace_back( new EnamineReaction( reaction_name, reaction_smiles, reagent_count ) );
		reaction_name_to_index_.insert( std::pair< std::string, core::Size >( reaction_name, reactions_.size() ) );

	}
}

void EnamineFragmentLibrary::load_reagents( std::string const& reagent_file_path ) {

	if ( reactions_.empty() ) {
		TR.Error << "Reagents should never be loaded before reactions. Please contact developer." << std::endl;
		utility_exit_with_message( "Programmer failed. Start becoming angry and write unfriendly messages." );
	}

	TR.Debug << "Loading reagents..." << std::endl;

	utility::io::izstream reagent_file( reagent_file_path );
	if ( !reagent_file ) {
		TR.Error << "Can't find reagent file " << reagent_file_path << std::endl;
		utility_exit_with_message( "Unable to open reagent file." );
	}

	std::string line;
	utility::vector1< std::string > header;
	utility::vector1< std::string > required_fields{ "SMILES", "synton_id", "synton#", "reaction_id" };
	core::Size smiles_i( 1 ), id_i( 1 ), pos_i( 1 ), reaction_i( 1 );

	while ( getline( reagent_file, line ) ) {

		utility::vector1< std::string > split_line( utility::split_whitespace( line ) );

		// skip empty lines
		if ( split_line.empty() ) { continue; }

		// set header
		if ( header.empty() ) {
			header = split_line;
			for ( std::string const& field : required_fields ) {
				if ( !header.has_value( field ) ) {
					TR.Error
						<< "Header line is either not present or does not contain the required field \""
						<< field
						<< "\""
						<< std::endl;
					utility_exit_with_message("Reagent file header is different from expectations.");
				}
			}
			smiles_i = header.index( required_fields.at( 1 ) );
			id_i = header.index( required_fields.at( 2 ) );
			pos_i = header.index( required_fields.at( 3 ) );
			reaction_i = header.index( required_fields.at( 4 ) );

			continue;
		}

		// allow uncommenting
		if ( split_line[ 1 ] == "#" ) { continue; }

		// skip flawed line
		if ( split_line.size() != header.size() ) {
			TR.Warning << "Discarded flawed reagent line: " << line << std::endl;
			continue;
		}

		std::string reagent_smiles = split_line[ smiles_i ];
		std::string reagent_name = split_line[ id_i ];
		core::Size position = utility::string2Size( split_line[ pos_i ] );
		core::Size reaction_index = reaction_name_to_index_.at( split_line[ reaction_i ] );

		// create reagent
		reagents_.emplace_back( new EnamineReagent( reagent_name, reagent_smiles, position, reaction_index ) );
		// add reagent to its reaction
		reactions_[ reaction_index ]->reagents_[ position ].emplace_back( reagents_.size() );

		if ( reagents_.size() % 100000 == 0 ) {
			TR.Debug << "Loaded " << reagents_.size() << " reagents." << std::endl;
		}

	}

	if ( reagents_.empty() ) {
		TR.Error << "Unable to load reagents." << std::endl;
		utility_exit_with_message( "Unable to load reagents." );
	}
}

core::conformation::ResidueOP EnamineFragmentLibrary::create_ligand( std::string const& smiles, bool create_rotamers ) const {

	// TODO Rosetta can't handle salts well, like an added Cl ion. Do I want to do something about that?


	RDKit::RWMol mutable_product( *RDKit::SmilesToMol( smiles ) );

	RDKit::MolOps::sanitizeMol( mutable_product );
	RDKit::MolOps::addHs( mutable_product );

	RDKit::DGeomHelpers::EmbedMolecule( mutable_product );

	core::chemical::rdkit::RDMolToRestype mol_to_res( mutable_product );
	core::chemical::MutableResidueTypeOP new_ligand = mol_to_res.generate_restype();

	if ( create_rotamers ) {
		core::Size n_rotamers = generate_rotamers( *new_ligand );
		TR.Debug << "Set " << n_rotamers << " rotamers." << std::endl;
	}

	core::chemical::ResidueTypeCOP restype = core::chemical::ResidueType::make( *new_ligand );

	core::conformation::ResidueOP tmp_residue( new core::conformation::Residue( restype, true ) );

	TR.Debug << "Done creating." << std::endl;

	return tmp_residue;
}

core::Size EnamineFragmentLibrary::generate_rotamers( core::chemical::MutableResidueType& new_ligand ) const {
	// This whole function is written by and copied from Rocco Moretti !

	// nconf logic adapted from Jean-Paul Ebejer's presentation
	// at the London RDKit User General Meeting
	// http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
	core::Size nconf( 50 );
	if ( new_ligand.nchi() > 12 ) {
		nconf = 300;
	} else if ( new_ligand.nchi() >= 8 ) {
		nconf = 200;
	}

	TR.Debug << "Try to create " << nconf << " rotamers." << std::endl;

	core::chemical::rotamers::StoredRotamerLibrarySpecificationOP rotamers_spec( new core::chemical::rotamers::StoredRotamerLibrarySpecification() );

	utility::vector1< std::map< std::string, core::Vector > > return_value;

	core::chemical::rdkit::RestypeToRDMol to_rdmol( new_ligand, /*neutralize = */ false, /*keep_hydro=*/ true );
	::RDKit::RWMOL_SPTR rdmol( to_rdmol.Mol() );

	core::chemical::VDIndexMapping const & rosetta_to_rdkit( to_rdmol.vd_to_index() );

	// Generate the rotamers
	::RDKit::INT_VECT conf_ids; // List of generated conformations
	::RDKit::DGeomHelpers::EmbedMultipleConfs(*rdmol, conf_ids, nconf,
		/*numThreads*/ 1, /*maxIterations*/ 30, /*seed*/ 11111, /*clearConfs*/ true, /*useRandomCoords*/ false, /*boxSizeMult*/ 2.0,   // We use a fixed seed such that conformer generation is repeatable.
		/*randNegEig*/ true, /*numZeroFail*/ 1, /*pruneRmsThresh*/ -1.0, /*coordMap*/ nullptr, /*optimizerForceTol*/ 1e-3, /*ignoreSmoothingFailures*/ false,
		/*enforceChirality*/ true, /*useExpTorsionAnglePrefs*/ true, /*useBasicKnowledge*/ true  // These are the ones we're changing
	);

	TR.Debug << "RDKit generated " << conf_ids.size() << " rotamers." << std::endl;

	// Now convert the coordinates to rotamers for Rosetta
	for ( core::Size ii(0); ii < conf_ids.size(); ++ii ) {

		::RDKit::Conformer const & conf( rdmol->getConformer(conf_ids[ii]) );
		std::map< std::string, core::Vector > single_rotamer_spec;

		for ( core::chemical::VD atm: new_ligand.all_atoms() ) {
			if ( rosetta_to_rdkit[ atm ] != rosetta_to_rdkit.invalid_entry() ) {
				::RDGeom::Point3D const & pos( conf.getAtomPos( rosetta_to_rdkit[ atm ] ) );
				//set the xyz coordinates in Rosetta
				single_rotamer_spec[ new_ligand.atom_name(atm) ] = core::Vector( pos.x, pos.y, pos.z );
			}
		}

		return_value.push_back(single_rotamer_spec );
	}

	rotamers_spec->set_rotamers(return_value );

	if ( return_value.empty() ) {
		TR.Warning << "Unable to generate rotamers." << std::endl;
		return 0;
	}

	new_ligand.rotamer_library_specification( rotamers_spec );

	return return_value.size();
}

LigandIdentifier EnamineFragmentLibrary::random_ligand() const {

	// set values for a completely random ligand
	LigandIdentifier lig( max_reagents_ + 1 );
	lig[ 1 ] = random_reaction();
	core::Size n_reagents = reactions_[ lig[ 1 ] ]->n_positions();
	for ( core::Size position = 1; position <= n_reagents; ++position ) {
		lig[ 1 + position ] = reactions_[ lig[ 1 ] ]->random_reagent_index( position );
	}

	// fill identifier with 0 values if the reaction is too short
	if ( n_reagents < max_reagents_ ) {
		for ( core::Size ii = n_reagents + 1; ii <= max_reagents_; ++ii ) {
			lig[ 1 + ii ] = 0;
		}
	}

	return lig;
}

ReagentSimilarityList EnamineFragmentLibrary::get_similar_reagents( core::Size reagent_id, core::Size reaction_id, core::Size position ) const {

	ReagentSimilarityList similar_reagents;

	// check if the reaction actually has that many positions. If it hasn't return a 0. It won't be used.
	if ( reactions_[ reaction_id ]->n_positions() < position ) {
		similar_reagents.emplace_back( 0, 0.0 );
		return similar_reagents;
	}

	// calculate all similarities between the given reagent and reagents at the given position in the given reaction
	for ( core::Size reagent : reactions_[ reaction_id ]->reagents_[ position ] ) {
		core::Real similarity( RDKit::TanimotoSimilarity( *( reagents_[ reagent_id ]->fingerprint_ ), *( reagents_[ reagent ]->fingerprint_ ) ) );
		similar_reagents.emplace_back( reagent, similarity );
	}

	// helper to sort the list
	struct {
		bool operator()( std::pair< core::Size, core::Real > const& a, std::pair< core::Size, core::Real > const& b ) const {
			return a.second > b.second;
		}
	} higherSimilarity;

	std::sort( similar_reagents.begin(), similar_reagents.end(), higherSimilarity );

	return similar_reagents;
}

core::Size EnamineFragmentLibrary::reactions_size() const {
	return reactions_.size();
}

core::Size EnamineFragmentLibrary::reagents_size() const {
	return reagents_.size();
}

core::Size EnamineFragmentLibrary::reagents_size( core::Size reaction_index ) const {
	core::Size size( 0 );
	EnamineReactionCOP reaction = reactions_[ reaction_index ];
	// iterate over all positions in the given reaction
	for ( core::Size position( 1 ); position <= reaction->reagents_.size(); ++position ) {
		size += reagents_size( reaction_index, position );
	}
	return size;
}

core::Size EnamineFragmentLibrary::reagents_size( core::Size reaction_index, core::Size position ) const {
	return reactions_[ reaction_index ]->reagents_[ position ].size();
}

core::pose::PoseOP EnamineFragmentLibrary::create_pose( core::conformation::Residue& ligand, char ligand_chain ) const {

	core::pose::PoseOP ligand_pose( new core::pose::Pose );
	ligand_pose->detached_copy( pose_ );

	core::pose::PDBInfoOP pdb_info( ligand_pose->pdb_info() );
	if ( !pdb_info ) {
		utility_exit_with_message( "Pose does not have PDBInfo - can't determine original chain" );
	}
	char original_chain = pdb_info->chain( 1 );

	ligand_pose->append_residue_by_jump( ligand, 1, "", "", true );

	core::Size ligand_position = 0;
	for ( core::Size ii( 1 ); ii <= pdb_info->nres(); ++ii ) {
		if ( pdb_info->chain( ii ) != original_chain ) {
			if ( ligand_position != 0 ) {
				TR.Warning << "Multiple residues are not on original chain " << original_chain << ". Found residue " << ligand_position << " on chain " << pdb_info->chain( ii ) << std::endl;
			}
			ligand_position = ii;
		}
	}
	if ( ligand_position == 0 ) {
		utility_exit_with_message( "Couldn't identify new ligand position." );
	}

	ligand_pose->pdb_info()->chain( ligand_position, ligand_chain );
	ligand_pose->update_pose_chains_from_pdb_chains();
	// this information doesn't persist into the pdb file, but it is useable in all other functions. So not really a problem but confusing
	TR.Debug << "Determined residue " << ligand_position << " to be the new ligand. Moved it to design chain " << ligand_chain << std::endl;

	return ligand_pose;
}

core::pose::PoseOP EnamineFragmentLibrary::create_ligand_pose( LigandIdentifier const& id, bool create_rotamers, char ligand_chain ) const {

	TR.Debug << "Turn id " << id << " into smiles." << std::endl;

	std::string smiles;
	if ( smiles_.empty() ) {
		smiles = run_reaction(id);
	} else {
		smiles = smiles_.at( id.at( 1 ) );
	}
	return create_ligand_pose( smiles, create_rotamers, ligand_chain );
}

core::pose::PoseOP EnamineFragmentLibrary::create_ligand_pose(const std::string &smiles, bool create_rotamers, char ligand_chain) const {
	TR.Debug << "Try to create ligand for " << smiles << std::endl;
	return create_pose( *create_ligand( smiles, create_rotamers ), ligand_chain );
}

void EnamineFragmentLibrary::set_pose( core::pose::Pose const& pose ) {
	pose_ = pose;
}

core::Size EnamineFragmentLibrary::random_reaction() const {
	return weighted_sampler_.random_sample();
}

core::Size EnamineFragmentLibrary::random_reaction( std::set<core::Size> const& exclude ) const {
	numeric::random::WeightedSampler tmp_sampler;
	// creates a weight for all reactions as 0
	tmp_sampler.resize( reactions_size() );
	for ( core::Size ii( 1 ); ii <= reactions_size(); ++ii ) {
		// if the reaction id is not excluded, adjust its weight. Otherwise, keep it at 0.
		if ( exclude.count( ii ) == 0 ) {
			tmp_sampler.set_weight( ii, reactions_[ ii ]->possible_molecules_ );
		}
	}
	return tmp_sampler.random_sample();
}

std::string const& EnamineFragmentLibrary::enamine_reaction( core::Size reaction_id ) const {
	return reactions_[ reaction_id ]->name_;
}

std::string const& EnamineFragmentLibrary::enamine_reagent( core::Size reagent_id ) const {
	return reagents_[ reagent_id ]->name_;
}

std::string EnamineFragmentLibrary::run_reaction( LigandIdentifier const& identifier ) const {
	core::Size reaction_index( identifier[ 1 ] );
	core::Size n_reaction_positions = reactions_[ reaction_index ]->n_positions();
	utility::vector1< core::Size > reagent_indices( identifier.begin() + 1, identifier.end() );

	// make sure that the given index actually refers to a reagent useable by the reaction
	for ( core::Size ii( 1 ); ii <= n_reaction_positions; ++ii ) {
		assert( reagent_indices[ ii ] >= 1 && reagent_indices[ ii ] <= reagents_.size() );
		if ( !reactions_[ reaction_index ]->reagents_[ ii ].has_value( reagent_indices[ ii ] ) ) {
			TR.Error << "Reaction " << reactions_[ reaction_index ]->name_ << " can't use " << reagents_[ reagent_indices[ ii ] ]->name_ << " at position " << ii << std::endl;
			utility_exit_with_message( "Invalid reagent selection" );
		}
	}

	TR.Debug << "Running reaction: " << reactions_[ reaction_index ]->name_ << " with";
	for ( core::Size ii( 1 ); ii <= n_reaction_positions; ++ii ) {
		TR.Debug << " reagent" << ii << " " << reagents_[ reagent_indices[ ii ] ]->name_;
	}
	TR.Debug << std::endl;

	// rdkit requires some weird vector and vector of vector shenanigans
	RDKit::MOL_SPTR_VECT react;
	for ( core::Size ii( 1 ); ii <= n_reaction_positions; ++ii ) {
		react.push_back( RDKit::RWMOL_SPTR( reagents_[ reagent_indices[ ii ] ]->mol_ ) );
	}

	auto products = reactions_[ reaction_index ]->reac_->runReactants( react );
	TR.Debug << "Len of outer vector: " << products.size() << std::endl;
	TR.Debug << "Len of inner vector: " << products[ 0 ].size() << std::endl;

	return RDKit::MolToSmiles( RDKit::RWMol( *products[ 0 ][ 0 ] ) );
}

std::string EnamineFragmentLibrary::identifier_to_smiles( LigandIdentifier const& identifier ) const {
	if ( smiles_.empty() ) {
		return run_reaction( identifier );
	} else {
		return smiles_.at( identifier.at( 1 ) );
	}
}

void EnamineFragmentLibrary::load_smiles( std::string const& path_to_data ) {
	utility::io::izstream smiles_file( path_to_data );
	TR << "Try to open file " << path_to_data << std::endl;
	if ( !smiles_file ) {
		TR.Error << "Can't find smiles file " << path_to_data << std::endl;
		utility_exit_with_message( "Unable to open smiles file." );
	}
	std::string line;
	while ( getline( smiles_file, line ) ) {
		smiles_.push_back( utility::string_split( line, ';' )[ 1 ] );
	}
	TR << "Loaded " << smiles_.size() << " molecules to test." << std::endl;
}

core::Size EnamineFragmentLibrary::smiles_to_score() const {
	return smiles_.size();
}

core::Size EnamineFragmentLibrary::max_positions() const {
	return max_reagents_;
}

core::Size EnamineFragmentLibrary::reaction_positions( core::Size reaction_id ) const {
	return reactions_[ reaction_id ]->n_positions();
}

core::Size EnamineFragmentLibrary::reaction_name_to_index(const std::string &reaction_name) const {
	if ( reaction_name_to_index_.count( reaction_name ) == 0 ) {
		return 0;
	}
	return reaction_name_to_index_.at( reaction_name );
}

core::Size EnamineFragmentLibrary::reagent_name_to_index(core::Size reaction_index, core::Size position, const std::string& reagent_name) const {
	for ( core::Size reagent_index : reactions_[ reaction_index ]->reagents_[ position ] ) {
		if ( reagent_name == reagents_[ reagent_index ]->name_ ) {
			return reagent_index;
		}
	}
	return 0;
}

RDKit::SparseIntVect<unsigned int>* EnamineFragmentLibrary::calculate_fingerprint(const LigandIdentifier &id) {
	if ( fingerprints_.count( id ) == 0 ) {
		std::string smiles = run_reaction(id);
		RDKit::RWMol* mol = RDKit::SmilesToMol( smiles );
		fingerprints_[ id ] = RDKit::MorganFingerprints::getFingerprint( *mol, 2 );
		delete mol;
	}
	return fingerprints_[ id ];
}

core::Real EnamineFragmentLibrary::similarity(const LigandIdentifier &id1, const LigandIdentifier &id2) {
	return RDKit::TanimotoSimilarity( *calculate_fingerprint( id1 ), *calculate_fingerprint( id2 ) );
}

// #########################################################################
// EnamineReaction
// #########################################################################

EnamineReaction::EnamineReaction( std::string const& name, std::string const& reaction_smiles, core::Size n_reagents )
:
	name_( name )
{
	reac_ = RDKit::RxnSmartsToChemicalReaction( reaction_smiles );
	reac_->initReactantMatchers();
	for ( core::Size ii( 1 ); ii <= n_reagents; ++ii ) {
		reagents_.emplace_back( );
	}
}

EnamineReaction::~EnamineReaction() {
	delete reac_;
}

core::Size EnamineReaction::random_reagent_index( core::Size position ) const {
	return reagents_[ position ][ numeric::random::random_range( 1, reagents_[ position ].size() ) ];
}

core::Size EnamineReaction::calculate_possible_molecules() {
	for ( core::Size ii = 1; ii <= reagents_.size(); ii++ ) {
		possible_molecules_ *= reagents_[ ii ].size();
	}
	if ( possible_molecules_ < 1 ) {
		possible_molecules_ = 1;
	}
	return possible_molecules_;
}

core::Size EnamineReaction::n_positions() const {
	return reagents_.size();
}

// #########################################################################
// EnamineReagent
// #########################################################################

EnamineReagent::EnamineReagent( std::string const& name, std::string const& reagent_smiles, core::Size reaction_position, core::Size reaction_index )
:
	name_( name ),
	reaction_position_( reaction_position ),
	reaction_index_( reaction_index ),
	mol_( RDKit::SmilesToMol( reagent_smiles ) )
{
	fingerprint_ = RDKit::MorganFingerprints::getFingerprint( *mol_, 2 );
}

EnamineReagent::~EnamineReagent() {
	delete fingerprint_;
}

}
}
