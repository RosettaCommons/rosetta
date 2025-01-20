// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Scorer.vv
/// @brief  Implementation of the %Scorer class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/Scorer.hh>

// utility headers
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>

// project headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/ligand_docking/ligand_scores.hh>
#include <protocols/ligand_docking/StartFrom.hh>

// C/C++ headers
#include <algorithm>
#include <fstream>
#include <iostream>

static basic::Tracer TR( "protocols.ligand_evolution.Scorer" ); // NOLINT(cert-err58-cpp)


namespace protocols {
namespace ligand_evolution {

Scorer::Scorer( FragmentLibrary& library, core::Size n_runs, char ligand_chain )
:
	ligand_chain_( ligand_chain ),
	n_runs_( n_runs ),
	library_( library ),
	best_pose_( new core::pose::Pose ),
	working_pose_( new core::pose::Pose )
{}

void Scorer::add_mover( moves::MoverOP const& mover ) {
	mover_.push_back( mover );
	TR.Debug << "Added mover " << mover->get_name() << " at position " << mover_.size() << std::endl;
}

void Scorer::score_population( Population& pop ) {

	// individual level
	for ( Individual& individual : pop.individuals() ) {
		score_individual( individual );
	}

	// population level
	core::Size total_sim_counter = 0;
	core::Real min_penalty = 99999.999;
	core::Real max_penalty = 0.00;
	core::Size sim_mol_counter = 0;
	pop.sort();
	// the best scoring molecules do not receive similarity penalties, but every following molecule gets one penalty for every similar mol with better score
	for ( core::Size ii(1); ii <= pop.size(); ++ii ) {
		core::Size sim_counter = 0;
		core::Real highest_sim = 0.00;
		for ( core::Size jj(1); jj < ii; ++jj ) {
			core::Real sim = library_.similarity( pop.individual(ii).identifier(), pop.individual(jj).identifier() );
			if ( sim > highest_sim ) highest_sim = sim;
			if ( sim > similarity_penalty_threshold_ ) {
                TR.Trace << "Individual " << pop.individual( ii ).id() << " at index " << ii << " has " << sim << " similarity to individual "
                << pop.individual( jj ).id() << " at index " << jj << std::endl;
				sim_counter++;
			}
		}
		if ( sim_counter > 0 ) {
			Individual& individual = pop.individual( ii );
			core::Real penalty = ((core::Real)sim_counter) * base_similarity_penalty_;
			if ( penalty < min_penalty ) min_penalty = penalty;
			if ( penalty > max_penalty ) max_penalty = penalty;
			core::Real new_score = individual.score() + penalty;
			individual.score( new_score );
			sim_mol_counter++;
			TR.Debug << "Individual " << individual.id() << " with ligand " << individual.identifier() << " got a similarity penalty of " << penalty << " for exceeding the similarity threshold " << sim_counter << " times. New score: " << new_score << std::endl;
		}
		total_sim_counter += sim_counter;
	}
	pop.sort();
	if ( total_sim_counter > 0 ) {
		TR << "Population recorded " << total_sim_counter << " similarities exceeding threshold " << similarity_penalty_threshold_ << ". Applied penalties between ";
		TR << min_penalty << " and " << max_penalty << " to " << sim_mol_counter << " of " << pop.size() << " molecules." << std::endl;
	}
}

void Scorer::score_individual( Individual& individual ) {
	score_ligand( individual.identifier() );
	for ( std::string const& term : score_terms_ ) {
		// add all score terms to individual
		individual.score( term, score_memory_.at( individual.identifier() ).at( term ) );
	}
	// set the main score
	individual.score( score_memory_.at( individual.identifier() ).at( main_score_term_ ) );
	if ( id_memory_.count( individual.identifier() ) == 0 ) {
		id_memory_[ individual.identifier() ] = std::set< core::Size > ();
	}
	id_memory_.at( individual.identifier() ).insert( individual.id() );
}

void Scorer::score_ligand( LigandIdentifier const& ligand, std::string const& smiles, core::Size save_n_scores ) {
	if ( check_memory( ligand ) ) {
		TR.Debug << "Found score for ligand " << ligand << ". Skip scoring..." << std::endl;
		// this should set the score for ligand, so that is_scored() returns true
	}
	if ( is_scored( ligand ) ) {
		return;
	}

	set_ligand( ligand, smiles );
	bool all_done = false;
	while ( !all_done ) {
		all_done = next_step( save_n_scores );
	}
}

void Scorer::set_ligand( LigandIdentifier const& ligand, std::string const& smiles ) {
	if ( !current_ligand_.empty() ) {
		TR.Warning << "Current ligand is not yet completely scored." << std::endl;
		return;
	}
	if ( is_scored( ligand ) ) {
		TR.Warning << "Ligand " << ligand << " is already scored and won't be set." << std::endl;
		return;
	}
	current_ligand_ = ligand;
	current_ligand_smiles_ = smiles;
	current_runs_ = 0;
}

bool Scorer::next_step( core::Size save_n_scores ) {
	if ( current_ligand_.empty() ) {
		TR.Error << "No ligand set." << std::endl;
		utility_exit_with_message( "Unable to perform a next step when no ligand was set." );
	}
	try {
		// check which steps are done, initiate next step, return
		// 1 create rotamers
		if ( basic_pose_ == nullptr ) {
			TR.Debug << "Create rotamers for " << current_ligand_ << std::endl;
			create_pose();
			return false;
		}
		// 2 run all movers
		if ( current_runs_ < n_runs_ && !score_next_ ) {
			current_runs_++;
			TR.Debug << "Apply movers for " << current_ligand_ << ". " << current_runs_ << " of " << n_runs_ << " runs." << std::endl;
			score_next_ = true;
			apply_movers();
			return false;
		}
		// 3 calculate all scores and save them if better than previous
		if ( score_next_ ) {
			TR.Debug << "Calculate scores for " << current_ligand_ << std::endl;
			score_next_ = false;
			calculate_scores( save_n_scores );
			return false;
		}
		// 4 return to 2 if current_runs_ < n_runs_
		// 5 dump pose to disk, empty current ligand, return true
		TR.Debug << "Done with ligand " << current_ligand_ << std::endl;
		// This seems weird, but if external scoring is performed (which is indicated by save_n_scores > 0) multiple poses are generated but it is unspecified
		// which pose is saved. It is the last one calculated with a score within the n best. But! This means if save_n_scores == 1, only the best post is kept
		// and can therefore be dumped to disk!
		if ( save_n_scores <= 1 ) {
			dump_pose();
		}
	} catch ( ... ) {
		TR.Error << "Encountered error during scoring of ligand " << current_ligand_ << " with smiles " << current_ligand_smiles_ << std::endl;
		TR.Error << "Break scoring and set scores to +9999.9 if not set yet." << std::endl;
		std::map< std::string, core::Real > scores;
		for ( std::string const& score_term : score_terms_ ) {
			scores[ score_term ] = 9999.9;
		}
		set_scores( current_ligand_, scores );
	}
	current_ligand_.clear();
	current_ligand_smiles_.clear();
	basic_pose_ = nullptr;
	best_pose_ = utility::pointer::make_shared< core::pose::Pose >();
	return true;
}

bool Scorer::is_scored( LigandIdentifier const& ligand ) const {
	return score_memory_.count( ligand ) == 1;
}

std::map< std::string, core::Real > const& Scorer::get_scores( LigandIdentifier const& ligand ) const {
	return score_memory_.at( ligand );
}

void Scorer::set_main_term( std::string const& score_term ) {
	if ( std::find( score_terms_.begin(), score_terms_.end(), score_term ) == score_terms_.end() ) {
		TR.Error << score_term << " is an invalid score term. Use one of these: " << score_terms_ << std::endl;
		utility_exit_with_message( "Tried to set invalid main score term." );
	}
	main_score_term_ = score_term;
}

void Scorer::set_pose_path( std::string const& path ) {
	pose_path_ = path;
}

void Scorer::set_score_function( core::scoring::ScoreFunctionOP score_function ) {
	score_function_ = score_function;
}

void Scorer::create_pose() {
	if ( basic_pose_ != nullptr ) {
		TR.Error << "Basic pose is already set. Something happened out of order." << std::endl;
		utility_exit_with_message( "Tried to override current basic pose." );
	}
	if ( current_ligand_smiles_.empty() ) {
		basic_pose_ = library_.create_ligand_pose(current_ligand_, true, ligand_chain_);
	} else {
		basic_pose_ = library_.create_ligand_pose( current_ligand_smiles_, true, ligand_chain_ );
	}
}

void Scorer::apply_movers() {
	working_pose_->detached_copy( *basic_pose_ );
	for ( moves::MoverOP const& mover : mover_ ) {
		mover->apply( *working_pose_ );
	}
}

void Scorer::calculate_scores( core::Size save_n_scores ) {

	// Calculate all scores and save them temporarily
	std::map< std::string, core::Real > scores;
	core::Size n_heavy_atoms = basic_pose_->residue( basic_pose_->size() ).nheavyatoms();
	for ( std::string const& term : score_terms_ ) {
		if ( term == "ligand_interface_delta" ) {
			std::map< std::string, core::Real > interface_delta = protocols::ligand_docking::get_interface_deltas( ligand_chain_, *working_pose_, score_function_ );
			std::string identifier = "interface_delta_";
			identifier.push_back( ligand_chain_ );
			core::Real score = interface_delta.at( identifier );
			TR.Debug << term << ":" << score << "\t";
			scores[ term ] = score;
		} else if ( term == "total_REU" ) {
			core::Real score = score_function_->score( *working_pose_ );
			TR.Debug << term << ":" << score << "\t";
			scores[ term ] = score;
		} else if ( term == "ligand_interface_delta_EFFICIENCY" ) {
			core::Real score = scores.at( "ligand_interface_delta" ) / core::Real( n_heavy_atoms );
			TR.Debug << term << ":" << score << "\t";
			scores[ term ] = score;
		} else if ( term == "lid_root2" ) {
			core::Real score = scores.at( "ligand_interface_delta" ) / std::pow( n_heavy_atoms, 1.0 / 2.0 );
			TR.Debug << term << ":" << score << "\t";
			scores[ term ] = score;
		} else if ( term == "lid_root3" ) {
			core::Real score = scores.at( "ligand_interface_delta" ) / std::pow( n_heavy_atoms, 1.0 / 3.0 );
			TR.Debug << term << ":" << score << "\t";
			scores[ term ] = score;
		} else {
			TR.Error << "You forgot to implement the calculation of " << term << std::endl;
			utility_exit_with_message( "Paul messed up and forgot something!" );
		}
	}
	TR.Debug << std::endl;

	if ( save_n_scores > 0 ) {
		// changes the current id to allow multiple scores to be saved
		if ( current_runs_ <= save_n_scores ) {
			// there are not enough scores saved so far
			current_ligand_.at( 2 ) = current_runs_;
		} else {
			// search the ligand with the worst score
			core::Size worst_index = 0;
			core::Real worst_score = -999999.9;
			for ( core::Size ii = 1; ii <= save_n_scores; ++ii ) {
				current_ligand_.at( 2 ) = ii;
				core::Real current_score = score_memory_.at( current_ligand_ ).at( main_score_term_ );
				if ( current_score > worst_score ) {
					worst_score = current_score;
					worst_index = ii;
				}
			}
			current_ligand_.at( 2 ) = worst_index;
			TR << "Worst score is " << current_ligand_ << std::endl;
		}
		TR.Warning << "Saving all scores. This should only be used for external smiles scoring" << std::endl;
	}

	// if the ligand is currently unscored, simply save all
	if ( !is_scored( current_ligand_ ) ) {
		score_memory_[ current_ligand_ ] = scores;
		best_pose_->detached_copy( *working_pose_ );
	} else {
		// check if the new pose is better than the current best
		if ( score_memory_.at( current_ligand_ ).at( main_score_term_ ) > scores.at( main_score_term_ ) ) {
			score_memory_[ current_ligand_ ] = scores;
			best_pose_->detached_copy( *working_pose_ );
		}
	}
}

void Scorer::dump_pose() {
	std::string name( pose_path_ );
	name += utility::join( current_ligand_, "_" );
	name += ".pdb";
	if ( !best_pose_->dump_pdb( name ) ) {
		TR.Error << "Unable to save pose to " << name << std::endl;
	}
}

void Scorer::set_scores( LigandIdentifier const& identifier, std::map< std::string, core::Real > const& scores ) {
	if ( is_scored( identifier ) ) {
		TR.Error << "Prevented override of scores for " << identifier << std::endl;
	} else {
		score_memory_[ identifier ] = scores;
	}
}

bool Scorer::has_ligand() const {
	return !current_ligand_.empty();
}

void Scorer::set_raw_scores( LigandIdentifier const& identifier, double const* raw_scores ) {
	std::map< std::string, core::Real > scores;
	for ( core::Size ii = 1; ii <= score_terms_.size(); ++ii ) {
		scores[ score_terms_[ ii ] ] = core::Real( raw_scores[ ii - 1 ] );
	}
	set_scores( identifier, scores );
}

core::Size Scorer::n_score_terms() const {
	return score_terms_.size();
}

utility::vector1< core::Real > Scorer::get_raw_scores( LigandIdentifier const& identifier ) const {
	utility::vector1< core::Real > raw_scores;
	for ( std::string const& term : score_terms_ ) {
		raw_scores.push_back( score_memory_.at( identifier ).at( term ) );
	}
	return raw_scores;
}

void Scorer::set_base_similarity_penalty( core::Real base_penalty ) {
	base_similarity_penalty_ = base_penalty;
}

void Scorer::save_results() const {

	utility::vector1< std::pair< core::Real, utility::vector1< std::string > > > lines;
	for ( auto const& mem : score_memory_ ) {
		lines.emplace_back( mem.second.at( main_score_term_ ), ligand_line( mem.first ) );
	}
	std::sort( lines.begin(), lines.end() );

	core::Size max_reagents = library_.max_positions();
	utility::vector1< std::string > header = { "id", "reaction" };
	for ( core::Size ii( 1 ); ii <= max_reagents; ++ii ) {
		header.emplace_back( "reagent" + utility::to_string( ii ) );
	}
	header.insert( header.end(), score_terms_.begin(), score_terms_.end() );
	header.push_back( "smiles" );

	utility::vector1< core::Size > entry_length( header.size() );
	for ( core::Size ii( 1 ); ii < header.size(); ++ii ) {
		entry_length.at( ii ) = header.at( ii ).size() + 4;
	}
	for ( std::pair< core::Real, utility::vector1< std::string > > const& line : lines ) {
		for ( core::Size ii( 1 ); ii < header.size(); ++ii ) {
			entry_length.at( ii ) = std::max( line.second.at( ii ).size() + 4, entry_length.at( ii ) );
		}
	}

	std::ios_base::openmode ios_mode = std::ios::out | std::ios::trunc;
	std::ofstream file;
	file.open( "ligands.tsv", ios_mode );
	if ( file.is_open() ) {
		for ( core::Size ii( 1 ); ii < header.size(); ++ii ) {
			file << utility::pad_right( header.at( ii ), entry_length.at( ii ), ' ' );
		}
		file << header.back() << std::endl;
		for ( std::pair< core::Real, utility::vector1< std::string > > const& line : lines ) {
			for ( core::Size ii( 1 ); ii < header.size(); ++ii ) {
				file << utility::pad_right( line.second.at( ii ), entry_length.at( ii ), ' ' );
			}
			file << line.second.back() << std::endl;
		}
	} else {
		TR.Error << "Unable to open ligands.tsv file to save results." << std::endl;
		utility_exit_with_message( "File Error" );
	}
	file.close();
}

utility::vector1< std::string > Scorer::ligand_line( LigandIdentifier const& identifier ) const {
	utility::vector1< std::string > entries;
	entries.push_back( utility::join( identifier, "_" ) );
	entries.push_back( library_.reaction_id( identifier[ 1 ] ) );
	for ( core::Size ii( 1 ); ii <= library_.reaction_positions( identifier[ 1 ] ); ++ii ) {
		entries.push_back( library_.reagent_id( identifier[ ii + 1 ] ) );
	}
	for ( core::Size ii( library_.reaction_positions( identifier[ 1 ] ) + 1 ); ii <= library_.max_positions(); ++ii ) {
		entries.push_back( "-" );
	}
	for ( std::string const& score : score_terms_ ) {
		entries.push_back( utility::to_string( score_memory_.at( identifier ).at( score ) ) );
	}
	entries.push_back( library_.identifier_to_smiles( identifier ) );
	return entries;
}

std::map< LigandIdentifier, std::set< core::Size > > const& Scorer::expose_id_memory() const {
	return id_memory_;
}

void Scorer::save_external_scoring_results( core::Size rank ) {
	TR.Warning << "External results should only be saved when external scores are calculated." << std::endl;
	// filename needs to include rank to prevent race conditions in mpi

	// one run per line: smiles;term1;term2;...;termN
	// multiple line will contain the same smiles

	std::ios_base::openmode ios_mode = std::ios::out | std::ios::app;
	std::ofstream file;
	file.open( "external_scores_" + std::to_string( rank ) + ".csv", ios_mode );
	if ( file.is_open() ) {
        std::string line;
        if ( int( file.tellp() ) == 0 ) {
            line = "smiles";
            for ( std::string const& term : score_terms_ ) {
                line += ";" +  term;
            }
            file << line << std::endl;
        }
		for ( std::pair< LigandIdentifier const, std::map< std::string, core::Real > > const& entry : score_memory_ ) {
			line = library_.identifier_to_smiles( entry.first );
			for ( std::string const& term : score_terms_ ) {
				line += ";" + std::to_string( entry.second.at( term ) );
			}
            file << line << std::endl;
		}
	} else {
		TR.Error << "Unable to open file external_scores_" + std::to_string( rank ) + ".csv to save results." << std::endl;
		utility_exit_with_message( "File Error" );
	}
	file.close();

	// Since this appends results, I need to delete all saved results after writing
	score_memory_.clear();
}

void Scorer::load_scores(std::string const& path) {

	// While I wish to use the format as the output files from single runs, I can't just save things into score_memory_ for two reasons:
	// 1) everything in score_memory_ gets printed. That is bad for run evaluations. If you intend to merge all results together, just cat the score files
	// 2) LigandIdentifiers use a numbering scheme to identify reagent and reaction. This is always working within the same run and mostly between multiple runs, but just mostly

	// plan:
	// make a new class member that holds the old scoring data
	// if a new ligand should be scored and is not yet part of score_memory_, check if it is part of the loaded memory
	// This check should be done via reagent and reaction names

	TR.Debug << path << " should be loaded." << std::endl;

	utility::io::izstream score_file( path );
	if ( !score_file ) {
		TR.Error << "Can't find score file at " << path << std::endl;
		utility_exit_with_message( "Unable to open score file." );
	}

    std::string line;
    utility::vector1< std::string > header;
    utility::vector1< std::string > required_fields{ "reaction", "reagent1", "reagent2" };

    std::map< std::string, core::Size > header_to_index;

    core::Size n_reagents ( 0 );

	core::Size counter( 0 );

	while ( getline( score_file, line ) ) {

		utility::vector1< std::string > split_line( utility::split_whitespace( line ) );

		if ( split_line.empty() ) continue;

        // allows comments
        if ( split_line[ 1 ] == "#" ) continue;

        if ( header.empty() ) {
            header = split_line;
            utility::vector1< std::string > missing_fields;
            for ( std::string const& field : required_fields ) {
                if ( !header.has_value( field ) ) missing_fields.push_back( field );
            }
            for ( std::string const& field : score_terms_ ) {
                if ( !header.has_value( field ) ) missing_fields.push_back( field );
            }
            if ( !missing_fields.empty() ) {
                TR.Error << "Header line is either not present before any data or does not contain the required field(s) " << missing_fields << std::endl;
                utility_exit_with_message("Score memory file header is different from expectations.");
            }

            for ( std::string const & field : header ) {
                if ( field.find( "reagent" ) != std::string::npos ) n_reagents++;
                header_to_index[ field ] = header.index( field );
            }

            continue;
        }

        if ( split_line.size() != header.size() ) {
            TR.Warning << "Line does not match header: " << line << std::endl;
            continue;
        }

        std::string reaction = split_line[ header_to_index[ "reaction" ] ];
        utility::vector1< std::string > scored_ligand;
        for ( core::Size ii( 1 ); ii <= n_reagents; ++ii  ) {
            std::string field = "reagent" + utility::to_string( ii );
            std::string reagent = split_line[ header_to_index[ field ] ];
            if ( reagent == "-" ) {
                scored_ligand.push_back( "0" );
            } else {
                scored_ligand.push_back( reagent );
            }
        }

        while ( scored_ligand.size() < library_.max_positions() ) {
            scored_ligand.push_back( "0" );
        }

        std::map< std::string, core::Real > scores;
        for ( std::string const& score_term : score_terms_ ) {
            scores[ score_term ] =  utility::string2Real( split_line[ header_to_index[ score_term ] ] );
        }

        if ( loaded_score_memory_.count( reaction ) == 0 ) {
            loaded_score_memory_[ reaction ] = std::map< utility::vector1< std::string >, std::map< std::string, core::Real > >();
        }
        loaded_score_memory_[ reaction ][ scored_ligand ] = scores;
        TR.Debug << "Saved scores for " << reaction << " " << scored_ligand << ". " << main_score_term_ << ": " << scores[ main_score_term_ ] << std::endl;
        ++counter;
	}

	TR << "Loaded " << counter << " scores from file. " << std::endl;

}

bool Scorer::check_memory(LigandIdentifier const& id) {

	// if scores have been loaded, check if the ligand is scored and copy scores
	if ( !loaded_score_memory_.empty() ) {
		TR.Debug << "Old scores were loaded. Create name based identifier and search for scores for " << id << std::endl;

		std::string reaction = library_.reaction_id(id[1]);
		utility::vector1< std::string > universal_identifier;
		for ( core::Size ii( 2 ); ii <= id.size(); ++ii ) {
			if ( id[ ii ] != 0 ) {
				universal_identifier.push_back( library_.reagent_id(id[ii]));
			} else {
				universal_identifier.push_back( "0" );
			}
		}
		TR.Debug << "Turned " << id << " into " << universal_identifier << " and reaction " << reaction
			<< std::endl;

		if ( loaded_score_memory_.count(reaction) == 1 ) {
            TR.Debug << "Found reaction " << reaction << ". ";
			if ( loaded_score_memory_[reaction].count(universal_identifier) == 1 ) {
				TR.Debug << "Found loaded scores for " << universal_identifier << " ." << std::endl;
				score_memory_[id] = loaded_score_memory_[reaction][universal_identifier];
				return true;
			} else {
				TR.Debug << "Unable to find synthon combination " << universal_identifier << " in loaded scores." << std::endl;
			}
		} else {
			TR.Debug << "Unable to find reaction " << reaction << " in loaded scores." << std::endl;
		}
	}
	TR.Debug << "No scores available." << std::endl;
	return false;
}

utility::vector1<LigandIdentifier> Scorer::get_best_loaded( core::Size size ) const {

	// collects all ids and scores to sort later
	utility::vector1< std::pair< core::Real, LigandIdentifier > > ligand_score_list;

	for ( std::pair< std::string const, std::map< utility::vector1< std::string >, std::map< std::string, core::Real > > > const& reaction_mem : loaded_score_memory_ ) {
		std::string const& reaction = reaction_mem.first;
		core::Size reaction_index;
		reaction_index = library_.reaction_name_to_index(reaction);
		if ( reaction_index == 0 ) {
			continue;
		}
		for ( std::pair< utility::vector1< std::string > const, std::map< std::string, core::Real > > ligand_mem : reaction_mem.second ) {
			utility::vector1< std::string > const& universal_id = ligand_mem.first;
			LigandIdentifier id{ reaction_index };
			core::Real score = ligand_mem.second.at( main_score_term_ );
			bool encountered_error = false;
			for ( core::Size pos( 1 ); pos <= universal_id.size(); ++pos ) {
				if ( universal_id[ pos ] != "0" ) {
					id.push_back(library_.reagent_name_to_index(reaction_index, pos,
						universal_id[pos]));
					if ( id.at( id.size() ) == 0 ) {
						encountered_error = true;
						break;
					}
				} else {
					id.push_back( 0 );
				}
			}
			if ( !encountered_error ) {
				ligand_score_list.emplace_back(score, id);
			}
		}
	}

	TR.Debug << "Mapped a total of " << ligand_score_list.size() << " entries to available fragments." << std::endl;

	// the first element should now be the lowest and therefore the best
	std::sort( ligand_score_list.begin(), ligand_score_list.end() );
	if ( !ligand_score_list.empty() ) {
		TR.Debug << "Best loaded score " << ligand_score_list[1].first << ", worst "
			<< ligand_score_list[ligand_score_list.size()].first << std::endl;
	}

	core::Size final_size;
	if ( size == 0 ) {
		final_size = ligand_score_list.size();
	} else {
		final_size = std::min( ligand_score_list.size(), size );
	}
	utility::vector1< LigandIdentifier > id_list;
	for ( core::Size ii( 1 ); ii <= final_size; ++ii ) {
		id_list.push_back( ligand_score_list[ ii ].second );
	}

	return id_list;
}

void Scorer::set_similarity_penalty_threshold(core::Real threshold) {
	similarity_penalty_threshold_ = threshold;
}

    void Scorer::initialize_from_options(EvolutionOptionsCOP options, protocols::rosetta_scripts::XmlObjectsCOP rosetta_script, core::Size rank) {

        set_main_term( options->get_main_term() );
        set_pose_path( options->get_pose_dir_path() );
        set_base_similarity_penalty( options->get_similarity_penalty() );
        set_similarity_penalty_threshold( options->get_similarity_penalty_threshold() );
        set_score_function( rosetta_script->get_score_function( options->get_main_scfx() ) );

        std::string const& score_memory_path = options->get_path_score_memory();
        if ( !score_memory_path.empty() && rank == 0 ) {
            // this function is only called by rank 0 because it requires knowledge about the fragment library
            load_scores( score_memory_path );
        }

        protocols::ligand_docking::StartFromOP start_from( new protocols::ligand_docking::StartFrom );
        start_from->add_coords( options->get_start_xyz() );
        start_from->chain( options->get_ligand_chain() );
        add_mover(start_from);
        add_mover( rosetta_script->get_mover("ParsedProtocol") );

    }

}
}
