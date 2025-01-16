// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/EvolutionManager.cc
/// @brief  Class definition for %EvolutionManager
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


// unit headers
#include <protocols/ligand_evolution/EvolutionManager.hh>

// package headers
#include <protocols/ligand_evolution/selectors/ElitistSelector.hh>
#include <protocols/ligand_evolution/selectors/RouletteSelector.hh>
#include <protocols/ligand_evolution/selectors/TournamentSelector.hh>

// project headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rosetta_scripts/XmlObjects.hh>

// utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/ligand_evolution.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>

// C/C++ headers
#include <fstream>
#include <iostream>

static basic::Tracer TR( "protocols.ligand_evolution.EvolutionManager" );


namespace protocols {
namespace ligand_evolution {

EvolutionManager::EvolutionManager( int rank )
:
	rank_( rank )
{}

void EvolutionManager::init() {

	// All defaults and sanity checks are supposed to happen inside EvolutionOptions!

	EvolutionOptionsOP evoopt;

	if ( !basic::options::option[ basic::options::OptionKeys::ligand_evolution::options ].active() ) {
		TR << "Path to evolutionary options file not set. Default settings will be used." << std::endl;
		evoopt = EvolutionOptionsOP(new EvolutionOptions());
	} else {
		std::string option_path = basic::options::option[basic::options::OptionKeys::ligand_evolution::options];
		TR << "Parsing evolutionary option file " << option_path << std::endl;
		evoopt = EvolutionOptionsOP(new EvolutionOptions(option_path));
	}

    // Movers and scoring functions from here will be used later. Parsing the script happens early to detect errors quickly
    TR << "Loading docking protocol from rosetta script " << evoopt->get_protocol_path() << std::endl;
    protocols::rosetta_scripts::XmlObjectsCOP rosetta_script = protocols::rosetta_scripts::XmlObjects::create_from_file( evoopt->get_protocol_path() );
    if ( !rosetta_script->list_movers().has_value( "ParsedProtocol" ) ) {
        TR.Error << "REvoLd expects to find a mover called 'ParsedProtocol' in your script. Specify it through the <PROTOCOLS> tag." << std::endl;
        TR.Error << "Available movers: " << rosetta_script->list_movers() << std::endl;
        utility_exit_with_message("No parsed protocol provided in XML script.");
    }

	max_generations_ = evoopt->get_max_generations();
    external_scoring_ = evoopt->get_external_scoring();

	library_.initialize_from_options( evoopt, external_scoring_, rank_ );

    scorer_ = ScorerOP( new Scorer( library_, evoopt->get_n_scoring_runs(), evoopt->get_ligand_chain()[ 0 ] ) );
    scorer_->initialize_from_options( evoopt, rosetta_script, rank_ );

	// These setups are not needed for external scoring runs, only evolutionary optimization
    // Results from external scoring are not intended to be loaded as pre-computed scores. External scoring uses only a smiles and an arbitrary identifier
	if ( external_scoring_ == 0 ) {
		if ( rank_ == 0 ) {

            init_evolution_protocol( evoopt );
            population_.initialize_from_evotoptions( *evoopt, library_, *scorer_ );

		} // if rank_==0
        init_workmanager();
	} // if external_scoring_==0
}

void EvolutionManager::init_workmanager() {
#ifdef USEMPI
        work_manager_ = WorkManagerOP( new WorkManager( scorer_, library_.max_positions() + 1, library_ ) );
#endif
}

void EvolutionManager::run( int mpi_size ) {

    if ( external_scoring_ ) {
        external_scoring( mpi_size );
        return;
    }
    if ( rank_ == 0 ) {

        TR << "Start evolutionary ligand optimization. Score initial population..." << std::endl;

        score();
        scorer_->save_results();
        write_population_information();
        population_.sort();

        core::Real global_best_score = population_.individual( 1 ).score();
        core::Size last_improve = 0;
        core::Real previous_best_score = global_best_score;
        calculate_quantiles();

        population_.next_generation( *selectors_[ main_selector_ ] );
        TR << print_scores() << std::endl;

        for ( core::Size generation( 1 ); generation <= max_generations_; ++generation ) {

            TR << "Producing offspring for generation " << generation << std::endl;

            // Step 1: Generate offspring
            utility::vector1< Individual > offspring;
            for ( utility::vector1< core::Size > const& offspring_option : offspring_options_ ) {

                TR.Debug << "Options: " << offspring_option << std::endl;

                // iterate through all options defining how to combine factories and selectors
                Selector& selector = *selectors_[ offspring_option[ 1 ]];
                OffspringFactory& factory = *factories_[ offspring_option[ 2 ]];
                core::Size selection_size = offspring_option[ 3 ];
                core::Size offspring_size = offspring_option[ 4 ];
                bool remove_from_pool = offspring_option[ 5 ];

                TR << "Select " << selection_size << " individuals with " << selector.name() << " for "
                    << factory.name()
                    << " to generate " << offspring_size << " new individuals. Selected ones ";
                if ( remove_from_pool ) {
                    TR << "will";
                } else {
                    TR << "won't";
                }
                TR << " be removed from pool." << std::endl;

                utility::vector1< Individual > tmp_offspring(
                    selector.apply( population_, selection_size, remove_from_pool ));
                tmp_offspring = factory.apply( tmp_offspring, offspring_size );
                offspring.insert( offspring.end(), tmp_offspring.begin(), tmp_offspring.end());

            }

            TR << "Generated " << offspring.size() << " new offspring. Start scoring..." << std::endl;

            // Step 2: Replace old individuals and score where needed
            population_.replace_population( offspring );
            score();

            population_.sort();
            core::Real new_best_score = population_.individual( 1 ).score();
            if ( new_best_score < global_best_score ) {
                last_improve = generation;
                previous_best_score = global_best_score;
                global_best_score = new_best_score;
            }

            // Step 3: Reduce population back down to its desired size with the main selector
            population_.next_generation( *selectors_[ main_selector_ ] );
            TR << print_scores() << std::endl;

            utility::vector1< core::Real > old_quantiles( quantiles );
            calculate_quantiles();

            TR << "******************************************************************************" << std::endl;
            TR << std::endl;
            TR << "Done with generation " << generation << ". Report:" << std::endl;
            TR << "\tglobal best:" << global_best_score << " (last improved in generation " << last_improve
                << " from " << previous_best_score << ")" << std::endl;
            TR << "\tbest: " << quantiles[ 1 ] << " (delta: " << quantiles[ 1 ] - old_quantiles[ 1 ] << ")"
                << std::endl;
            TR << "\t25%-quantile: " << quantiles[ 2 ] << " (delta: " << quantiles[ 2 ] - old_quantiles[ 2 ] << ")"
                << std::endl;
            TR << "\tmedian: " << quantiles[ 3 ] << " (delta: " << quantiles[ 3 ] - old_quantiles[ 3 ] << ")"
                << std::endl;
            TR << "\t75%-quantile: " << quantiles[ 4 ] << " (delta: " << quantiles[ 4 ] - old_quantiles[ 4 ] << ")"
                << std::endl;
            TR << std::endl;
            TR << "******************************************************************************" << std::endl;
            scorer_->save_results();
            write_population_information();

        } // for loop
    } // if rank == 0
#ifdef USEMPI
    else {
        work_manager_->work_loop();
    }

    work_manager_->clean_up();
#endif
}

void EvolutionManager::score() {
#ifdef USEMPI
        work_manager_->score( population_ );
#else
	scorer_->score_population( population_ );
#endif
}

std::string EvolutionManager::print_scores() const {
	std::string str = "Scores:";
	for ( Individual const& individual : population_.individuals() ) {
		str += " " + utility::to_string( individual.score() );
	}
	return str;
}

void EvolutionManager::calculate_quantiles() {
	core::Size size = population_.size();

	if ( size >= 4 ) {
		population_.sort();

		// best individual
		quantiles[1] = population_.individual(1).score();

		// 0.25 quantile
		quantiles[2] = population_.individual(size / 4).score();
		if ( size % 4 != 0 ) {
			quantiles[2] += population_.individual(size / 4 + 1).score();
			quantiles[2] /= 2;
		}

		// 0.5 quantile
		quantiles[3] = population_.individual(size / 2).score();
		if ( size % 2 != 0 ) {
			quantiles[3] += population_.individual(size / 2 + 1).score();
			quantiles[3] /= 2;
		}

		// 0.75 quantile
		quantiles[4] = population_.individual(size * 0.75).score();
		if ( (size * 3) % 4 != 0 ) {
			quantiles[4] += population_.individual(size * 0.75 + 1).score();
			quantiles[4] /= 2;
		}
	} else {
		quantiles[ 1 ] = 999.999;
		quantiles[ 2 ] = 999.999;
		quantiles[ 3 ] = 999.999;
		quantiles[ 4 ] = 999.999;
	}
}

void EvolutionManager::write_population_information() const {
	// open file
	std::ios_base::openmode ios_mode = std::ios::out | std::ios::trunc;
	std::ofstream file;
	file.open( "population.tsv", ios_mode );

	// write ids with labels (ligand id)
    file << "#\tMapping ligand identifiers to individual ids\n";
    file << "#\tAs the same ligand might be represented by multiple individuals\n";
    file << "#ligand_identifier: indivdual_id1 indivdual_id2 indivdual_id3 ...";
    file << std::endl;
	for ( auto const& ligands : scorer_->expose_id_memory() ) {
		std::string ligand_label = utility::join( ligands.first, "_" );
		file << ligand_label << ":";
		for ( core::Size id : ligands.second ) {
			file << " " << id;
		}
		file << std::endl;
	}
	file << std::endl;

	// write edges
    file << "#\tEdges in the family tree of molecules\n";
    file << "#\tChild ID->Parent ID\n";
    file << "#\t0: Part of the initial population\n";
    file << "#\tA molecule can have 1 parent if it originated from mutation or 2 in the case of crossover\n";
    file << "#Child Parent";
    file << std::endl;
	for ( std::pair< core::Size, core::Size > const& edge : population_.expose_inheritance_graph() ) {
		file << edge.first << " " << edge.second << std::endl;
	}
	file << std::endl;

	// write generation information
    file << "#\tList of individual ID for each generation\n";
    file << "#\t2 lines represent each generation, starting with the random initial population\n";
    file << "#\tThe first line is all individuals of that generation, the second is which remained after applying the main selector";
    file << std::endl;
	for ( utility::vector1< core::Size > const& generation : population_.expose_generation_log() ) {
		file << utility::join( generation, " " ) << std::endl;
	}

	// close file
	file.close();
}

    void EvolutionManager::init_evolution_protocol(EvolutionOptionsCOP options) {
        // ------------------------------------------------------------
        // Selector setup
        // ------------------------------------------------------------
        for ( std::string const& name: options->get_selector_names() ) {
            std::string const& type = options->get_selector_type(name);
            if ( type == "elitist" ) {
                selectors_.emplace_back(new ElitistSelector());
            } else if ( type == "roulette" ) {
                RouletteSelectorOP roulette(new RouletteSelector());
                roulette->consider_positive(
                        core::Size(options->get_selector_parameter(name, "consider_positive")));
                selectors_.push_back(roulette);
            } else if ( type == "tournament" ) {
                core::Size tournament_size = core::Size(
                        options->get_selector_parameter(name, "tournament_size"));
                core::Real win_accept = options->get_selector_parameter(name, "acceptance_chance");
                selectors_.emplace_back(new TournamentSelector(tournament_size, win_accept));
            } else {
                TR.Error << "Unknown selector type " << type << ". Please program setup in EvolutionManager."
                         << std::endl;
                utility_exit_with_message(
                        "EvolutionManager encountered unknown type of selector and can't set it up.");
            }
            selector_map_[name] = selectors_.size();
        }
        main_selector_ = selector_map_.at(options->get_main_selector());

        // ------------------------------------------------------------
        // Factory setup
        // ------------------------------------------------------------
        for ( std::string const& name: options->get_factory_names() ) {
            std::string const& type = options->get_factory_type(name);
            if ( type == "crossover" ) {
                factories_.emplace_back(new Crossover(library_));
            } else if ( type == "identity" ) {
                factories_.emplace_back(new IdentityFactory());
            } else if ( type == "mutator" ) {
                core::Real reaction_weight = options->get_factory_parameter(name, "reaction_weight");
                core::Real reagent_weight = options->get_factory_parameter(name, "reagent_weight");
                core::Real min_similarity = options->get_factory_parameter(name, "min_similarity");
                core::Real max_similarity = options->get_factory_parameter(name, "max_similarity");
                factories_.emplace_back(
                        new Mutator(library_, {reaction_weight, reagent_weight}, min_similarity,
                                    max_similarity));
            } else {
                TR.Error << "Unknown factory type " << type << ". Please program setup in EvolutionManager."
                         << std::endl;
                utility_exit_with_message(
                        "EvolutionManager encountered unknown type of factory and can't set it up.");
            }
            factory_map_[name] = factories_.size();
        }

        // ------------------------------------------------------------
        // Link and protocol setup
        // ------------------------------------------------------------
        // This is arguably a little confusing, but it is a leftover, and I already have to code way too much for this option system
        for ( std::pair<std::string, std::string> const& link: options->get_selector_factory_links() ) {
            std::string const& selector = link.first;
            std::string const& factory = link.second;
            offspring_options_.push_back({
                                                 selector_map_.at(selector),
                                                 factory_map_.at(factory),
                                                 core::Size(options->get_selector_parameter(selector, "size")),
                                                 core::Size(options->get_factory_parameter(factory, "size")),
                                                 core::Size(options->get_selector_parameter(selector, "remove"))
                                         });
        }
    }

    void EvolutionManager::external_scoring( int mpi_size ) {

        core::Size smiles_to_score = library_.n_unscored_smiles() / mpi_size;
        for ( core::Size ii( 1 + smiles_to_score * rank_ ); ii <= smiles_to_score * ( rank_ + 1 ); ++ii ) {
            try {
                scorer_->score_ligand( { ii, 0 }, "", external_scoring_ );
                if ( ii % 10 == 0 ) {
                    scorer_->save_external_scoring_results( core::Size( rank_ + 1 ));
                }
            } catch( ... ) {
                // print message and continue scoring
                TR.Error << "Unable to score " << library_.identifier_to_smiles( { ii, 0 } ) << std::endl;
            }
        }
        core::Size leftovers = library_.n_unscored_smiles() % mpi_size;
        // if the total numbers of smiles is not dividable by num of processors, let each one pick one smiles from the end
        if ( leftovers != 0 && core::Size( rank_ ) < leftovers ) {
            core::Size left_index = smiles_to_score * mpi_size + rank_ + 1;
            try {
                scorer_->score_ligand( { left_index, 0 }, "", external_scoring_ );
            } catch( ... ) {
                // nothing, just continue scoring
            }
        }
        scorer_->save_external_scoring_results( core::Size( rank_ + 1 ));

    }
}
}
