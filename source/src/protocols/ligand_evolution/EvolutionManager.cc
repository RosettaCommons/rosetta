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
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/ligand_docking/FinalMinimizer.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/StartFrom.hh>
#include <protocols/ligand_docking/Transform.hh>

#include <protocols/qsar/scoring_grid/ClassicGrid.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

// utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/ligand_evolution.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>
#include <utility/io/izstream.hh>

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

	// ------------------------------------------------------------
	// Option setup
	// ------------------------------------------------------------
	EvolutionOptionsOP evoopt;

	if ( !basic::options::option[ basic::options::OptionKeys::ligand_evolution::options ].active() ) {
		TR << "Path to evolutionary options file not set. Default settings will be used." << std::endl;
		evoopt = EvolutionOptionsOP(new EvolutionOptions());
	} else {
		std::string option_path = basic::options::option[basic::options::OptionKeys::ligand_evolution::options].value();
		TR << "Parsing evolutionary option file " << option_path << "..." << std::endl;
		evoopt = EvolutionOptionsOP(new EvolutionOptions(option_path));
	}

	// ------------------------------------------------------------
	// EvolutionManager setup
	// ------------------------------------------------------------
	max_generations_ = evoopt->get_max_generations();
	external_scoring_ = evoopt->get_external_scoring();

	// ------------------------------------------------------------
	// library setup
	// ------------------------------------------------------------
    // TODO switch to posestream to make it more universal -> MetaPoseInputStream
	utility::vector1<std::string> filenames = basic::options::option[basic::options::OptionKeys::in::file::s].value();
	if ( filenames.empty() ) {
		TR.Error << "No pdb provided. Please use the -in::file::s option" << std::endl;
		utility_exit_with_message("No pdb option provided.");
	}
	core::pose::PoseOP pose = core::import_pose::pose_from_file(filenames[1]);
	library_.set_pose(*pose);

	if ( external_scoring_ ) {
		library_.load_smiles(evoopt->get_path_to_external_smiles());
	} else {
		library_.load_data(evoopt->get_path_to_reactions(), evoopt->get_path_to_reagents(), rank_);
	}

	// ------------------------------------------------------------
	// Scoring function setup
	// ------------------------------------------------------------
	std::map< std::string, core::scoring::ScoreFunctionOP > scfx_memory;
	for ( const std::string& sfx_name : evoopt->get_scfx_names() ) {
		core::scoring::ScoreFunctionOP score_fct = core::scoring::get_score_function();
		score_fct->initialize_from_file( evoopt->get_scfx_wts( sfx_name ) );
		for ( const std::pair< std::string, core::Real >& reweigh : evoopt->get_scfx_reweighs( sfx_name ) ) {
			score_fct->set_weight( core::scoring::score_type_from_name( reweigh.first ), reweigh.second );
		}
		scfx_memory[sfx_name] = score_fct;
	}

	// ------------------------------------------------------------
	// Scorer setup
	// ------------------------------------------------------------
	scorer_ = ScorerOP( new Scorer( library_, evoopt->get_n_scoring_runs(), evoopt->get_ligand_chain()[ 0 ] ) );
	scorer_->set_main_term( evoopt->get_main_term() );
	scorer_->set_pose_path( evoopt->get_pose_dir_path() );
	scorer_->set_base_similarity_penalty( evoopt->get_similarity_penalty() );
	scorer_->set_similarity_penalty_threshold( evoopt->get_similarity_penalty_threshold() );
	scorer_->set_score_function( scfx_memory[ evoopt->get_main_scfx() ] );

	const std::string& score_memory_path = evoopt->get_path_score_memory();
	if ( !score_memory_path.empty() && rank_ == 0 ) {
		// this function is only called by rank 0 because it requires knowledge about the fragment library
		scorer_->load_scores( score_memory_path );
	}

	// ------------------------------------------------------------
	// Mover setup
	// ------------------------------------------------------------
    // TODO create a rosetta script xml within my option xml, parse that and retrieve the protocol mover
	for ( const std::string& mover_name : evoopt->get_mover_protocol() ) {
		const std::string& mover_type = evoopt->get_mover_type( mover_name );
		if ( mover_type == "start_from" ) {
			protocols::ligand_docking::StartFromOP start_from( new protocols::ligand_docking::StartFrom );
			start_from->chain( evoopt->get_ligand_chain() );
			core::Real x_coord = evoopt->get_mover_parameter( mover_name, "x" );
			core::Real y_coord = evoopt->get_mover_parameter( mover_name, "y" );
			core::Real z_coord = evoopt->get_mover_parameter( mover_name, "z" );
			start_from->add_coords( { x_coord, y_coord, z_coord } );
			scorer_->add_mover( start_from );
		} else if ( mover_type == "transform" ) {
			protocols::qsar::scoring_grid::GridSetOP gs( new protocols::qsar::scoring_grid::GridSet );
			gs->chain( evoopt->get_ligand_chain()[ 0 ] );   // makes no difference
			gs->width( evoopt->get_mover_parameter( mover_name, "grid_size" ) );
			protocols::qsar::scoring_grid::GridBaseOP cs( new protocols::qsar::scoring_grid::ClassicGrid );
			gs->add_grid( "classic", cs, 1.0 );
			core::Real box_size = evoopt->get_mover_parameter( mover_name, "box_size" );
			core::Real move_distance = evoopt->get_mover_parameter( mover_name, "max_move_distance" );
			core::Real angle = evoopt->get_mover_parameter( mover_name, "max_rotation_angle" );
			core::Size cycles = core::Size( evoopt->get_mover_parameter( mover_name, "cycles" ) );
			core::Real temperature = evoopt->get_mover_parameter( mover_name, "temperature" );
			protocols::moves::MoverOP transform( new protocols::ligand_docking::Transform( gs, evoopt->get_ligand_chain(), box_size, move_distance, angle, cycles, temperature ) );
			scorer_->add_mover( transform );
		} else if ( mover_type == "high_res_docker" ) {
			protocols::ligand_docking::LigandAreaOP inhibitor_dock_sc( new protocols::ligand_docking::LigandArea() );
			inhibitor_dock_sc->chain_ = evoopt->get_ligand_chain()[ 0 ];
			inhibitor_dock_sc->cutoff_ = 6.0;
			inhibitor_dock_sc->add_nbr_radius_ = true;
			inhibitor_dock_sc->all_atom_mode_ = true;
			inhibitor_dock_sc->minimize_ligand_ = 10.0;
			utility::vector1< protocols::ligand_docking::LigandAreaOP > inhibitor_dock_sc_vec{ inhibitor_dock_sc };
			protocols::ligand_docking::InterfaceBuilderOP side_chain_for_docking( new protocols::ligand_docking::InterfaceBuilder( inhibitor_dock_sc_vec ) );
			protocols::ligand_docking::MoveMapBuilderOP docking( new protocols::ligand_docking::MoveMapBuilder( side_chain_for_docking, nullptr, true ) );
			core::Size cycles = core::Size( evoopt->get_mover_parameter( mover_name, "cycles" ) );
			core::Size repack = core::Size( evoopt->get_mover_parameter( mover_name, "repack_every_nth" ) );
			protocols::ligand_docking::HighResDockerOP high_res_docker(
				new protocols::ligand_docking::HighResDocker(
				cycles, repack, scfx_memory.at( evoopt->get_mover_scfx( mover_name ) ), docking
				) );
			scorer_->add_mover( high_res_docker );
		} else if ( mover_type == "final_minimizer" ) {
			protocols::ligand_docking::LigandAreaOP inhibitor_final_sc( new protocols::ligand_docking::LigandArea() );
			inhibitor_final_sc->chain_ = evoopt->get_ligand_chain()[ 0 ];
			inhibitor_final_sc->cutoff_ = 6.0;
			inhibitor_final_sc->add_nbr_radius_ = true;
			inhibitor_final_sc->all_atom_mode_ = true;
			utility::vector1< protocols::ligand_docking::LigandAreaOP > inhibitor_final_sc_vec{ inhibitor_final_sc };
			protocols::ligand_docking::LigandAreaOP inhibitor_final_bb( new protocols::ligand_docking::LigandArea() );
			inhibitor_final_bb->chain_ = evoopt->get_ligand_chain()[ 0 ];
			inhibitor_final_bb->cutoff_ = 7.0;
			inhibitor_final_bb->add_nbr_radius_ = false;
			inhibitor_final_bb->all_atom_mode_ = true;
			inhibitor_final_bb->Calpha_restraints_ = 0.3;
			utility::vector1< protocols::ligand_docking::LigandAreaOP > inhibitor_dock_bb_vec{ inhibitor_final_bb };
			protocols::ligand_docking::InterfaceBuilderOP side_chain_for_final( new protocols::ligand_docking::InterfaceBuilder( inhibitor_final_sc_vec ) );
			protocols::ligand_docking::InterfaceBuilderOP backbone( new protocols::ligand_docking::InterfaceBuilder( inhibitor_dock_bb_vec, 3 ) );
			protocols::ligand_docking::MoveMapBuilderOP finalMM( new protocols::ligand_docking::MoveMapBuilder( side_chain_for_final, backbone, true ) );
			protocols::ligand_docking::FinalMinimizerOP final_minimizer( new protocols::ligand_docking::FinalMinimizer( scfx_memory[ evoopt->get_mover_scfx( mover_name ) ], finalMM ) );
			scorer_->add_mover( final_minimizer );
		} else {
			TR.Error << "Unknown mover type " << mover_type << ". Please program setup in EvolutionManager."
				<< std::endl;
			utility_exit_with_message(
				"EvolutionManager encountered unknown type of mover and can't set it up.");
		}
	}

	// These setups are not needed for external scoring runs, only evolutionary optimization
    // TODO make sure input format is the same as the outputted format
	if ( external_scoring_ == 0 ) {
		if ( rank_ == 0 ) {

			// ------------------------------------------------------------
			// Selector setup
			// ------------------------------------------------------------
			for ( const std::string &name: evoopt->get_selector_names() ) {
				const std::string &type = evoopt->get_selector_type(name);
				if ( type == "elitist" ) {
					selectors_.emplace_back(new ElitistSelector());
				} else if ( type == "roulette" ) {
					RouletteSelectorOP roulette(new RouletteSelector());
					roulette->consider_positive(
						core::Size(evoopt->get_selector_parameter(name, "consider_positive")));
					selectors_.push_back(roulette);
				} else if ( type == "tournament" ) {
					core::Size tournament_size = core::Size(
						evoopt->get_selector_parameter(name, "tournament_size"));
					core::Real win_accept = evoopt->get_selector_parameter(name, "acceptance_chance");
					selectors_.emplace_back(new TournamentSelector(tournament_size, win_accept));
				} else {
					TR.Error << "Unknown selector type " << type << ". Please program setup in EvolutionManager."
						<< std::endl;
					utility_exit_with_message(
						"EvolutionManager encountered unknown type of selector and can't set it up.");
				}
				selector_map_[name] = selectors_.size();
			}
			main_selector_ = selector_map_.at(evoopt->get_main_selector());

			// ------------------------------------------------------------
			// Factory setup
			// ------------------------------------------------------------
			for ( const std::string &name: evoopt->get_factory_names() ) {
				const std::string &type = evoopt->get_factory_type(name);
				if ( type == "crossover" ) {
					factories_.emplace_back(new Crossover(library_));
				} else if ( type == "identity" ) {
					factories_.emplace_back(new IdentityFactory());
				} else if ( type == "mutator" ) {
					core::Real reaction_weight = evoopt->get_factory_parameter(name, "reaction_weight");
					core::Real reagent_weight = evoopt->get_factory_parameter(name, "reagent_weight");
					core::Real min_similarity = evoopt->get_factory_parameter(name, "min_similarity");
					core::Real max_similarity = evoopt->get_factory_parameter(name, "max_similarity");
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
			for ( const std::pair<std::string, std::string> &link: evoopt->get_selector_factory_links() ) {
				const std::string &selector = link.first;
				const std::string &factory = link.second;
				offspring_options_.push_back({
					selector_map_.at(selector),
					factory_map_.at(factory),
					core::Size(evoopt->get_selector_parameter(selector, "size")),
					core::Size(evoopt->get_factory_parameter(factory, "size")),
					core::Size(evoopt->get_selector_parameter(selector, "remove"))
					});
			}

			// ------------------------------------------------------------
			// Population setup
			// ------------------------------------------------------------
            // TODO here and in general, create init functions for all options and pass them the options object - makes cleaner code and moves logic into context
			population_.set_supported_size(evoopt->get_population_supported_size());
			for ( const std::pair<const std::string, std::map<std::string, core::Size> > &init_opt: evoopt->get_pop_init_options() ) {
                // TODO change from west const to east const everywhere
				const std::string &init_type = init_opt.first;
				const std::map<std::string, core::Size> &type_options = init_opt.second;
				if ( init_type == "random" ) {
					population_.add_random(type_options.at("size"), library_);
				} else if ( init_type == "best_loaded" ) {
					utility::vector1<LigandIdentifier> best_individuals = scorer_->get_best_loaded(
						type_options.at("selection"));
					numeric::random::WeightedReservoirSampler<LigandIdentifier> sampler(type_options.at("size"));
					for ( const LigandIdentifier &indi: best_individuals ) {
						sampler.consider_sample(indi, 1.0);
					}
					TR.Debug << "Found " << best_individuals.size() << " best individuals in loaded scores. ";
					best_individuals.clear();
					sampler.samples(&best_individuals);
					if ( best_individuals.empty() ) {
						TR.Warning
							<< "No best individuals where selected from loaded scores. This is due to them not being present in available fragments."
							<< std::endl;
					} else {
						TR.Debug << " Add " << best_individuals.size()
							<< " randomly selected to initial population."
							<< std::endl;
						population_.add_individuals(best_individuals);
					}
				}
			}
		} // if rank_==0
		// ------------------------------------------------------------
		// MPI setup
		// ------------------------------------------------------------
        // TODO rename to something init workmanager
		init_mpi();
	} // if external_scoring_==0
}

void EvolutionManager::init_mpi() {
#ifdef USEMPI
        work_manager_ = WorkManagerOP( new WorkManager( scorer_, library_.max_positions() + 1, library_ ) );
#endif
}

void EvolutionManager::run( int mpi_size ) {

	if ( !external_scoring_ ) {
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
                // TODO add headers to output file to make it a bit more accessible
				write_population_information();

			} // for loop
		} // if rank == 0
#ifdef USEMPI
        else {
            work_manager_->work_loop();
        }

        work_manager_->clean_up();
#endif
	} else {    // if external scoring run
        // todo turn into a separate function which gets called at the beginning of run and returned - makes it clearer what the main protocol is
		core::Size smiles_to_score = library_.smiles_to_score() / mpi_size;
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
		core::Size leftovers = library_.smiles_to_score() % mpi_size;
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
	for ( std::pair< core::Size, core::Size > const& edge : population_.expose_inheritance_graph() ) {
		file << edge.first << " " << edge.second << std::endl;
	}
	file << std::endl;

	// write generation information
	for ( utility::vector1< core::Size > const& generation : population_.expose_generation_log() ) {
		file << utility::join( generation, " " ) << std::endl;
	}

	// close file
	file.close();
}
}
}
