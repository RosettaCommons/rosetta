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
#include <protocols/ligand_evolution/EvolutionOptions.hh>

// package headers

// project headers
#include <protocols/qsar/scoring_grid/ClassicGrid.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

// utility headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C/C++ headers
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/detail/file_parser_error.hpp>
#include <iostream>

static basic::Tracer TR( "protocols.ligand_evolution.EvolutionOptions" );


namespace protocols {
namespace ligand_evolution {

EvolutionOptions::EvolutionOptions() {
	check_validity();
}

EvolutionOptions::EvolutionOptions( const std::string& path_to_option_file) {

	utility::io::izstream option_file( path_to_option_file );
	if ( !option_file ) {
		TR.Error << "Can't find evolutionary option file at " << path_to_option_file << std::endl;
		utility_exit_with_message( "Unable to open option file." );
	}

	parse_option_file( path_to_option_file );

	check_validity();

	TR.Debug << "Options passed all checks." << std::endl;
}

void EvolutionOptions::check_validity() {

	core::Size error_counter = 0;

	if ( initialized_ ) {
		TR.Error << "Option system has been checked before and declared valid. This function should not be called again." << std::endl;
		error_counter++;
	}

	if ( external_scoring_ > 0 ) {

		utility::io::izstream ligand_smiles_file( path_to_external_smiles_ );
		if ( !ligand_smiles_file ) {
			TR.Error << "Can't find ligand smiles list to score at " << path_to_external_smiles_ << "." << std::endl;
			error_counter++;
		}

	} else {

		if ( generations_ <= 0 ) {
			TR.Error << "The algorithm is required to run for at least one generation." << std::endl;
			error_counter++;
		}

		utility::io::izstream reaction_file( path_to_reactions_ );
		if ( !reaction_file ) {
			TR.Error << "Can't find reactions file at " << path_to_reactions_ << "." << std::endl;
			error_counter++;
		}

		utility::io::izstream reagents_file( path_to_reagents_ );
		if ( !reagents_file ) {
			TR.Error << "Can't find reagents file at " << path_to_reagents_ << "." << std::endl;
			error_counter++;
		}

		if ( main_selector_.empty() ) {
			TR.Error << "No main selector provided." << std::endl;
			error_counter++;
		}

		if ( supported_population_size_ <= 0 ) {
			TR.Error << "Supported population size has to be greater than zero." << std::endl;
			error_counter++;
		}

		error_counter += check_selectors();
		error_counter += check_factories();
		error_counter += check_pop_init();

	}

	error_counter += check_scoring_functions();
	error_counter += check_scorer_setup();
	error_counter += check_movers();

	if ( error_counter != 0 ) {
		utility_exit_with_message( "Found a total of " + std::to_string( error_counter ) + " errors. See details above." );
	}

	initialized_ = true;

}

core::Size EvolutionOptions::check_scoring_functions() const {

	core::Size error_counter = 0;

	for ( const auto& sfx_option : score_function_wts_ ) {

		std::string sfx_name = sfx_option.first;
		std::string wts_file = sfx_option.second;

		try {
			core::scoring::find_weights_file( wts_file );
		} catch( ... ) {
			TR.Error << "Scorefunction weight file " << wts_file << " for scorefunction " << sfx_name << " cannot be found." << std::endl;
			error_counter++;
		}

	}

	for ( const auto& reweighs : score_function_reweighs_ ) {
		std::string sfx_name = reweighs.first;
		if ( score_function_wts_.count( sfx_name ) == 0 ) {
			TR.Warning << "Reweighs are defined for scorefunction named " << sfx_name  << ", but no such score function exists." << std::endl;
		}
		for ( const auto& term_reweighs : reweighs.second ) {
			std::string term_name = term_reweighs.first;
			try {
				core::scoring::score_type_from_name( term_name );
			} catch ( ... ) {
				TR.Error << "Score type " << term_name << " reweighed for scorefunction " << sfx_name << " does not exist." << std::endl;
				error_counter++;
			}
		}
	}

	return error_counter;
}

core::Size EvolutionOptions::check_scorer_setup() const {

	core::Size error_counter = 0;

	if ( !path_to_score_memory_.empty() ) {
		utility::io::izstream score_memory_file( path_to_score_memory_ );
		if ( !score_memory_file ) {
			TR.Error << "A score memory file was specified at " << path_to_score_memory_ << " but can not be found." << std::endl;
			error_counter++;
		}
	}

	if ( score_runs_ < 1 ) {
		TR.Error << "Score runs are set to less than one." << std::endl;
		error_counter++;
	}

	if ( ligand_chain_.empty() ) {
		TR.Error << "No ligand chain specified." << std::endl;
		error_counter++;
	}

	if ( main_score_term_.empty() ) {
		TR.Error << "No score term to optimize defined." << std::endl;
		error_counter++;
	} else if ( std::find( score_terms_.begin(), score_terms_.end(), main_score_term_ ) == score_terms_.end() ) {
		TR.Error << main_score_term_ << " is not available. Use one of " << score_terms_ << std::endl;
		error_counter++;
	}

	if ( pose_dump_directory_.empty() ) {
		TR.Error << "No directory to save poses has been defined." << std::endl;
		error_counter++;
	} else {
		utility::io::ozstream dummy_pdb( pose_dump_directory_ + "/dummy.pdb" );
		if ( !dummy_pdb ) {
			TR.Error << "Can't create a dummy pdb file at " << pose_dump_directory_ << "." << std::endl;
			error_counter++;
		}
	}

	if ( scoring_function_.empty() ) {
		TR.Error << "No main scoring function for scorer defined to calculate base pose energy." << std::endl;
		error_counter++;
	} else if ( score_function_wts_.count( scoring_function_ ) == 0 ) {
		TR.Error << scoring_function_ << " is not defined." << std::endl;
		error_counter++;
	}

	return error_counter;
}

core::Size EvolutionOptions::check_selectors() const {

	core::Size error_counter = 0;

	std::map< std::string, bool > is_used;

	const utility::vector1< std::string > selector_types{ "elitist", "roulette", "tournament" };

	for ( auto& selec_ops : selector_options_ ) {
		std::string selector_name = selec_ops.first;
		std::string selector_type = selec_ops.second.first;
		is_used[ selector_name ] = false;

		if ( std::find( selector_types.begin(), selector_types.end(), selector_type ) == selector_types.end() ) {
			TR.Error << selector_type << " is an unknown type of selector. Available are " << selector_types << std::endl;
			error_counter++;
		}

		if ( selec_ops.second.second.count( "size" ) == 0 ) {
			TR.Error << "Size option has to be set for all Selectors but is absent for " << selector_name << "." << std::endl;
			error_counter++;
		} else if ( selec_ops.second.second.at( "size" ) <= 0 ) {
			TR.Error << "Size option is less or equal 0 for " << selector_name << "." << std::endl;
			error_counter++;
		}

		if ( selector_type == "tournament" ) {

			if ( selec_ops.second.second.count( "tournament_size" ) == 0 ) {
				TR.Error << "tournament_size option is not set for " << selector_name << "." << std::endl;
				error_counter++;
			} else if ( selec_ops.second.second.at( "tournament_size" ) < 0 ) {
				TR.Error << "tournament_size option is less than 1 for " << selector_name << "." << std::endl;
				error_counter++;
			}

			if ( selec_ops.second.second.count( "acceptance_chance" ) == 0 ) {
				TR.Error << "acceptance_chance option is not set for " << selector_name << "." << std::endl;
				error_counter++;
			} else if ( selec_ops.second.second.at( "acceptance_chance" ) <= 0 ) {
				TR.Error << "acceptance_chance option is less or equal to 0 for " << selector_name << "." << std::endl;
				error_counter++;
			}
		}

		if ( selector_type == "roulette" ) {
			if ( selec_ops.second.second.count( "consider_positive" ) == 0 ) {
				TR.Error << "consider_positive option is not set for " << selector_name << "." << std::endl;
				error_counter++;
			}
		}

	}

	core::Size leftover_popsize = supported_population_size_;
	for ( const auto& link : selector_factory_links_ ) {
		std::string selector_name = link.first;
		std::string factory_name = link.second;
		if ( is_used.count( selector_name ) == 0 ) {
			TR.Error << "Selector " << selector_name << " linked to " << factory_name << " is not defined." << std::endl;
			error_counter++;
		} else {
			is_used[ selector_name ] = true;
			core::Size selection_size = core::Size( selector_options_.at( selector_name ).second.at( "size" ) );
			if ( selection_size > leftover_popsize ) {
				TR.Error << "Selector " << selector_name << " tries to select " << selection_size << ", but at this point only " << leftover_popsize << " individuals are available." << std::endl;
				error_counter++;
			}
			if ( bool( selector_options_.at( selector_name ).second.at( "remove" ) ) ) {
				leftover_popsize -= selection_size;
			}
		}
	}

	if ( is_used.count( main_selector_ ) == 0 ) {
		TR.Error << "Main Selector " << main_selector_ << " is not defined." << std::endl;
		error_counter++;
	} else {
		is_used[ main_selector_ ] = true;
	}

	for ( const auto& usage : is_used ) {
		if ( !usage.second ) {
			TR.Warning << "Selector " << usage.first << " is defined but never used." << std::endl;
		}
	}

	// TODO This does not check if options provided are unsupported or wrong. Maybe I should add this.

	return error_counter;
}

core::Size EvolutionOptions::check_factories() const {

	core::Size error_counter = 0;

	std::map< std::string, bool > is_used;

	const utility::vector1< std::string > factory_types{ "mutator", "crossover", "identity" };

	for ( const auto& fac_ops : factory_options_ ) {

		std::string factory_name = fac_ops.first;
		std::string factory_type = fac_ops.second.first;

		is_used[ factory_name ] = false;

		if ( std::find( factory_types.begin(), factory_types.end(), factory_type ) == factory_types.end() ) {
			TR.Error << factory_type << " is an unknown type of factory. Available are " << factory_types << std::endl;
			error_counter++;
		}

		if ( fac_ops.second.second.count( "size" ) == 0 ) {
			TR.Error << "size option is not set for factory " << factory_name << "." << std::endl;
			error_counter++;
		} else if ( fac_ops.second.second.at( "size" ) < 1.0 ) {
			TR.Error << "size option is set to less than one for factory " << factory_name << "." << std::endl;
			error_counter++;
		}

		if ( factory_type == "mutator" ) {

			bool reaction_weight_set = true;
			if ( fac_ops.second.second.count( "reaction_weight" ) == 0 ) {
				TR.Error << "reaction_weight option is not set for factory " << factory_name << "." << std::endl;
				reaction_weight_set = false;
				error_counter++;
			} else if ( fac_ops.second.second.at( "reaction_weight" ) < 0.0 ) {
				TR.Warning << "reaction_weight option is set to less than zero for factory " << factory_name << "." << std::endl;
			}

			bool reagent_weight_set = true;
			if ( fac_ops.second.second.count( "reagent_weight" ) == 0 ) {
				TR.Error << "reagent_weight option is not set for factory " << factory_name << "." << std::endl;
				reagent_weight_set = false;
				error_counter++;
			} else if ( fac_ops.second.second.at( "reagent_weight" ) < 0.0 ) {
				TR.Warning << "reagent_weight option is set to less than zero for factory " << factory_name << ". This is treated as being set to 0." << std::endl;
			}

			if ( reaction_weight_set && reagent_weight_set && fac_ops.second.second.at( "reagent_weight" ) <= 0.0 && fac_ops.second.second.at( "reaction_weight" ) <= 0 ) {
				TR.Error << "All weights for " << factory_name << " are set to less or equal than zero, causing it to mutate nothing." << std::endl;
				error_counter++;
			}

			bool min_sim_set = true;
			if ( fac_ops.second.second.count( "min_similarity" ) == 0 ) {
				TR.Error << "min_similarity option is not set for factory " << factory_name << "." << std::endl;
				min_sim_set = false;
				error_counter++;
			} else if ( fac_ops.second.second.at( "min_similarity" ) < 0.0 ) {
				TR.Error << "min_similarity option is set to less than zero for factory " << factory_name << "." << std::endl;
				error_counter++;
			} else if ( fac_ops.second.second.at( "min_similarity" ) > 1.0 ) {
				TR.Error << "min_similarity option is set to greater than one for factory " << factory_name << "." << std::endl;
				error_counter++;
			}

			bool max_sim_set = true;
			if ( fac_ops.second.second.count( "max_similarity" ) == 0 ) {
				TR.Error << "max_similarity option is not set for factory " << factory_name << "." << std::endl;
				max_sim_set = false;
				error_counter++;
			} else if ( fac_ops.second.second.at( "max_similarity" ) < 0.0 ) {
				TR.Error << "max_similarity option is set to less than zero for factory " << factory_name << "." << std::endl;
				error_counter++;
			} else if ( fac_ops.second.second.at( "max_similarity" ) > 1.0 ) {
				TR.Error << "max_similarity option is set to greater than one for factory " << factory_name << "." << std::endl;
				error_counter++;
			}

			if ( min_sim_set && max_sim_set && fac_ops.second.second.at( "max_similarity" ) <= fac_ops.second.second.at( "min_similarity" ) ) {
				TR.Error << "min_similarity is set to greater or equal than max_similarity for factory " << factory_name << "." << std::endl;
				error_counter++;
			}
		}

	}

	core::Size generated_popsize = 0;
	for ( const auto& link : selector_factory_links_ ) {
		std::string selector_name = link.first;
		std::string factory_name = link.second;
		if ( factory_options_.count( factory_name ) == 0 ) {
			TR.Error << factory_name << " linked to " << selector_name << " is not defined." << std::endl;
			error_counter++;
		} else {
			is_used[ factory_name ] = true;
			generated_popsize += core::Size( factory_options_.at( factory_name ).second.at( "size" ) );
		}
	}

	if ( generated_popsize < supported_population_size_ ) {
		TR.Error << "All offspring factories combined produce only " << generated_popsize << " new individuals, but " << supported_population_size_ << " are supported." << std::endl;
		error_counter++;
	}

	for ( const auto& usage : is_used ) {
		if ( !usage.second ) {
			TR.Warning << usage.first << " is defined but never used." << std::endl;
		}
	}

	// TODO This does not check if options provided are unsupported or wrong. Maybe I should add this.

	return error_counter;
}

core::Size EvolutionOptions::check_movers() const {

	core::Size error_counter = 0;

	std::map< std::string, bool > is_used;

	const utility::vector1< std::string > mover_types{ "start_from", "transform", "high_res_docker", "final_minimizer" };

	for ( const auto& mov_opt : mover_options_ ) {
		std::string mover_name = mov_opt.first;
		std::string mover_type = mov_opt.second.first;
		is_used[ mover_name ] = false;

		if ( std::find( mover_types.begin(), mover_types.end(), mover_type ) == mover_types.end() ) {
			TR.Error << mover_type << " is not supported. Available are " << mover_types;
			error_counter++;
		}

		if ( mover_type == "start_from" ) {
			if ( mov_opt.second.second.count( "x" ) == 0 ) {
				TR.Error << "x option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			}
			if ( mov_opt.second.second.count( "y" ) == 0 ) {
				TR.Error << "y option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			}
			if ( mov_opt.second.second.count( "z" ) == 0 ) {
				TR.Error << "z option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			}
		} else if ( mover_type == "transform" ) {
			if ( mov_opt.second.second.count( "box_size" ) == 0 ) {
				TR.Error << "box_size option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "box_size" ) < 1.0 ) {
				TR.Error << "box_size option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}

			if ( mov_opt.second.second.count( "max_move_distance" ) == 0 ) {
				TR.Error << "max_move_distance option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "max_move_distance" ) <= 0.0 ) {
				TR.Error << "max_move_distance option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}

			if ( mov_opt.second.second.count( "max_rotation_angle" ) == 0 ) {
				TR.Error << "max_rotation_angle option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "max_rotation_angle" ) <= 0.0 ) {
				TR.Error << "max_rotation_angle option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}

			if ( mov_opt.second.second.count( "cycles" ) == 0 ) {
				TR.Error << "cycles option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "cycles" ) < 1.0 ) {
				TR.Error << "cycles option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}

			if ( mov_opt.second.second.count( "temperature" ) == 0 ) {
				TR.Error << "temperature option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "temperature" ) <= 0.0 ) {
				TR.Error << "temperature option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}
		} else if ( mover_type == "high_res_docker" ) {

			if ( mov_opt.second.second.count( "cycles" ) == 0 ) {
				TR.Error << "cycles option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "cycles" ) <= 0.0 ) {
				TR.Error << "cycles option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}

			if ( mov_opt.second.second.count( "repack_every_nth" ) == 0 ) {
				TR.Error << "repack_every_nth option is not set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( mov_opt.second.second.at( "repack_every_nth" ) < 0.0 ) {
				TR.Error << "repack_every_nth option is set too small to be feasible for mover " << mover_name << "." << std::endl;
				error_counter++;
			}

		} else if ( mover_type == "final_minimizer" || mover_type == "high_res_docker" ) {
			if ( mover_sfx_links_.count( mover_name ) == 0 ) {
				TR.Error << "No score function is set for mover " << mover_name << "." << std::endl;
				error_counter++;
			} else if ( score_function_wts_.count( mover_sfx_links_.at( mover_name ) ) == 0 ) {
				TR.Error << mover_sfx_links_.at( mover_name ) << " is used by " << mover_name << ", but undefined." << std::endl;
				error_counter++;
			}
		}
	}

	for ( const std::string& mover_name : mover_protocol_ ) {
		if ( is_used.count( mover_name ) == 0 ) {
			TR.Error << mover_name << " is used in protocol, but is not defined." << std::endl;
			error_counter++;
		} else {
			is_used[ mover_name ] = true;
		}
	}

	if ( mover_protocol_.empty() ) {
		TR.Warning << "No movers will be applied. Make sure this is intended behavior." << std::endl;
	}

	for ( const auto& usage : is_used ) {
		if ( !usage.second ) {
			TR.Warning << usage.first << " is defined, but never used." << std::endl;
		}
	}

	// TODO This does not check if options provided are unsupported or wrong. Maybe I should add this.

	return error_counter;
}

core::Size EvolutionOptions::check_pop_init() const {

	core::Size error_counter = 0;
	core::Size combined_init_size = 0;

	for ( const auto& init_opt : pop_init_options_ ) {
		const std::string& type = init_opt.first;
		core::Size size = 0;
		if ( init_opt.second.count( "size" ) == 0 ) {
			TR.Error << "Size attribute is not defined for pop init type " << type << std::endl;
			error_counter++;
		} else if ( init_opt.second.at( "size" ) <= 0 ) {
			TR.Error << "Size attribute for pop init type " << type << " is less or equal than zero." << std::endl;
			error_counter++;
		} else {
			size = init_opt.second.at( "size" );
		}
		combined_init_size += size;
		if ( type == "random" ) {
			for ( const auto& random_opt : init_opt.second ) {
				if ( random_opt.first != "size" ) {
					TR.Error << "Random pop init only supports size attribute, not " << random_opt.first << std::endl;
					error_counter++;
				}
			}
		} else if ( type == "best_loaded" ) {
			if ( path_to_score_memory_.empty() ) {
				TR.Error << type << " requires scores to be loaded, but no path to memory was defined." << std::endl;
				error_counter++;
			}
			if ( init_opt.second.count( "selection" ) == 0 ) {
				TR.Error << "Selection attribute needs to be defined for init type " << type << std::endl;
				error_counter++;
			} else if ( init_opt.second.at( "selection" ) <= 0 ) {
				TR.Error << "Selection attribute needs to greater 0 for init type " << type << std::endl;
				error_counter++;
			} else if ( init_opt.second.at( "selection" ) <= size ) {
				TR.Error << "Selection attribute needs to greater than init size for init type " << type << std::endl;
				error_counter++;
			}
			for ( const auto& random_opt : init_opt.second ) {
				if ( random_opt.first != "size" && random_opt.first != "selection" ) {
					TR.Error << "Random pop init only supports size attribute, not " << random_opt.first << std::endl;
					error_counter++;
				}
			}
		} else {
			TR.Error << "Unknown pop init type " << type << std::endl;
			error_counter++;
		}
	}

	if ( combined_init_size < supported_population_size_ ) {
		TR.Error << "Combined initial population size has to be at least as high as the supported population size." << std::endl;
		error_counter++;
	}

	if ( combined_init_size <= 0 ) {
		TR.Error << "Combined initial population size has to be greater than zero." << std::endl;
		error_counter++;
	}

	return error_counter;
}

core::Size EvolutionOptions::get_max_generations() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return generations_;
}

core::Size EvolutionOptions::get_external_scoring() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return external_scoring_;
}

utility::vector1< std::string > EvolutionOptions::get_selector_names() const {
	utility::vector1< std::string > names;
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	for ( const auto& option : selector_options_ ) {
		names.emplace_back( option.first );
	}
	return names;
}

const std::string& EvolutionOptions::get_selector_type( const std::string& name ) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return selector_options_.at( name ).first;
}

core::Real EvolutionOptions::get_selector_parameter( const std::string& name, const std::string& parameter) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	if ( selector_options_.count( name ) == 0 ) {
		utility_exit_with_message( name + " is not defined as selector." );
	} else if ( selector_options_.at( name ).second.count( parameter ) == 0 ) {
		utility_exit_with_message( parameter + " for selector " + name + " is not set." );
	}
	return selector_options_.at( name ).second.at( parameter );
}

utility::vector1<std::string> EvolutionOptions::get_factory_names() const {
	utility::vector1< std::string > names;
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	for ( const auto& option : factory_options_ ) {
		names.emplace_back( option.first );
	}
	return names;
}

const std::string& EvolutionOptions::get_factory_type( const std::string& name ) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return factory_options_.at( name ).first;
}

core::Real EvolutionOptions::get_factory_parameter( const std::string& name, const std::string& parameter) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	if ( factory_options_.count( name ) == 0 ) {
		utility_exit_with_message( name + " is not defined as factory." );
	} else if ( factory_options_.at( name ).second.count( parameter ) == 0 ) {
		utility_exit_with_message( parameter + " for factory " + name + " is not set." );
	}
	return factory_options_.at( name ).second.at( parameter );
}

const std::string &EvolutionOptions::get_path_to_external_smiles() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return path_to_external_smiles_;
}

const std::string &EvolutionOptions::get_path_to_reactions() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return path_to_reactions_;
}

const std::string &EvolutionOptions::get_path_to_reagents() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return path_to_reagents_;
}

core::Size EvolutionOptions::get_population_supported_size() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return supported_population_size_;
}

const utility::vector1<std::pair<std::string, std::string> > &EvolutionOptions::get_selector_factory_links() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return selector_factory_links_;
}

const std::string &EvolutionOptions::get_main_selector() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return main_selector_;
}

utility::vector1<std::string> EvolutionOptions::get_scfx_names() const {
	utility::vector1< std::string > names;
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	for ( const auto& option : score_function_wts_ ) {
		names.emplace_back( option.first );
	}
	return names;
}

const std::string &EvolutionOptions::get_scfx_wts(const std::string &name) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return score_function_wts_.at( name );
}

const utility::vector1< std::pair< std::string, core::Real > >& EvolutionOptions::get_scfx_reweighs(const std::string &name) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return score_function_reweighs_.at( name );
}

const std::string &EvolutionOptions::get_main_scfx() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return scoring_function_;
}

const std::string &EvolutionOptions::get_path_score_memory() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return path_to_score_memory_;
}

const std::string &EvolutionOptions::get_ligand_chain() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return ligand_chain_;
}

core::Size EvolutionOptions::get_n_scoring_runs() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return score_runs_;
}

const std::string &EvolutionOptions::get_pose_dir_path() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return pose_dump_directory_;
}

core::Real EvolutionOptions::get_similarity_penalty() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return similarity_penalty_;
}

const std::string &EvolutionOptions::get_main_term() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return main_score_term_;
}

const utility::vector1<std::string> &EvolutionOptions::get_mover_protocol() const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return mover_protocol_;
}

const std::string &EvolutionOptions::get_mover_type(const std::string &name) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return mover_options_.at( name ).first;
}

core::Real EvolutionOptions::get_mover_parameter(const std::string &name, const std::string& param) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	if ( mover_options_.count( name ) == 0 ) {
		utility_exit_with_message( name + " is not defined as mover." );
	} else if ( mover_options_.at( name ).second.count( param ) == 0 ) {
		utility_exit_with_message( param + " for mover " + name + " is not set." );
	}
	return mover_options_.at( name ).second.at( param );
}

const std::string &EvolutionOptions::get_mover_scfx(const std::string &name) const {
	if ( !initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is not initialized." );
	}
	return mover_sfx_links_.at( name );
}

void EvolutionOptions::parse_option_file( const std::string& option_path ) {
	if ( initialized_ ) {
		utility_exit_with_message( "Evolutionary option system is already initialized, option file should not be parsed again." );
	}

	// if any of those two get set, all default settings for these are removed
	bool evolution_protocol_defined = false;
	bool mover_protocol_defined = false;
	bool pop_init_defined = false;

	core::Size error_counter = 0;
	boost::property_tree::ptree tag_tree;

	try {
		boost::property_tree::read_xml( option_path, tag_tree );
	} catch ( boost::property_tree::xml_parser_error& e ) {
		TR.Error << "Failed to parse xml file. " << e.what() << std::endl;
		utility_exit_with_message("Error whilst parsing evolutionary xml file.");
	} catch ( ... ) {
		utility_exit_with_message("Unexpected error whilst opening or parsing evolutionary xml file.");
	}

	// that is how you get the attributes
	// important - skip attribute thingy if tag.first == <xmlcomment>
	for ( const std::pair< const std::string, boost::property_tree::ptree >& tag : tag_tree ) {
		const std::string& option_name = tag.first;
		if ( option_name == "<xmlcomment>" ) {
			continue;
		}
		if ( option_name == "Evolution" ) {
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "generations" ) {
					generations_ = utility::string2Size( attribute_value );
				} else if ( attribute_name == "external_scoring" ) {
					external_scoring_ = utility::string2Size( attribute_value );
				} else if ( attribute_name == "ligand_chain" ) {
					ligand_chain_ = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter++;
				}
			}
		} else if ( option_name == "Population" ) {
			for ( const boost::property_tree::ptree::value_type &v: tag.second.get_child("<xmlattr>") ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "main_selector" ) {
					main_selector_ = attribute_value;
				} else if ( attribute_name == "supported_size" ) {
					supported_population_size_ = utility::string2Size(attribute_value);
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag "
						<< option_name << std::endl;
					error_counter++;
				}
			}
		} else if ( option_name == "PopulationInit" ) {
			if ( !pop_init_defined ) {
				pop_init_defined = true;
				TR.Warning << "Population initialization is defined. All defaults are removed." << std::endl;
				pop_init_options_.clear();
			}
			std::string type;
			std::map< std::string, core::Size > opt_map;
			for ( const boost::property_tree::ptree::value_type &v: tag.second.get_child("<xmlattr>") ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "init_type" ) {
					type = attribute_value;
				} else if ( attribute_name == "size" || attribute_name == "selection" ) {
					opt_map[ attribute_name ] = utility::string2Size(attribute_value);
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag "
						<< option_name << std::endl;
					error_counter++;
				}
			}
			pop_init_options_[ type ] = opt_map;
		} else if ( option_name == "Library" ) {
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "external_smiles" ) {
					path_to_external_smiles_ = attribute_value;
				} else if ( attribute_name == "reaction_file" ) {
					path_to_reactions_ = attribute_value;
				} else if ( attribute_name == "reagent_file" ) {
					path_to_reagents_ = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter++;
				}
			}
		} else if ( option_name == "ScoreFunction" ) {
			std::string scfx_name;
			std::string scfx_wts;
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "name" ) {
					scfx_name = attribute_value;
				} else if ( attribute_name == "weight_file" ) {
					scfx_wts = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag "
						<< option_name << std::endl;
					error_counter++;
				}
			}
			if ( scfx_name.empty() || scfx_wts.empty() ) {
				TR.Error << option_name << " requires both name and weight_file attribute to be set." << std::endl;
				error_counter++;
			} else {
				score_function_wts_[ scfx_name ] = scfx_wts;
			}
		} else if ( option_name == "ScoreFunctionReweigh" ) {
			std::string scfx_name;
			std::string term;
			std::string weight;
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "name" ) {
					scfx_name = attribute_value;
				} else if ( attribute_name == "term" ) {
					term = attribute_value;
				} else if ( attribute_name == "weight" ) {
					weight = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag "
						<< option_name << std::endl;
					error_counter++;
				}
			}
			if ( scfx_name.empty() || term.empty() || weight.empty() ) {
				TR.Error << option_name << " requires name, term and weight attribute to be set." << std::endl;
				error_counter++;
			} else {
				score_function_reweighs_[ scfx_name ].emplace_back( term, utility::string2Real( weight ) );
			}
		} else if ( option_name == "Scorer" ) {
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "n_runs" ) {
					score_runs_ = utility::string2Size( attribute_value );
				} else if ( attribute_name == "main_term" ) {
					main_score_term_ = attribute_value;
				} else if ( attribute_name == "pose_dir" ) {
					pose_dump_directory_ = attribute_value;
				} else if ( attribute_name == "similarity_penalty" ) {
					similarity_penalty_ = utility::string2Real( attribute_value );
				} else if ( attribute_name == "similarity_penalty_threshold" ) {
					similarity_penalty_threshold_ = utility::string2Real( attribute_value );
				} else if ( attribute_name == "score_function" ) {
					scoring_function_ = attribute_value;
				} else if ( attribute_name == "memory" ) {
					path_to_score_memory_ = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter++;
				}
			}
		} else if ( option_name == "Mover" ) {
			std::string name;
			std::string type;
			std::map< std::string, std::string > params;
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "name" ) {
					name = attribute_value;
				} else if ( attribute_name == "type" ) {
					type = attribute_value;
				} else {
					params[ attribute_name ] = attribute_value;
				}
			}
			if ( name.empty() || type.empty() ) {
				TR.Error << option_name << " requires name and type attribute to be set." << std::endl;
				error_counter++;
			} else {
				if ( type == "start_from" ) {
					mover_options_[ name ] = mover_options_.at( "std_start" );
				} else if ( type == "transform" ) {
					mover_options_[ name ] = mover_options_.at( "std_transform" );
				} else if ( type == "high_res_docker" ) {
					mover_options_[ name ] = mover_options_.at( "std_hres" );
				} else if ( type == "final_minimizer" ) {
					mover_options_[ name ] = mover_options_.at( "std_final_min" );
				} else {
					TR.Error << "Unknown type " << type << " for mover " << name << std::endl;
					error_counter++;
				}
				if ( params.empty() ) {
					TR.Warning << "Mover " << name << " receives " << type << " default parameters only." << std::endl;
				} else {
					for ( const std::pair< std::string, std::string > p : params ) {
						if ( p.first == "score_function" ) {
							mover_sfx_links_[ name ] = p.second;
						} else {
							mover_options_.at( name ).second[ p.first ] = utility::string2Real( p.second );
						}
					}
				}
			}
		} else if ( option_name == "Selector" ) {
			std::string name;
			std::string type;
			std::map< std::string, std::string > params;
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "name" ) {
					name = attribute_value;
				} else if ( attribute_name == "type" ) {
					type = attribute_value;
				} else {
					params[ attribute_name ] = attribute_value;
				}
			}
			if ( name.empty() || type.empty() ) {
				TR.Error << option_name << " requires name and type attribute to be set." << std::endl;
				error_counter++;
			} else {
				if ( type == "elitist" ) {
					selector_options_[ name ] = selector_options_.at( "std_elitist" );
				} else if ( type == "tournament" ) {
					selector_options_[ name ] = selector_options_.at( "std_tournament" );
				} else if ( type == "roulette" ) {
					selector_options_[ name ] = selector_options_.at( "std_roulette" );
				} else {
					TR.Error << "Unknown type " << type << " for selector " << name << std::endl;
					error_counter++;
				}
				if ( params.empty() ) {
					TR.Warning << "Selector " << name << " receives " << type << " default parameters only." << std::endl;
				} else {
					for ( const std::pair< std::string, std::string > p : params ) {
						if ( p.first == "remove" ) {
							selector_options_.at( name ).second[ p.first ] = ( p.second == "True" );
						} else {
							selector_options_.at(name).second[p.first] = utility::string2Real(p.second);
						}
					}
				}
			}
		} else if ( option_name == "Factory" ) {
			std::string name;
			std::string type;
			std::map< std::string, std::string > params;
			for ( const boost::property_tree::ptree::value_type& v : tag.second.get_child( "<xmlattr>" ) ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "name" ) {
					name = attribute_value;
				} else if ( attribute_name == "type" ) {
					type = attribute_value;
				} else {
					params[ attribute_name ] = attribute_value;
				}
			}
			if ( name.empty() || type.empty() ) {
				TR.Error << option_name << " requires name and type attribute to be set." << std::endl;
				error_counter++;
			} else {
				if ( type == "mutator" ) {
					factory_options_[ name ] = factory_options_.at( "std_mutator" );
				} else if ( type == "crossover" ) {
					factory_options_[ name ] = factory_options_.at( "std_crossover" );
				} else if ( type == "identity" ) {
					factory_options_[ name ] = factory_options_.at( "std_identity" );
				} else {
					TR.Error << "Unknown type " << type << " for factory " << name << std::endl;
					error_counter++;
				}
				if ( params.empty() ) {
					TR.Warning << "Factory " << name << " receives " << type << " default parameters only." << std::endl;
				} else {
					for ( const std::pair< std::string, std::string > p : params ) {
						factory_options_.at( name ).second[ p.first ] = utility::string2Real( p.second );
					}
				}
			}
		} else if ( option_name == "EvolutionProtocol" ) {
			if ( !evolution_protocol_defined ) {
				evolution_protocol_defined = true;
				selector_factory_links_.clear();
				TR.Warning << option_name << " is defined. All default settings for it are deleted." << std::endl;
			}
			std::string selector;
			std::string factory;
			for ( const boost::property_tree::ptree::value_type &v: tag.second.get_child("<xmlattr>") ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "selector" ) {
					selector = attribute_value;
				} else if ( attribute_name == "factory" ) {
					factory = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter++;
				}
			}
			if ( selector.empty() || factory.empty() ) {
				TR.Error << option_name << " requires selector and factory attribute to be set." << std::endl;
				error_counter++;
			} else {
				selector_factory_links_.emplace_back( selector, factory );
			}
		} else if ( option_name == "MoverProtocol" ) {
			if ( !mover_protocol_defined ) {
				mover_protocol_defined = true;
				mover_protocol_.clear();
				TR.Warning << option_name << " is defined. All default settings for it are deleted." << std::endl;
			}
			std::string mover;
			for ( const boost::property_tree::ptree::value_type &v: tag.second.get_child("<xmlattr>") ) {
				const std::string attribute_name = v.first;
				const std::string attribute_value = v.second.data();
				if ( attribute_name == "mover" ) {
					mover = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter++;
				}
			}
			if ( mover.empty() ) {
				TR.Error << option_name << " requires mover attribute to be set." << std::endl;
				error_counter++;
			} else {
				mover_protocol_.emplace_back( mover );
			}
		} else {
			TR.Error << "Unable to parse tag " << option_name << std::endl;
			error_counter++;
		}
	}

	if ( error_counter != 0 ) {
		TR.Error << "Found a total of " << error_counter << " errors whilst parsing file " << option_path << ". See details above." << std::endl;
		utility_exit_with_message( "Unable to parse option file." );
	}

}

const std::map< std::string, std::map< std::string, core::Size > >& EvolutionOptions::get_pop_init_options() const {
	return pop_init_options_;
}

core::Real EvolutionOptions::get_similarity_penalty_threshold() const {
	return similarity_penalty_threshold_;
}
}
}
