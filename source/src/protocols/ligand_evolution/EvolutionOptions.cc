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
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/ligand_evolution.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
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
    parse_cmdline();
	check_validity();
    TR.Debug << "Options passed all checks." << std::endl;
}

EvolutionOptions::EvolutionOptions( std::string const& path_to_option_file) {

	utility::io::izstream option_file( path_to_option_file );
	if ( !option_file ) {
		TR.Error << "Can't find evolutionary option file at " << path_to_option_file << std::endl;
		utility_exit_with_message( "Unable to open option file." );
	}

    parse_cmdline();
	parse_option_file( path_to_option_file );

    // TODO I removed some functionality from the evo option file to commandline. Make sure this is reflected here.
	check_validity();

	TR.Debug << "Options passed all checks." << std::endl;
}

void EvolutionOptions::parse_cmdline() {

    if ( !basic::options::option[ basic::options::OptionKeys::parser::protocol ].active() ) {
        TR.Error << "No docking protocol provided. Use -parser:protocol to define a path to a rosetta script" << std::endl;
        error_counter_++;
    } else {
        protocol_path_ = basic::options::option[ basic::options::OptionKeys::parser::protocol ];
        utility::io::izstream protocol_file( protocol_path_ );
        if ( !protocol_file ) {
            TR.Error << "Can't find protocol file " << protocol_path_ << "." << std::endl;
            error_counter_++;
        }
    }

    pose_stream_ = core::import_pose::pose_stream::streams_from_cmd_line();
    if ( !pose_stream_.has_another_pose() ) {
        TR.Error << "No pdb provided. Please use one of the commandline input functions like -in::file::s" << std::endl;
        error_counter_++;
    }

    if ( !basic::options::option[ basic::options::OptionKeys::ligand_evolution::xyz ].active() ){
        TR.Error << "No position for ligand placement specified. Use -ligand_evolution::xyz to define a position" << std::endl;
        error_counter_++;
    } else {
        TR.Debug << "Passed start position: "
                 << " " << basic::options::option[basic::options::OptionKeys::ligand_evolution::xyz][1]
                 << " " << basic::options::option[basic::options::OptionKeys::ligand_evolution::xyz][2]
                 << " " << basic::options::option[basic::options::OptionKeys::ligand_evolution::xyz][3]
                 << std::endl;
        xyz_.x(basic::options::option[basic::options::OptionKeys::ligand_evolution::xyz][1]);
        xyz_.y(basic::options::option[basic::options::OptionKeys::ligand_evolution::xyz][2]);
        xyz_.z(basic::options::option[basic::options::OptionKeys::ligand_evolution::xyz][3]);
    }

    if ( basic::options::option[ basic::options::OptionKeys::ligand_evolution::external_scoring ].active() != basic::options::option[ basic::options::OptionKeys::ligand_evolution::smiles_file ].active() ){
        TR.Error << "Only one external scoring option is provided. Make sure you specified external_scoring and smiles_file" << std::endl;
        error_counter_++;
    } else if ( basic::options::option[ basic::options::OptionKeys::ligand_evolution::external_scoring ].active() && basic::options::option[ basic::options::OptionKeys::ligand_evolution::smiles_file ].active() ) {
        external_scoring_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::external_scoring ];
        path_to_external_smiles_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::smiles_file ].value();
    }

    if ( !basic::options::option[ basic::options::OptionKeys::ligand_evolution::main_scfx ].active() ){
        TR.Error << "No main scoring function defined. Use -ligand_evolution::main_scfx to define a position." << std::endl;
        error_counter_++;
    } else {
        scoring_function_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::main_scfx ];
    }

    if ( !basic::options::option[ basic::options::OptionKeys::ligand_evolution::reaction_file ].active() ){
        if ( external_scoring_ == 0 ) {
            TR.Error << "No reaction file defined. Use -ligand_evolution::reagent_file to specify your library."
                     << std::endl;
            error_counter_++;
        }
    } else {
        path_to_reactions_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::reaction_file ].value();
    }

    if ( !basic::options::option[ basic::options::OptionKeys::ligand_evolution::reagent_file ].active() ){
        if ( external_scoring_ == 0 ) {
            TR.Error << "No reagent file defined. Use -ligand_evolution::reagent_file to specify your library."
                     << std::endl;
            error_counter_++;
        }
    } else {
        path_to_reagents_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::reagent_file ].value();
    }

    // defaults are set for these options
    score_runs_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::n_scoring_runs ];
    ligand_chain_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::ligand_chain ];
    pose_dump_directory_ = basic::options::option[ basic::options::OptionKeys::ligand_evolution::pose_output_directory ].value();

    // these options are optional
    if ( basic::options::option[ basic::options::OptionKeys::ligand_evolution::score_mem_path ].active() ) {
        path_to_score_memory_ = basic::options::option[basic::options::OptionKeys::ligand_evolution::score_mem_path].value();
    }

}

void EvolutionOptions::check_validity() {

	if ( external_scoring_ > 0 ) {

		utility::io::izstream ligand_smiles_file( path_to_external_smiles_ );
		if ( !ligand_smiles_file ) {
			TR.Error << "Can't find ligand smiles list to score at " << path_to_external_smiles_ << "." << std::endl;
			error_counter_++;
		}

	} else {

		if ( generations_ <= 0 ) {
			TR.Error << "The algorithm is required to run for at least one generation." << std::endl;
			error_counter_++;
		}

		utility::io::izstream reaction_file( path_to_reactions_ );
		if ( !reaction_file ) {
			TR.Error << "Can't find reactions file at " << path_to_reactions_ << "." << std::endl;
			error_counter_++;
		}

		utility::io::izstream reagents_file( path_to_reagents_ );
		if ( !reagents_file ) {
			TR.Error << "Can't find reagents file at " << path_to_reagents_ << "." << std::endl;
			error_counter_++;
		}

		if ( main_selector_.empty() ) {
			TR.Error << "No main selector provided." << std::endl;
			error_counter_++;
		}

		if ( supported_population_size_ <= 0 ) {
			TR.Error << "Supported population size has to be greater than zero." << std::endl;
			error_counter_++;
		}

		error_counter_ += check_selectors();
		error_counter_ += check_factories();
		error_counter_ += check_pop_init();

	}

	error_counter_ += check_scorer_setup();

	if ( error_counter_ != 0 ) {
		utility_exit_with_message( "Found a total of " + std::to_string( error_counter_ ) + " errors. See details above." );
	}

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
	}

	return error_counter;
}

core::Size EvolutionOptions::check_selectors() const {

	core::Size error_counter = 0;

	std::map< std::string, bool > is_used;

	utility::vector1< std::string > const selector_types{ "elitist", "roulette", "tournament" };

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
	for ( auto const& link : selector_factory_links_ ) {
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

	for ( auto const& usage : is_used ) {
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

	utility::vector1< std::string > const factory_types{ "mutator", "crossover", "identity" };

	for ( auto const& fac_ops : factory_options_ ) {

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
	for ( auto const& link : selector_factory_links_ ) {
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

	for ( auto const& usage : is_used ) {
		if ( !usage.second ) {
			TR.Warning << usage.first << " is defined but never used." << std::endl;
		}
	}

	// TODO This does not check if options provided are unsupported or wrong. Maybe I should add this.

	return error_counter;
}

core::Size EvolutionOptions::check_pop_init() const {

	core::Size error_counter = 0;
	core::Size combined_init_size = 0;

	for ( auto const& init_opt : pop_init_options_ ) {
		std::string const& type = init_opt.first;
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
			for ( auto const& random_opt : init_opt.second ) {
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
			for ( auto const& random_opt : init_opt.second ) {
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
	return generations_;
}

core::Size EvolutionOptions::get_external_scoring() const {
	return external_scoring_;
}

utility::vector1< std::string > EvolutionOptions::get_selector_names() const {
	utility::vector1< std::string > names;
	for ( auto const& option : selector_options_ ) {
		names.emplace_back( option.first );
	}
	return names;
}

std::string const& EvolutionOptions::get_selector_type( std::string const& name ) const {
	return selector_options_.at( name ).first;
}

core::Real EvolutionOptions::get_selector_parameter( std::string const& name, std::string const& parameter) const {
	if ( selector_options_.count( name ) == 0 ) {
		utility_exit_with_message( name + " is not defined as selector." );
	} else if ( selector_options_.at( name ).second.count( parameter ) == 0 ) {
		utility_exit_with_message( parameter + " for selector " + name + " is not set." );
	}
	return selector_options_.at( name ).second.at( parameter );
}

utility::vector1<std::string> EvolutionOptions::get_factory_names() const {
	utility::vector1< std::string > names;
	for ( auto const& option : factory_options_ ) {
		names.emplace_back( option.first );
	}
	return names;
}

std::string const& EvolutionOptions::get_factory_type( std::string const& name ) const {
	return factory_options_.at( name ).first;
}

core::Real EvolutionOptions::get_factory_parameter( std::string const& name, std::string const& parameter) const {
	if ( factory_options_.count( name ) == 0 ) {
		utility_exit_with_message( name + " is not defined as factory." );
	} else if ( factory_options_.at( name ).second.count( parameter ) == 0 ) {
		utility_exit_with_message( parameter + " for factory " + name + " is not set." );
	}
	return factory_options_.at( name ).second.at( parameter );
}

std::string const& EvolutionOptions::get_path_to_external_smiles() const {
	return path_to_external_smiles_;
}

std::string const& EvolutionOptions::get_path_to_reactions() const {
	return path_to_reactions_;
}

std::string const& EvolutionOptions::get_path_to_reagents() const {
	return path_to_reagents_;
}

core::Size EvolutionOptions::get_population_supported_size() const {
	return supported_population_size_;
}

utility::vector1<std::pair<std::string, std::string> > const& EvolutionOptions::get_selector_factory_links() const {
	return selector_factory_links_;
}

std::string const& EvolutionOptions::get_main_selector() const {
	return main_selector_;
}

std::string const& EvolutionOptions::get_main_scfx() const {
	return scoring_function_;
}

std::string const& EvolutionOptions::get_path_score_memory() const {
	return path_to_score_memory_;
}

std::string const& EvolutionOptions::get_ligand_chain() const {
	return ligand_chain_;
}

core::Size EvolutionOptions::get_n_scoring_runs() const {
	return score_runs_;
}

std::string const& EvolutionOptions::get_pose_dir_path() const {
	return pose_dump_directory_;
}

core::Real EvolutionOptions::get_similarity_penalty() const {
	return similarity_penalty_;
}

std::string const& EvolutionOptions::get_main_term() const {
	return main_score_term_;
}

void EvolutionOptions::parse_option_file( std::string const& option_path ) {

	// if any of those two get set, all default settings for these are removed
	bool evolution_protocol_defined = false;
	bool mover_protocol_defined = false;
	bool pop_init_defined = false;
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
	for ( std::pair< std::string const, boost::property_tree::ptree > const& tag : tag_tree ) {
		std::string const& option_name = tag.first;
		if ( option_name == "<xmlcomment>" ) {
			continue;
		}
		if ( option_name == "Evolution" ) {
			for ( boost::property_tree::ptree::value_type const& v : tag.second.get_child( "<xmlattr>" ) ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
				if ( attribute_name == "generations" ) {
					generations_ = utility::string2Size( attribute_value );
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter_++;
				}
			}
		} else if ( option_name == "Population" ) {
			for ( boost::property_tree::ptree::value_type const& v: tag.second.get_child("<xmlattr>") ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
				if ( attribute_name == "main_selector" ) {
					main_selector_ = attribute_value;
				} else if ( attribute_name == "supported_size" ) {
					supported_population_size_ = utility::string2Size(attribute_value);
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag "
						<< option_name << std::endl;
					error_counter_++;
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
			for ( boost::property_tree::ptree::value_type const& v: tag.second.get_child("<xmlattr>") ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
				if ( attribute_name == "init_type" ) {
					type = attribute_value;
				} else if ( attribute_name == "size" || attribute_name == "selection" ) {
					opt_map[ attribute_name ] = utility::string2Size(attribute_value);
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag "
						<< option_name << std::endl;
					error_counter_++;
				}
			}
			pop_init_options_[ type ] = opt_map;
		} else if ( option_name == "Scorer" ) {
			for ( boost::property_tree::ptree::value_type const& v : tag.second.get_child( "<xmlattr>" ) ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
				if ( attribute_name == "main_term" ) {
					main_score_term_ = attribute_value;
				} else if ( attribute_name == "similarity_penalty" ) {
					similarity_penalty_ = utility::string2Real( attribute_value );
				} else if ( attribute_name == "similarity_penalty_threshold" ) {
					similarity_penalty_threshold_ = utility::string2Real( attribute_value );
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter_++;
				}
			}
		} else if ( option_name == "Selector" ) {
			std::string name;
			std::string type;
			std::map< std::string, std::string > params;
			for ( boost::property_tree::ptree::value_type const& v : tag.second.get_child( "<xmlattr>" ) ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
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
				error_counter_++;
			} else {
				if ( type == "elitist" ) {
					selector_options_[ name ] = selector_options_.at( "std_elitist" );
				} else if ( type == "tournament" ) {
					selector_options_[ name ] = selector_options_.at( "std_tournament" );
				} else if ( type == "roulette" ) {
					selector_options_[ name ] = selector_options_.at( "std_roulette" );
				} else {
					TR.Error << "Unknown type " << type << " for selector " << name << std::endl;
					error_counter_++;
				}
				if ( params.empty() ) {
					TR.Warning << "Selector " << name << " receives " << type << " default parameters only." << std::endl;
				} else {
					for ( std::pair< std::string, std::string > const& p : params ) {
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
			for ( boost::property_tree::ptree::value_type const& v : tag.second.get_child( "<xmlattr>" ) ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
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
				error_counter_++;
			} else {
				if ( type == "mutator" ) {
					factory_options_[ name ] = factory_options_.at( "std_mutator" );
				} else if ( type == "crossover" ) {
					factory_options_[ name ] = factory_options_.at( "std_crossover" );
				} else if ( type == "identity" ) {
					factory_options_[ name ] = factory_options_.at( "std_identity" );
				} else {
					TR.Error << "Unknown type " << type << " for factory " << name << std::endl;
					error_counter_++;
				}
				if ( params.empty() ) {
					TR.Warning << "Factory " << name << " receives " << type << " default parameters only." << std::endl;
				} else {
					for ( std::pair< std::string, std::string > const& p : params ) {
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
			for ( boost::property_tree::ptree::value_type const& v: tag.second.get_child("<xmlattr>") ) {
				std::string const& attribute_name = v.first;
				std::string const& attribute_value = v.second.data();
				if ( attribute_name == "selector" ) {
					selector = attribute_value;
				} else if ( attribute_name == "factory" ) {
					factory = attribute_value;
				} else {
					TR.Error << "Unable to parse " << attribute_name << "=\"" << attribute_value << "\" for tag " << option_name << std::endl;
					error_counter_++;
				}
			}
			if ( selector.empty() || factory.empty() ) {
				TR.Error << option_name << " requires selector and factory attribute to be set." << std::endl;
				error_counter_++;
			} else {
				selector_factory_links_.emplace_back( selector, factory );
			}
		} else {
			TR.Error << "Unable to parse tag " << option_name << std::endl;
			error_counter_++;
		}
	}

	if ( error_counter_ != 0 ) {
		TR.Error << "Found a total of " << error_counter_ << " errors whilst parsing file " << option_path << ". See details above." << std::endl;
		utility_exit_with_message( "Unable to parse option file." );
	}

}

std::map< std::string, std::map< std::string, core::Size > > const& EvolutionOptions::get_pop_init_options() const {
	return pop_init_options_;
}

core::Real EvolutionOptions::get_similarity_penalty_threshold() const {
	return similarity_penalty_threshold_;
}

    std::string const& EvolutionOptions::get_protocol_path() const {
        return protocol_path_;
    }

    core::pose::PoseOP EvolutionOptions::get_pose_from_stream() {
        core::pose::PoseOP pose( new core::pose::Pose );
        pose_stream_.fill_pose(*pose);
        std::string pose_descriptor = pose_stream_.get_last_pose_descriptor_string();
        if( pose_stream_.has_another_pose() ) {
            TR.Warning << "More than one input structures are described. Selecting the first one encountered:" << std::endl;
            TR.Warning << pose_descriptor << std::endl;
        }
        return pose;
    }

    numeric::xyzVector<core::Real> EvolutionOptions::get_start_xyz() const {
        return xyz_;
    }
}
}
