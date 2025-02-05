// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Population.cc
/// @brief  Implementation of the %Population class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/Population.hh>

// package headers
#include <protocols/ligand_evolution/FragmentLibrary.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <numeric/random/WeightedReservoirSampler.hh>

// C/C++ headers
#include <algorithm>

static basic::Tracer TR( "protocols.ligand_evolution.Population" );


namespace protocols {
namespace ligand_evolution {

utility::vector1< Individual > const& Population::individuals() const {
	return individuals_;
}

void Population::sort() {
	assert( individuals_.size() > 0 );
	if ( is_sorted() ) {
		return;
	}
	struct {
		bool operator()( Individual const& a, Individual const& b ) const {
			return a.score() < b.score();
		}
	} lowerEnergy;
	std::sort( individuals_.begin(), individuals_.end(), lowerEnergy );
	sorting_guaranteed_ = true;
}

core::Size Population::size() const {
	return individuals_.size();
}

Individual const& Population::individual( core::Size index ) const {
	assert( individuals_.size() > 0 );
	return individuals_[ index ];
}

void Population::add_random( core::Size n_random_individuals, FragmentLibrary const& lib ) {
	for ( core::Size ii( 1 ); ii <= n_random_individuals; ++ii ) {
		Individual individual( lib.random_ligand(), { 0 }, "random" );
		add_individual( individual );
	}
}

void Population::add_individuals(utility::vector1< LigandIdentifier > const& new_individuals ) {
	utility::vector1< Individual > indis;
	for ( LigandIdentifier const& id : new_individuals ) {
		indis.push_back( Individual( id, { 0 }, "manually added" ) );
	}
	add_individuals( indis );
}

void Population::add_individuals( utility::vector1< Individual > const& new_individuals ) {
	for ( Individual const& individual : new_individuals ) {
		add_individual( individual );
	}
	TR.Debug << "Added " << size() << " new individuals" << std::endl;
}

void Population::next_generation( Selector const& selector ) {

	TR.Debug << "Ended generation " << generation_ << " with " << size() << " individuals. " << std::endl;
	log_generation();
	replace_population( selector.apply(*this, supported_size_, false ) );
	log_generation();
	++generation_;
	sort();

	TR.Debug << "Ready for generation " << generation_ << " with " << size() << " individuals." << std::endl;
}

core::Size Population::generation() const {
	return generation_;
}

bool Population::check_sorting() {
	sorting_guaranteed_ = true;
	core::Real last_score = -99999999.9;
	// iterate over all individuals and check if they are in order
	for ( Individual const& individual : individuals_ ) {
		if ( !individual.is_scored() ) {
			// an unscored individual can't be sorted
			sorting_guaranteed_ = false;
			break;
		} else if ( individual.score() >= last_score ) {
			// this is fine
			last_score = individual.score();
		} else {
			// the individual must have a lower and therefore better score than the previous one
			sorting_guaranteed_ = false;
			break;
		}
	}
	return sorting_guaranteed_;
}

bool Population::is_sorted() {
	if ( sorting_guaranteed_ ) {
		return true;
	} else {
		return check_sorting();
	}
}

void Population::add_individual( Individual const& individual ) {
	sorting_guaranteed_ = false;
	individuals_.push_back( individual );
	// try to set the id of the individual. This will be successful only if no id has been set before and return true if successful
	if ( individuals_.back().id( next_id_ ) ) {
		expand_inheritance_graph( individuals_.back() );
		next_id_++;
	}
}

utility::vector1<Individual>& Population::individuals() {
	sorting_guaranteed_ = false;
	return individuals_;
}

utility::vector1< Individual > Population::remove_individuals( utility::vector1< core::Size > indices ) {

	if ( indices.empty() ) {
		TR.Warning << "Empty list of indices was hand over to remove. Do nothing." << std::endl;
		return utility::vector1< Individual >();
	}

	std::sort( indices.begin(), indices.end() );

	utility::vector1< Individual > return_individuals;
	utility::vector1< Individual > remaining_individuals;

	core::Size current_index = 1;
	for ( core::Size ii( 1 ); ii <= size(); ++ii ) {
		if ( ii == indices[ current_index ] ) {
			return_individuals.push_back( individuals_[ ii ] );
			current_index = std::min( current_index + 1, indices.size() ) ;
		} else {
			remaining_individuals.push_back( individuals_[ ii ] );
		}
	}
	individuals_ = remaining_individuals;

	return return_individuals;
}

void Population::replace_population( utility::vector1< Individual > const& individuals ) {
	individuals_.clear();
	add_individuals( individuals );
}

utility::vector1< Individual > Population::individuals( utility::vector1< core::Size > const& indices ) {

	utility::vector1< Individual > return_individuals;
	for ( core::Size ii : indices ) {
		return_individuals.push_back( individuals_[ ii ] );
	}

	return return_individuals;
}

void Population::log_generation() {
	sort();
	utility::vector1< core::Size > ids( size() );
	for ( core::Size ii( 1 ); ii <= size(); ++ii ) {
		ids.at( ii ) = individuals_.at( ii ).id();
	}
	generation_log_.push_back( ids );
}

utility::vector1< utility::vector1< core::Size > > const& Population::expose_generation_log() const {
	return generation_log_;
}

void Population::expand_inheritance_graph( Individual const& individual ) {
	for ( core::Size parent : individual.parents() ) {
		inheritance_graph_.emplace_back( individual.id(), parent );
	}
}

utility::vector1< std::pair< core::Size, core::Size > > const& Population::expose_inheritance_graph() const {
	return inheritance_graph_;
}

void Population::set_supported_size(core::Size supported_size) {
	if ( supported_size_ != 0 ) {
		TR.Warning << "Supported size for population gets changed after being set before. This should not happen!" << std::endl;
	}
	supported_size_ = supported_size;
}

Individual &Population::individual(core::Size index) {
	sorting_guaranteed_ = false;
	return individuals_[ index ];
}

void Population::initialize_from_evotoptions( EvolutionOptions const& options, FragmentLibrary const& library, Scorer const& scorer ) {
	set_supported_size(options.get_population_supported_size());
	for ( std::pair<std::string const, std::map<std::string, core::Size> > const& init_opt: options.get_pop_init_options() ) {
		std::string const &init_type = init_opt.first;
		std::map<std::string, core::Size> const &type_options = init_opt.second;
		if ( init_type == "random" ) {
			add_random(type_options.at("size"), library);
		} else if ( init_type == "best_loaded" ) {
			utility::vector1<LigandIdentifier> best_individuals = scorer.get_best_loaded(
				type_options.at("selection"));
			numeric::random::WeightedReservoirSampler<LigandIdentifier> sampler(type_options.at("size"));
			for ( LigandIdentifier const &indi: best_individuals ) {
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
				add_individuals(best_individuals);
			}
		}
	}
}

}
}
