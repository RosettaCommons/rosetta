// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file GeneticAlgorithm.hh
/// @brief template class for genetic algorithm protocols
/// @author ashworth, based on template "pseudo"code by Colin Smith

#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>

#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/genetic_algorithm/FitnessFunction.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

#include <numeric/random/random.fwd.hh>

#include <boost/unordered_map.hpp>

#include <algorithm> // std::copy
#include <iostream>
#include <set>

#include <utility/exit.hh>


namespace protocols {
namespace genetic_algorithm {

static thread_local basic::Tracer TR( "protocols.genetic_algorithm" );

GeneticAlgorithmBase::GeneticAlgorithmBase() :
	utility::pointer::ReferenceCount(),
	entity_randomizer_(0),
	entity_template_(0),
	current_generation_(1),
	max_generations_(0),
	max_population_size_(0),
	number_to_propagate_(1),
	fraction_by_recombination_(0.5),
	checkpoint_prefix_(""),
	checkpoint_write_interval_(0),
	checkpoint_gzip_(false),
	checkpoint_rename_(false)
{}

GeneticAlgorithmBase::~GeneticAlgorithmBase()
{
	// ensure that checkpoint files are not accidentally reused
	if (checkpoint_rename_) rename_checkpoint_files();
}

Entity::OP
GeneticAlgorithmBase::add_entity( EntityElements const & traits )
{
	EntityOP entity;
	TraitEntityHashMap::iterator hash_it( entity_cache_.find( traits ) );
	//Vec1Hash h;
	if ( hash_it == entity_cache_.end() ) {
		// make a new entity if it does not exist in the cache
		entity = new_entity();
		entity->set_traits( traits );
		// add the new entity to the cache
		entity_cache_[entity->traits()] = entity;
		//TR << "(traits) Adding new entity " << *entity << " " << h( traits ) << std::endl;
	} else {
		// otherwise use the cached entity
		entity = hash_it->second;
		//TR << "(traits) Using cached entity" << h( traits ) << std::endl;
	}
	//TR << "Adding entity " << entity.get() << " to the current generation " << std::endl;
	generations_[current_generation_].push_back( entity );
	return entity;
}

Entity::OP
GeneticAlgorithmBase::add_entity( EntityOP entity )
{
	runtime_assert( entity != 0 );
	TraitEntityHashMap::iterator hash_it( entity_cache_.find( entity->traits() ) );
	//Vec1Hash h;
	if ( hash_it == entity_cache_.end() ) {
		// add the entity to the cache if it does not exist
		entity_cache_[entity->traits()] = entity;
		//TR << "(entity) Adding new entity " << *entity << " " << h( entity->traits() ) << std::endl;
	} else {
		// otherwise substitute the equivalent cached entity for the one given
		entity = hash_it->second;
		//TR << "(entity) Using cached entity" << h( entity->traits() ) << std::endl;
	}
	//TR << "Adding entity " << entity.get() << " to the current generation " << std::endl;
	generations_[current_generation_].push_back( entity );
	return entity;
}

Entity::OP
GeneticAlgorithmBase::add_parent_entity( EntityElements const & traits )
{
	EntityOP entity;
	TraitEntityHashMap::iterator hash_it( entity_cache_.find( traits ) );
	if ( hash_it == entity_cache_.end() ) {
		// make a new entity if it does not exist in the cache
		entity = new_entity();
		entity->set_traits( traits );
		//std::cout << __LINE__ << " Parent entity was not found in entity_cache_" << std::endl;
		entity_cache_[ traits ] = entity;
	} else {
		// otherwise use the cached entity
		entity = hash_it->second;
	}
	parent_entities_.push_back( entity );
	return entity;
}

Entity::OP
GeneticAlgorithmBase::add_parent_entity( EntityOP entity )
{
	runtime_assert( entity != 0 );
	TraitEntityHashMap::iterator hash_it( entity_cache_.find( entity->traits() ) );
	if ( hash_it != entity_cache_.end() ) {
		// use an equivalent hashed entity if it exists
		entity = hash_it->second;
		//std::cout << __LINE__ << " Parent entity was not found in entity_cache_" << std::endl;
		entity_cache_[ entity->traits() ] = entity;
	}
	parent_entities_.push_back( entity );
	return entity;
}

void GeneticAlgorithmBase::clear_parents() { parent_entities_.clear(); }


void
GeneticAlgorithmBase::add_parents_from_current_generation()
{
	runtime_assert(generations_[current_generation_].size());

	parent_entities_.insert(
		parent_entities_.end(),
		generations_[current_generation_].begin(),
		generations_[current_generation_].end()
	);
}

///@brief add the best entities from the previous generation
void
GeneticAlgorithmBase::propagate_best_from_previous_generation( core::Size size /* = 1 */, bool unique /* = true */ )
{
	runtime_assert(current_generation_ > 1);

	utility::vector1< EntityOP > sorted_entities(generations_[current_generation_-1]);
	std::sort( sorted_entities.begin(), sorted_entities.end(), lt_OP_deref< Entity > );

	// this will hopefully help keep population diversity higher
	if (unique && size > 1) {
		pop_iter new_last(std::unique(sorted_entities.begin(), sorted_entities.end(), eq_OP_deref< Entity >));
		sorted_entities.erase(new_last, sorted_entities.end());
	}

	for (core::Size i = 1; i <= size && i <= sorted_entities.size(); ++i) {
		add_entity(sorted_entities[i]);
	}
}

core::Real
GeneticAlgorithmBase::best_fitness_from_current_generation() const
{
	core::Real best( 0.0 ); bool first_found = false;
	for ( Size ii = 1; ii <= generations_[ current_generation_].size(); ++ii ) {
		if ( generations_[ current_generation_ ][ ii ] ) {
			if ( ! first_found || best > generations_[ current_generation_ ][ ii ]->fitness() ) {
				best = generations_[ current_generation_ ][ ii ]->fitness();
				first_found = true;
			}
		}
	}
	return best;
}

void GeneticAlgorithmBase::set_rand( EntityRandomizerOP r ) { entity_randomizer_ = r; }

void
GeneticAlgorithmBase::set_max_generations( core::Size s )
{
	max_generations_ = s;
	if (generations_.size() < max_generations_) generations_.resize(max_generations_);
}

void
GeneticAlgorithmBase::fill_with_random_entities( core::Size size /* = 0 */ )
{
	if ( size == 0 ) size = max_population_size_;
	while ( generations_[current_generation_].size() < size ) {
		add_entity( entity_randomizer_->random_entity() );
	}
}

void
GeneticAlgorithmBase::fill_with_perturbations_of_existing_entities( core::Size size /* = 0 */ )
{
	core::Size start_size = generations_[ current_generation_ ].size();
	for ( core::Size ii = 1; ii <= start_size; ++ii ) {
		TR << "fill_with_perturbations_of_existing_entities " << ii;
		generations_[ current_generation_ ][ ii ]->show( TR );
		TR << std::endl;
	}
	if ( size == 0 ) size = max_population_size_;
	TR << "creating an additional " << size - generations_[ current_generation_ ].size() << " sequences" << std::endl;

	while ( generations_[ current_generation_ ].size() < size ) {
		core::Size seed_sequence = numeric::random::uniform() * start_size + 1;
		EntityOP child = generations_[ current_generation_ ][ seed_sequence ]->clone();
		entity_randomizer_->mutate( *child );
		TR << "fill_with_perturbations_of_existing_entities ";
		child->show( TR );
		TR << std::endl;
		add_entity( child );
	}
}

///@brief add entities that are recombinants of fit parents
void
GeneticAlgorithmBase::fill_by_crossover( core::Size size /* = 0 */ )
{
	runtime_assert(parent_entities_.size());

	if ( size == 0 ) size = max_population_size_;

	while ( generations_[current_generation_].size() < size ) {
		// pick two random parents
		EntityOP child1 = tournament_select( parent_entities_ ).clone();
		EntityOP child2 = tournament_select( parent_entities_ ).clone();
		entity_randomizer_->crossover( *child1, *child2 );
		add_entity( child1 );
		// rosetta++ only produced one child upon recombination
		//if ( generations_[current_generation_].size() < size ) add_entity( child2 );
	}
}

///@brief add entities that are mutants of fit parents
void
GeneticAlgorithmBase::fill_by_mutation( core::Size size /* = 0 */ )
{
	runtime_assert(parent_entities_.size());

	if ( size == 0 ) size = max_population_size_;

	while ( generations_[current_generation_].size() < size ) {
		EntityOP child = tournament_select( parent_entities_ ).clone();
		entity_randomizer_->mutate( *child );
		add_entity( child );
	}
}


/// @brief progress to the next generation and generate new entities
/// @details
/// This method performs the following steps:
///
/// 1. If parent entities were not already specified, sets the parents to the current generation
/// 2. Increments the generation counter
/// 3. Copies the best scoring entities from the previous generation to the new one
/// 4. Generates new entities by crossover and/or mutation
/// 5. Clears the parent entities now that they have been used
void
GeneticAlgorithmBase::evolve_next_generation()
{
	// if parents were not already added, use the current generation as parents
	if (parent_entities_.size() == 0) add_parents_from_current_generation();
	// increment the generation counter
	++current_generation_;
	// copy the best scoring entities to the next generation
	propagate_best_from_previous_generation(number_to_propagate_);
	// add entities via crossover and mutation
	core::Size pop_size = generations_[current_generation_].size();
	fill_by_crossover( static_cast<core::Size>( fraction_by_recombination_*(max_population_size_-pop_size)+pop_size ) );
	fill_by_mutation( max_population_size_ );
	// reset parent entities now that they've been used
	parent_entities_.clear();
}

bool
GeneticAlgorithmBase::current_generation_complete()
{
	for (core::Size i = 1; i <= generations_[current_generation_].size(); ++i) {
		if (!generations_[current_generation_][i]->fitness_valid()) return false;
	}

	return true;
}

bool
GeneticAlgorithmBase::complete()
{
	if (current_generation_ < max_generations_) return false;
	if (current_generation_ == max_generations_ && !current_generation_complete()) return false;
	return true;
}


///@brief returns variable number of best (const) entities via vector of pointers to them
Entity::CAPs
GeneticAlgorithmBase::best_entities( core::Size num ) // nonconst method to permit sort
{
	std::sort( generations_[current_generation_].begin(), generations_[current_generation_].end(), lt_OP_deref< Entity > );
	EntityCAPs best_entities;
	pop_const_iter entity( generations_[current_generation_].begin() );
	while ( best_entities.size() < num && entity != generations_[current_generation_].end() ) {
		best_entities.push_back( EntityCAP(*entity) );
		++entity;
	}
	return best_entities;
}


///@brief pick two random entities from an unordered vector, return the one whose fitness is better
Entity const &
GeneticAlgorithmBase::tournament_select(
	utility::vector1< EntityCOP > const & pvec
) const
{
	using numeric::random::uniform;
  Entity const & e1( *pvec[ static_cast<Size>( uniform() * pvec.size() ) + 1 ] );
  Entity const & e2( *pvec[ static_cast<Size>( uniform() * pvec.size() ) + 1 ] );
  return ( e1.fitness() < e2.fitness() ? e1 : e2 );
}

GeneticAlgorithmBase::TraitEntityHashMap &
GeneticAlgorithmBase::entity_cache() { return entity_cache_; }

GeneticAlgorithmBase::TraitEntityHashMap const &
GeneticAlgorithmBase::entity_cache() const { return entity_cache_; }

utility::vector1<utility::vector1< EntityOP > > const &
GeneticAlgorithmBase::generations() const { return generations_; }

///@brief true const (read-only) access to entity population
Entity::COPs
GeneticAlgorithmBase::population( core::Size gen_num ) const
{
	return generations_[gen_num];
}

void
GeneticAlgorithmBase::print_generation_statistics(
	std::ostream & os,
	core::Size gen_num
) const
{
	std::set<EntityOP> earlier_generations_set;
	std::set<EntityOP> previous_generation_set;
	std::set<EntityOP> current_generation_set;

	// make a set of entities in the previous generation
	if (gen_num-1 >= 1) {
		previous_generation_set.insert(generations_[gen_num-1].begin(), generations_[gen_num-1].end());
	}

	// make a set of entities in earlier generations but not in the previous one
	for (core::Size i = 1; i+2 <= gen_num; ++i) {
		for (pop_const_iter iter(generations_[i].begin()), iter_end(generations_[i].end()); iter != iter_end;
		     ++iter) {
			if (previous_generation_set.find(*iter) == previous_generation_set.end()) {
				earlier_generations_set.insert(*iter);
			}
		}
	}

	core::Size resurrected_entities(0);
	core::Size heldover_entities(0);
	core::Size duplicate_new_entities(0);
	core::Size new_entities(0);
	EntityOP best_new_entity(NULL);

	for (pop_const_iter iter(generations_[gen_num].begin()), iter_end(generations_[gen_num].end()); iter != iter_end;
		     ++iter) {
		if (previous_generation_set.find(*iter) != previous_generation_set.end()) {
			++heldover_entities;
		} else if (earlier_generations_set.find(*iter) != earlier_generations_set.end()) {
			++resurrected_entities;
		} else if (current_generation_set.find(*iter) != current_generation_set.end()) {
			++duplicate_new_entities;
		} else {
			current_generation_set.insert(*iter);
			++new_entities;
			if (!best_new_entity || ((*iter)->fitness_valid() && (*iter)->fitness() < best_new_entity->fitness())) {
				best_new_entity = *iter;
			}
		}
	}

	utility::vector1<EntityOP> sorted_generation(generations_[gen_num]);
	std::sort( sorted_generation.begin(), sorted_generation.end(), lt_OP_deref< Entity > );
	core::Real gen_size(static_cast<core::Real>(sorted_generation.size()));

	os << "Distinct new entities: " << new_entities << std::endl;
	os << "Duplicate new entities: " << duplicate_new_entities << std::endl;
	os << "Entities from previous generation: " << heldover_entities << std::endl;
	os << "Entities resurrected from earlier generations: " << resurrected_entities << std::endl;
	os << "Fitness Percentiles: 0%=" << sorted_generation.front()->fitness()
	   << " 25%=" << sorted_generation[static_cast<core::Size>(ceil(gen_size*.25))]->fitness()
	   << " 50%=" << sorted_generation[static_cast<core::Size>(ceil(gen_size*.50))]->fitness()
	   << " 75%=" << sorted_generation[static_cast<core::Size>(ceil(gen_size*.75))]->fitness()
	   << " 100%=" << sorted_generation.back()->fitness() << std::endl;
	os << "Best new entity:" << '\n';
	os << *best_new_entity << std::endl;
}

void
GeneticAlgorithmBase::print_population( std::ostream & os ) const
{
	for ( pop_const_iter it( generations_[current_generation_].begin() ), end( generations_[current_generation_].end() );
				it != end; ++it ) {
		os << **it << '\n';
	}
	os << std::flush;
}


void
GeneticAlgorithmBase::print_cache( std::ostream & os ) const
{
	for ( TraitEntityHashMap::const_iterator it( entity_cache_.begin() ), end( entity_cache_.end() );
				it != end; ++it ) {
		if (it->second->fitness_valid()) {
			it->second->write_checkpoint(os);
			os << '\n';
		}
	}
	os << std::flush;
}

std::string
GeneticAlgorithmBase::entities_checkpoint_filename(std::string suffix /* = "" */) const
{
	if ( checkpoint_prefix_ == "" ) return "";
	std::string filename(checkpoint_prefix_ + ".ga.entities");
	filename += suffix;
	if (checkpoint_gzip_) filename += ".gz";
	return filename;
}

bool
GeneticAlgorithmBase::read_entities_checkpoint( bool overwrite /* = false */ )
{
	// if cache is not empty, then loading from checkpoint file is assumed not necessary by default
	if ( !entity_cache_.empty() && !overwrite ) return false;
	if ( checkpoint_prefix_ == "" ) return false;
	std::string const filename(entities_checkpoint_filename());
	utility::io::izstream file;
	file.open( filename.c_str() );
	if ( !file ) return false;
	TR(basic::t_info) << "Reading cached entities from file " << filename << '\n';

	std::string line;
	core::Size counter(0);
	EntityOP entity = new_entity();
	while ( entity->read_checkpoint(file) ) {
		TR(basic::t_debug) << *entity << '\n';
		entity_cache_[ entity->traits() ] = entity;
		++counter;
		entity = new_entity();
	}
	TR(basic::t_debug) << std::flush;
	file.close();

	TR(basic::t_info) << "Read " << counter << " cached fitnesses" << std::endl;
	return true;
}

bool
GeneticAlgorithmBase::write_entities_checkpoint() const
{
	if ( checkpoint_prefix_ == "" ) return false;
	std::string const filename(entities_checkpoint_filename());
	std::string const filename_tmp(entities_checkpoint_filename(".tmp"));
	utility::io::ozstream file( filename_tmp.c_str() );
	if ( !file ) {
		std::cerr << "trouble opening file " << filename << " for writing" << '\n';
		return false;
	}

	print_cache( file );

	file.close();

	// atomically replace the previous checkpoint file
	if ( std::rename(filename_tmp.c_str(), filename.c_str()) ) {
		std::cerr << "trouble renaming file " << filename_tmp << " to " << filename << std::endl;
		return false;
	}

	return true;
}

std::string
GeneticAlgorithmBase::generations_checkpoint_filename(std::string suffix /* = "" */) const
{
	if ( checkpoint_prefix_ == "" ) return "";
	std::string filename(checkpoint_prefix_ + ".ga.generations");
	filename += suffix;
	if (checkpoint_gzip_) filename += ".gz";
	return filename;
}


/// This seems to duplicate the functionality of the Entity's write_checkpoint function...
bool
GeneticAlgorithmBase::write_generations_checkpoint() const
{
	if ( checkpoint_prefix_ == "" ) return false;
	std::string const filename(generations_checkpoint_filename());
	std::string const filename_tmp(generations_checkpoint_filename(".tmp"));
	utility::io::ozstream file( filename_tmp.c_str() );
	if ( !file ) {
		std::cerr << "trouble opening file " << filename << " for writing" << std::endl;
		return false;
	}

	for (core::Size i = 1; i <= generations_.size(); ++i) {
		utility::vector1< EntityOP > const & generation(generations_[i]);
		if (generation.size()) {
			file << "generation " << i << '\n';
			for (core::Size j = 1; j <= generation.size(); ++j) {
				EntityElements const & traits(generation[j]->traits());
				for (core::Size k = 1; k <= traits.size(); ++k) {
					if (k != 1) file << ' ';
					file << traits[k]->to_string();
				}
				file << '\n';
			}
		}
	}

	file.close();

	// atomically replace the previous checkpoint file
	if ( std::rename(filename_tmp.c_str(), filename.c_str()) ) {
		std::cerr << "trouble renaming file " << filename_tmp << " to " << filename << std::endl;
		return false;
	}

	return true;
}

bool
GeneticAlgorithmBase::read_generations_checkpoint()
{
	if ( checkpoint_prefix_ == "" ) return false;
	std::string const filename(generations_checkpoint_filename());
	utility::io::izstream file( filename.c_str() );
	if ( !file ) {
		std::cerr << "trouble opening file " << filename << " for reading" << std::endl;
		return false;
	}

	core::Size gen_num = 0;
	std::string line;
	while (file.getline(line)) {
		std::istringstream iss(line);
		std::string word;
		if (!(iss >> word)) return false;
		if ( word == "generation" ) {
			if (!(iss >> gen_num)) return false;
			runtime_assert(gen_num > 0);
			runtime_assert(gen_num <= max_generations_);
			if (generations_.size() < gen_num) generations_.resize(gen_num);
			current_generation_ = gen_num;
		} else {
			// make sure a generation number has been set
			runtime_assert(gen_num > 0);
			EntityElements traits;
			traits.push_back( EntityElementFactory::get_instance()->element_from_string(word));
			while (iss >> word) {
				traits.push_back( EntityElementFactory::get_instance()->element_from_string(word));
			}
			add_entity(traits);
		}
	}

	file.close();

	return false;
}

bool
GeneticAlgorithmBase::read_checkpoint()
{
	// entities should be read first
	if ( !read_entities_checkpoint() ) return false;
	// only read generations if the entities were read successfully
	return (read_generations_checkpoint());
}

void
GeneticAlgorithmBase::rename_checkpoint_files() const
{
	if ( checkpoint_prefix_ == "" ) return;

	std::string suffix(".old");

	if ( utility::file::file_exists( entities_checkpoint_filename() ) ) {
		std::rename( entities_checkpoint_filename().c_str(), entities_checkpoint_filename(suffix).c_str() );
	}

	if ( utility::file::file_exists( generations_checkpoint_filename() ) ) {
		std::rename( generations_checkpoint_filename().c_str(), generations_checkpoint_filename(suffix).c_str() );
	}
}

EntityCOP GeneticAlgorithmBase::entity_template() const { return entity_template_; }
void GeneticAlgorithmBase::set_entity_template( EntityCOP entity ) { entity_template_ = entity; }


Entity::OP
GeneticAlgorithmBase::new_entity()
{
	if (entity_template_) {
		return entity_template_->clone();
	} else {
		return new Entity;
	}
}

GeneticAlgorithmBase::pop_iter       GeneticAlgorithmBase::current_generation_begin()
{ return generations_[current_generation_].begin(); }

GeneticAlgorithmBase::pop_iter       GeneticAlgorithmBase::current_generation_end()
{ return generations_[current_generation_].end(); }

GeneticAlgorithmBase::pop_const_iter GeneticAlgorithmBase::current_generation_begin() const
{ return generations_[current_generation_].begin(); }

GeneticAlgorithmBase::pop_const_iter GeneticAlgorithmBase::current_generation_end() const
{ return generations_[current_generation_].end(); }

GeneticAlgorithmBase::EntityOP       GeneticAlgorithmBase::curr_gen_entity( core::Size index )
{ return generations_[current_generation_][ index ]; }

core::Size GeneticAlgorithmBase::checkpoint_write_interval() const
{ return checkpoint_write_interval_; }


GeneticAlgorithm::GeneticAlgorithm() :
	GeneticAlgorithmBase(),
	fitness_function_(0)
{}

GeneticAlgorithm::~GeneticAlgorithm() {}

void
GeneticAlgorithm::set_func( FitnessFunctionOP f ) { fitness_function_ = f; }

void
GeneticAlgorithm::evaluate_fitnesses()
{
	write_generations_checkpoint();
	core::Size num_uncached(0);
	for ( pop_iter it( current_generation_begin() ), end( current_generation_end() );
	      it != end; ++it ) {
		if ( ! (*it)->fitness_valid() ) {
			//TR << "Evaluating fitness for entity " << (*it).get() << std::endl;
			//if ( (*it) != entity_cache_[ (*it)->traits() ] ) {
			//	TR << "WEIRD: Evaluating fitness for entity that is not in the cache" << std::endl;
			//}
			// entity's fitness not valid
			fitness_function_->evaluate( **it );
			// allow intermediate checkpointing for very large populations
			++num_uncached;
			if ( checkpoint_write_interval() && num_uncached % checkpoint_write_interval() == 0 ) {
				write_entities_checkpoint();
			}
		}
	}
	write_entities_checkpoint();
}



} // namespace genetic_algorithm
} // namespace protocols

