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

#ifndef INCLUDED_protocols_genetic_algorithm_GeneticAlgorithm_hh
#define INCLUDED_protocols_genetic_algorithm_GeneticAlgorithm_hh

// Unit headers
#include <protocols/genetic_algorithm/GeneticAlgorithm.fwd.hh>

#include <protocols/genetic_algorithm/Entity.fwd.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/genetic_algorithm/FitnessFunction.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// AUTO-REMOVED #include <utility/file/file_sys_util.hh> // file_exists
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

// AUTO-REMOVED #include <numeric/random/random.fwd.hh>

#include <boost/unordered_map.hpp>

#include <algorithm> // std::copy

//Auto Headers
namespace protocols {
namespace genetic_algorithm {

class GeneticAlgorithm : public utility::pointer::ReferenceCount {

public:
	typedef FitnessFunction::OP FitnessFunctionOP;
	typedef EntityRandomizer::OP EntityRandomizerOP;
	typedef Entity::OP EntityOP;
	typedef Entity::COP EntityCOP;
	typedef Entity::COPs EntityCOPs;
	typedef Entity::CAP EntityCAP;
	typedef Entity::CAPs EntityCAPs;
	typedef utility::vector1< EntityOP >::iterator pop_iter;
	typedef utility::vector1< EntityOP >::const_iterator pop_const_iter;
	typedef utility::vector1< EntityCOP >::const_iterator pop_const_const_iter;
	typedef boost::unordered_map< EntityElements , EntityOP, Vec1Hash, EntityElementsEqual > TraitEntityHashMap;

public:
	GeneticAlgorithm();
	virtual ~GeneticAlgorithm();

	virtual EntityOP add_entity( EntityElements const & traits );
	virtual EntityOP add_entity( EntityOP entity );
	virtual EntityOP add_parent_entity( EntityElements const & traits );
	virtual EntityOP add_parent_entity( EntityOP entity );
	virtual void clear_parents();
	virtual void add_parents_from_current_generation();
	virtual void propagate_best_from_previous_generation( core::Size size = 1, bool unique = true );
	virtual void fill_with_random_entities( core::Size size = 0 );
	virtual void fill_with_perturbations_of_existing_entities( core::Size size = 0 );
	virtual void fill_by_crossover( core::Size size = 0 );
	virtual void fill_by_mutation( core::Size size = 0 );
	virtual void evaluate_fitnesses();
	virtual void evolve_next_generation();
	virtual bool current_generation_complete();
	virtual bool complete();
	virtual core::Real best_fitness_from_current_generation() const;

	virtual void set_func( FitnessFunctionOP f );
	virtual void set_rand( EntityRandomizerOP r );
	virtual core::Size current_generation() const { return current_generation_; }
	virtual core::Size max_generations() const { return max_generations_; }
	virtual void set_max_generations( core::Size s );
	virtual void set_max_pop_size( core::Size s ) { max_population_size_ = s; }
	virtual void set_num_to_propagate( core::Size s ) { number_to_propagate_ = s; }
	virtual void set_frac_by_recomb( core::Real f ) { fraction_by_recombination_ = f; }
	virtual void set_checkpoint_prefix( std::string const & p ) { checkpoint_prefix_ = p; }
	virtual void set_checkpoint_write_interval( core::Size i ) { checkpoint_write_interval_ = i; }
	virtual void set_checkpoint_gzip( bool b ) { checkpoint_gzip_ = b; }
	virtual void set_checkpoint_rename( bool b ) { checkpoint_rename_ = b; }

	///@brief non-const to permit sort
	virtual EntityCAPs best_entities( core::Size num );
	virtual Entity const & tournament_select( utility::vector1< EntityCOP > const & pvec ) const;
	virtual TraitEntityHashMap & entity_cache();
	virtual TraitEntityHashMap const & entity_cache() const;
	virtual utility::vector1<utility::vector1< EntityOP > > const & generations() const;
	///@brief true const (read-only) access to entity population: new vector of const pointers
	virtual EntityCOPs population( core::Size gen_num ) const;
	virtual void print_generation_statistics( std::ostream & os, core::Size gen_num ) const;
	virtual void print_population( std::ostream & ) const;
	virtual void print_cache( std::ostream & ) const;
	virtual std::string entities_checkpoint_filename(std::string suffix = "") const;
	///@brief for checkpointing fitness cache
	virtual bool read_entities_checkpoint( bool overwrite = false );
	///@brief for checkpointing fitness cache
	virtual bool write_entities_checkpoint() const;
	virtual std::string generations_checkpoint_filename(std::string suffix = "") const;
	virtual bool write_generations_checkpoint() const;
	virtual bool read_generations_checkpoint();
	virtual bool read_checkpoint();
	///@brief allows the prevention of accidental reuse of checkpoint files
	virtual void rename_checkpoint_files() const;
	virtual EntityCOP entity_template() const;
	virtual void set_entity_template(EntityCOP entity);
	virtual EntityOP new_entity();

protected:

private:
	utility::vector1< utility::vector1< EntityOP > > generations_;
	utility::vector1< EntityCOP > parent_entities_;
	TraitEntityHashMap entity_cache_;
	FitnessFunctionOP fitness_function_;
	EntityRandomizerOP entity_randomizer_;
	EntityCOP entity_template_;
	core::Size current_generation_;
	core::Size max_generations_;
	core::Size max_population_size_;
	core::Size number_to_propagate_;
	core::Real fraction_by_recombination_;
	std::string checkpoint_prefix_;
	core::Size checkpoint_write_interval_;
	bool checkpoint_gzip_;
	bool checkpoint_rename_;

};


} // namespace genetic_algorithm
} // namespace protocols

#endif
