// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/andrew/apl_msd.cc
/// @brief  Multistate design executable as written by apl.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Pack Daemon headers
#include <devel/pack_daemon/EntityCorrespondence.hh>
#include <devel/pack_daemon/MultistateAggregateFunction.hh>
#include <devel/pack_daemon/MultistateFitnessFunction.hh>
#include <devel/pack_daemon/PackDaemon.hh>

/// Core headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/multistate_design/util.hh>

// Utility headers
#include <utility/string_util.hh>

// option key includes
#include <basic/options/keys/ms.OptionKeys.gen.hh>

using basic::t_info;
static THREAD_LOCAL basic::Tracer TR( "app.andrew.apl_msd", t_info );

class SimpleDGBindAggregateFunction : public devel::pack_daemon::MultistateAggregateFunction
{
public:
	typedef devel::pack_daemon::MultistateAggregateFunction parent;

public:
	SimpleDGBindAggregateFunction() : parent() {}
	virtual ~SimpleDGBindAggregateFunction() {}

	virtual core::Real evaluate( utility::vector1< core::Real > const & vec, Entity const & ) {
		if ( vec.size() != 2 ) {
			utility_exit_with_message( "SimpleDGBindAggregateFunction expects exactly 2 states" );
		}
		return vec[ 1 ] - vec[ 2 ];
	}

	virtual StateIndices select_relevant_states( StateEnergies const &, Entity const & ) {
		StateIndices two_vec;
		two_vec.push_back( 1 );
		two_vec.push_back( 2 );
		return two_vec;
	}
};

int main( int argc, char ** argv )
{
	try {
	using namespace devel::pack_daemon;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::genetic_algorithm;

	devel::init( argc, argv );

	std::string entity_resfile( "2wo2_entity.resfile" );
	std::string secondary_resfile( "2wo2_secondary.resfile" );
	std::string bound_pdb( "2wo2.pdb" );
	std::string unbound_pdb( "2wo2_sep.pdb" );
	std::string entity_correspondence_file( "entity_map.txt" );

	DaemonSetOP ds = new DaemonSet;
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	ds->set_score_function( *sfxn );

	ds->set_entity_resfile( entity_resfile );

	ds->add_pack_daemon( 1, bound_pdb,   entity_correspondence_file, secondary_resfile );
	ds->add_pack_daemon( 2, unbound_pdb, entity_correspondence_file, secondary_resfile );

	ds->setup_daemons();

	devel::pack_daemon::MultistateFitnessFunctionOP func = new devel::pack_daemon::MultistateFitnessFunction;
	func->daemon_set( ds );
	func->aggregate_function( new SimpleDGBindAggregateFunction );
	func->set_history_size( 3 );

	// <stolen code>
	// set up genetic algorithm
	GeneticAlgorithm ga;
	ga.set_max_generations( option[ ms::generations ] );
	ga.set_max_pop_size( option[ ms::pop_size ]() );
	ga.set_num_to_propagate( static_cast<core::Size>(0.5*option[ ms::pop_size ]) );
	ga.set_frac_by_recomb( option[ ms::fraction_by_recombination ]() );

	// set up sequence randomizer
	PositionSpecificRandomizer::OP rand = new PositionSpecificRandomizer;
	rand->set_mutation_rate( option[ ms::mutate_rate ]() );
	for ( Size ii = 1; ii <= ds->entity_task()->size(); ++ii ) {
		EntityElements ii_elements =
			protocols::multistate_design::list_amino_acid_options(
			ii, ds->entity_task()->residue_task( ii ) );
		rand->append_choices( ii_elements );
	}
	// done setting up randomizer
	ga.set_rand( rand );
	ga.set_func( func );
	// </stolen code>

	/// Now initialize the GA with a population from which to begin exploration.
	/// Initialize this population completely randomly so that the native sequence
	/// is not arrived at unfairly.

	ga.fill_with_random_entities();
	// clear parents for the next generation
	// ga.clear_parents(); // do I need this?
	// loop over generations
	while ( !ga.complete() ) {
		if (ga.current_generation_complete()) ga.evolve_next_generation();
		TR(t_info) << "Generation " << ga.current_generation() << ":" << std::endl;
		ga.evaluate_fitnesses();
		ga.print_population( TR(t_info) );
	}

	// <stolen code>
	// sort local copy of sequence/fitness cache
	typedef GeneticAlgorithm::TraitEntityHashMap TraitEntityHashMap;
	TraitEntityHashMap const & cache( ga.entity_cache() );
	utility::vector1<Entity::OP> sortable;
	for ( TraitEntityHashMap::const_iterator it( cache.begin() ), end( cache.end() ); it != end; ++it ) {
		if ( it->second == 0 ) continue; /// Why would this be?
		sortable.push_back( it->second );
	}
	std::sort( sortable.begin(), sortable.end(), lt_OP_deref< Entity > );

	TR(t_info) << "Evaluated " << sortable.size() << " sequences out of " << rand->library_size()
	           << " possible.\nBest sequences:\n";

	int counter = 0;
	for ( utility::vector1<Entity::OP>::const_iterator it( sortable.begin() ), end( sortable.end() );
				it != end && counter < option[ ms::numresults ](); ++it, ++counter ) {
		// This next line of code is in fact original and not stolen from JA.
		std::list< std::pair< core::Size, core::pose::PoseOP > > pose_list  = func->recover_relevant_poses_for_entity( **it );
		utility::vector1< core::pose::PoseOP > poses( 2 );
		for ( std::list< std::pair< core::Size, core::pose::PoseOP > >::const_iterator
				iter = pose_list.begin(), iter_end = pose_list.end(); iter != iter_end; ++iter ) {
			poses[ iter->first ] = iter->second;
		}
		poses[ 1 ]->dump_pdb( "msa_design_" + utility::to_string( counter ) + "_bound.pdb" );
		poses[ 2 ]->dump_pdb( "msa_design_" + utility::to_string( counter ) + "_unbound.pdb" );
		break;
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
