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

#ifdef USEMPI
#include <mpi.h>
#endif

/// Pack Daemon headers
#include <devel/pack_daemon/EntityCorrespondence.hh>
#include <devel/pack_daemon/DynamicAggregateFunction.hh>
#include <devel/pack_daemon/MultistateAggregateFunction.hh>
#include <devel/pack_daemon/MultistateFitnessFunction.hh>
#include <devel/pack_daemon/PackDaemon.hh>

/// Core headers
#include <devel/init.hh>
//#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>


#include <basic/Tracer.hh>

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/multistate_design/util.hh>

// Utility headers
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/numeric.functions.hh>

// option key includes
#include <basic/options/keys/ms.OptionKeys.gen.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


OPT_1GRP_KEY( String, msd, entity_resfile )
OPT_1GRP_KEY( String, msd, fitness_file )
OPT_1GRP_KEY( Integer, msd, double_lazy_ig_mem_limit )
OPT_1GRP_KEY( Boolean, msd, dont_score_bbhbonds )
OPT_1GRP_KEY( Boolean, msd, exclude_background_energies )

/*  Option( 'double_lazy_ig_mem_limit', 'Integer',
desc="The amount of memory, in MB, that each double-lazy interaction graph should be allowed \
to allocate toward rotamer pair energies.  Using this flag will not trigger the \
use of the double-lazy interaction graph, and this flag is not read in the PackerTask's \
initialize_from_command_line routine.  For use in multistate design",
default='0',
),*/


using basic::t_info;
using basic::t_debug;
static basic::Tracer TR( "app.andrew.mpi_msd", t_info );

/*class SimpleDGBindAggregateFunction : public devel::pack_daemon::MultistateAggregateFunction
{
public:
typedef devel::pack_daemon::MultistateAggregateFunction parent;

public:
SimpleDGBindAggregateFunction() : parent() {}
virtual ~SimpleDGBindAggregateFunction() {}

virtual core::Real evaluate( utility::vector1< core::Real > const & vec, Entity const &  ) {
if ( vec.size() != 2 ) {
utility_exit_with_message( "SimpleDGBindAggregateFunction expects exactly 2 states" );
}
return vec[ 1 ] - vec[ 2 ];
}

virtual StateIndices select_relevant_states( StateEnergies const & ) {
StateIndices two_vec;
two_vec.push_back( 1 );
two_vec.push_back( 2 );
return two_vec;
}
};*/

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief different set of choices at each position in Entity's traits
class Mutate1Randomizer : public protocols::genetic_algorithm::PositionSpecificRandomizer {
public:
	typedef protocols::genetic_algorithm::PositionSpecificRandomizer parent;

public:

	virtual ~Mutate1Randomizer() {}

	virtual void mutate( protocols::genetic_algorithm::Entity & entity )
	{
		if ( mutation_rate() == 1.0 ) {
			for ( Size ii = 1; ii <= entity.traits().size(); ++ii ) {
				core::Size const n_mutation_choices( choices()[ ii ].size() );
				Size new_element_ind = static_cast< core::Size >( numeric::random::uniform() * n_mutation_choices ) + 1;
				entity.set_entity_element( ii, choices()[ ii ][ new_element_ind ] );
			}
		} else {
			Size pos_to_mutate = static_cast< Size > ( entity.traits().size() * numeric::random::uniform() ) + 1;
			core::Size const n_mutation_choices( choices()[ pos_to_mutate ].size() );
			Size new_element_ind = static_cast< core::Size >( numeric::random::uniform() * n_mutation_choices ) + 1;
			entity.set_entity_element( pos_to_mutate, choices()[ pos_to_mutate ][ new_element_ind ] );
		}
	}

};


int main( int argc, char ** argv )
{
	try {

		using namespace utility;
		using namespace devel::pack_daemon;
		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::genetic_algorithm;

		NEW_OPT( msd::entity_resfile, "Resfile for the entity elements which are shared between the multiple states", "" );
		NEW_OPT( msd::fitness_file,   "Fitness function file specifying the multiple states and the objective function to optimize", "" );
		NEW_OPT( msd::double_lazy_ig_mem_limit, "The amount of memory, in MB, that each double-lazy interaction graph should be allowed to allocate toward rotamer pair energies.", 0 );
		NEW_OPT( msd::dont_score_bbhbonds, "Disable the default activation of the decompose_bb_hb_into_pair_energies flag for hbonds", false );
		NEW_OPT( msd::exclude_background_energies, "Disable the default activation of the inclusion of background one-body and background/background two-body interaction energies in the state energies (which until now held only the packer energies)", false );

		devel::init( argc, argv );

		//std::string secondary_resfile( "2wo2_secondary.resfile" );
		//std::string bound_pdb( "2wo2.pdb" );
		//std::string unbound_pdb( "2wo2_sep.pdb" );
		//std::string entity_correspondence_file( "entity_map.txt" );
		//std::string daf_filename = "objective_function.daf";

		if ( ! option[ msd::entity_resfile ].user() ) {
			utility_exit_with_message("The entity resfile must be specified for the mpi_msd application with the -msd::entity_resfile <filename> flag" );
		}

		if ( ! option[ msd::fitness_file ].user() ) {
			utility_exit_with_message("The fitness-function file must be specified for the mpi_msd application with the -msd::fitness_file <filename> flag" );
		}

		std::string entity_resfile( option[ msd::entity_resfile ] );
		std::string daf_filename( option[ msd::fitness_file ] );

		DaemonSetOP ds = new DaemonSet;
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

		if ( ! option[ msd::dont_score_bbhbonds ] ) {
			/// Count bb/bb hydrogen bonds in the packer energy; otherwise, the MSD
			/// code cannot say that one set of bb/bb contacts is worse than
			/// another set of bb/bb contacts
			if ( mpi_rank() == 0 ) { TR << "Activating decompose_bb_hb_into_pair_energies in the score function" << std::endl; }
			using namespace core;
			scoring::methods::EnergyMethodOptionsOP emopts( new scoring::methods::EnergyMethodOptions( sfxn->energy_method_options() ));
			emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
			sfxn->set_energy_method_options( *emopts );
		}
		if ( option[ msd::exclude_background_energies ] ) {
			ds->set_include_background_energies( false );
		}


		ds->set_score_function( *sfxn );
		ds->set_entity_resfile( entity_resfile );
		if ( option[ msd::double_lazy_ig_mem_limit ].user() ) {
			std::cout << "Setting dlig nmeg limit" << option[ msd::double_lazy_ig_mem_limit ] << std::endl;
			ds->set_dlig_nmeg_limit( option[ msd::double_lazy_ig_mem_limit ] );
		}

		//if ( mpi_rank() == 0 ) {
		// ds->add_pack_daemon( 1, bound_pdb,   entity_correspondence_file, secondary_resfile );
		//} else {
		// ds->add_pack_daemon( 2, unbound_pdb, entity_correspondence_file, secondary_resfile );
		//}


		if ( mpi_rank() == 0 ) {
			devel::pack_daemon::MPIMultistateFitnessFunctionOP func = new devel::pack_daemon::MPIMultistateFitnessFunction;
			devel::pack_daemon::DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
			daf->set_num_entity_elements( ds->entity_task()->size() );
			utility::io::izstream daf_file( daf_filename );
			try {
				daf->initialize_from_input_file( ds, daf_file );
			} catch ( utility::excn::EXCN_Msg_Exception & e ) {
				std::cerr << "Caught exception" << std::endl;
				std::cerr << e.msg() << std::endl;
				exit(1);
			}
			func->daemon_set( ds );
			func->set_num_pack_daemons( daf->num_states() );
			func->aggregate_function( daf );
			func->set_history_size( option[ ms::numresults ]() );

			// <stolen code>
			// set up genetic algorithm
			GeneticAlgorithm ga;
			ga.set_max_generations( option[ ms::generations ] );
			ga.set_max_pop_size( option[ ms::pop_size ]() );
			ga.set_num_to_propagate( static_cast<core::Size>(0.5*option[ ms::pop_size ]) );
			ga.set_frac_by_recomb( option[ ms::fraction_by_recombination ]() );

			// set up sequence randomizer
			Mutate1Randomizer::OP rand = new Mutate1Randomizer;
			// reset the default value of 1.0, but the mutation-rate variable is not used by the Mutate1Randomizer!
			rand->set_mutation_rate( 0.0 /*option[ ms::mutate_rate ]()*/ );

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
				clock_t starttime = clock();
				if ( ga.current_generation_complete() ) ga.evolve_next_generation();
				ga.evaluate_fitnesses();
				if ( TR.visible( t_debug ) ) {
					TR(t_debug) << "Generation " << ga.current_generation() << ":" << std::endl;
					ga.print_population( TR(t_debug) );
				}
				clock_t stoptime = clock();
				TR << "Generation " << ga.current_generation() << " took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC
					<< " seconds; best fitness = " << ga.best_fitness_from_current_generation() << std::endl;
			}

			TR(t_info) << "Final Generation " << ga.current_generation() << ":" << std::endl;
			ga.print_population( TR(t_info) );


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
				<< " possible." << std::endl;

			int counter = 1;
			for ( utility::vector1<Entity::OP>::const_iterator it( sortable.begin() ), end( sortable.end() );
					it != end && counter <= option[ ms::numresults ](); ++it, ++counter ) {
				if ( ! *it ) {
					--counter;
					continue;
				}
				TR << "Top set #" << counter << ". Writing pdbs for entity: " << **it << std::endl;
				// This next line of code is in fact original and not stolen from JA.
				typedef std::list< std::pair< core::Size, core::pose::PoseOP > > SizePosePairList;
				SizePosePairList pose_list  = func->recover_relevant_poses_for_entity( **it );
				/* Old code from my hard-coded test case that involved exactly four poses.
				utility::vector1< core::pose::PoseOP > poses( 4 );
				for ( std::list< std::pair< core::Size, core::pose::PoseOP > >::const_iterator
				iter = pose_list.begin(), iter_end = pose_list.end(); iter != iter_end; ++iter ) {
				poses[ iter->first ] = iter->second;
				}
				poses[ 1 ]->dump_pdb( "msa_design_" + utility::to_string( counter ) + "_bound1.pdb" );
				poses[ 2 ]->dump_pdb( "msa_design_" + utility::to_string( counter ) + "_bound2.pdb" );
				poses[ 3 ]->dump_pdb( "msa_design_" + utility::to_string( counter ) + "_unbound1.pdb" );
				poses[ 4 ]->dump_pdb( "msa_design_" + utility::to_string( counter ) + "_unbound2.pdb" );
				//break;
				*/
				for ( SizePosePairList::const_iterator iter = pose_list.begin(), iter_end = pose_list.end();
						iter != iter_end; ++iter ) {
					std::string output_name = "msd_output_";
					output_name += utility::to_string( counter );
					output_name += "_" + daf->state_name( iter->first ) + ".pdb";
					TR << "Writing structure " << output_name << " with score: " << (*sfxn)( *(iter->second) ) << std::endl;
					utility::io::ozstream outfile( output_name );
					core::io::pdb::dump_pdb( *(iter->second), outfile );
					core::io::pdb::extract_scores( *(iter->second), outfile );
				}
			}

			func->send_spin_down_signal();
		} else {
			ds->activate_daemon_mode();
		}

#ifdef USEMPI
	MPI_Finalize();
#endif

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


