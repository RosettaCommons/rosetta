// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/andrew/apl_msd.cc
/// @brief  Multistate design executable as written by apl.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

/// Pack Daemon headers
#include <protocols/pack_daemon/EntityCorrespondence.hh>
#include <protocols/pack_daemon/DynamicAggregateFunction.hh>
#include <protocols/pack_daemon/MultistateAggregateFunction.hh>
#include <protocols/pack_daemon/MultistateFitnessFunction.hh>
#include <protocols/pack_daemon/PackDaemon.hh>
#include <protocols/pack_daemon/util.hh>

/// Core headers
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/MetricValue.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>

#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/genetic_algorithm/Mutate1Randomizer.hh>
#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/multistate_design/util.hh>

//#include <protocols/toolbox/pose_metric_calculators/HPatchCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

// Devel headers
//This comment brackets a snippet of code that is automatically stripped from the release, usually for semi-forbidden interactions with the unreleased devel library
//If you want non-devel code stripped from the release, see the release machinery in tools/release and contact the release manager (Steven Lewis smlewi@gmail.com at this time)
///DONOTRELEASE_TOP
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
///DONOTRELEASE_BOTTOM
//end automatic stripping comment

// Utility headers
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>

// option key includes
#include <basic/options/keys/ms.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


OPT_1GRP_KEY( String, msd, entity_resfile )
OPT_1GRP_KEY( String, msd, fitness_file )
OPT_1GRP_KEY( StringVector, msd, seed_sequences )
OPT_1GRP_KEY( Boolean, msd, fill_gen1_from_seed_sequences )
OPT_1GRP_KEY( Integer, msd, double_lazy_ig_mem_limit )
OPT_1GRP_KEY( Boolean, msd, dont_score_bbhbonds )
OPT_1GRP_KEY( Boolean, msd, exclude_background_energies )
OPT_1GRP_KEY( String, msd, seed_sequence_from_input_pdb )
OPT_1GRP_KEY( String, msd, seed_sequence_using_correspondence_file )

/*		Option( 'double_lazy_ig_mem_limit', 'Integer',
			desc="The amount of memory, in MB, that each double-lazy interaction graph should be allowed \
				to allocate toward rotamer pair energies.  Using this flag will not trigger the \
				use of the double-lazy interaction graph, and this flag is not read in the PackerTask's \
				initialize_from_command_line routine.  For use in multistate design",
			default='0',
		),*/


using basic::t_info;
using basic::t_debug;
static basic::Tracer TR("apps.public.design.mpi_msd",t_info);

/*class SimpleDGBindAggregateFunction : public protocols::pack_daemon::MultistateAggregateFunction
{
public:
	typedef protocols::pack_daemon::MultistateAggregateFunction parent;

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


//protocols::genetic_algorithm::EntityElements
//entity_elements_from_1letterstring(
//	std::string const & input
//)
//{
//	protocols::genetic_algorithm::EntityElements elements( input.size() );
//	for ( core::Size ii = 0, count = 1; ii < input.size(); ++ii, ++count ) {
//		std::ostringstream output;
//		output << "AA:" << count << ":" << input[ ii ];
//		elements[ count ] = protocols::genetic_algorithm::EntityElementFactory::get_instance()->element_from_string( output.str() );
//	}
//	return elements;
//}

std::string
read_native_sequence_for_entity_elements( core::Size n_designed_positions )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace protocols::pack_daemon;

	if ( ! option[ msd::seed_sequence_using_correspondence_file ].user() ) {
		utility_exit_with_message( "Must provide correspondence file to read the native sequence" );
	}

	std::string pdb_name = option[ msd::seed_sequence_from_input_pdb ];
	std::string correspondence_file_name = option[ msd::seed_sequence_using_correspondence_file ];

	/// Read in the pdb
	pose::Pose pose;
	import_pose::pose_from_pdb( pose, pdb_name );

	utility::io::izstream correspondence_file( correspondence_file_name );
	if ( ! correspondence_file ) {
		utility_exit_with_message( "Could not open correspondence file named: " + correspondence_file_name );
	}

	EntityCorrespondenceOP ec = new EntityCorrespondence;
	ec->set_pose( new pose::Pose( pose ));
	ec->set_num_entities( n_designed_positions );
	ec->initialize_from_correspondence_file( correspondence_file );


	std::map< Size, chemical::AA > aa_for_design_position;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size ii_entity = ec->entity_for_residue( ii );
		if ( ii_entity  == 0 ) continue;
		if ( aa_for_design_position.find( ii_entity ) != aa_for_design_position.end() ) {
			utility_exit_with_message( "Repeat entity element for native pdb: " + pdb_name + " with correspondence file " +
				correspondence_file_name + ".  Entity correspondence file should only include each residue once");
		}
		aa_for_design_position[ ii_entity ] = pose.residue_type( ii ).aa();
	}
	std::string aa_string( n_designed_positions, 'X' );
	for ( Size ii = 1; ii <= n_designed_positions; ++ii ) {
		if ( aa_for_design_position.find( ii ) == aa_for_design_position.end() ) {
			utility_exit_with_message( "Did not find residue assigned to correspond to entity element " +
				utility::to_string( ii ) + " while reading correspondence file " + correspondence_file_name );
		}
		aa_string[ ii-1 ] = oneletter_code_from_aa( aa_for_design_position[ ii ] );
	}
	return aa_string;


}

///////////////////////////////////////////////////////////////////////////////
///////////////////////   HPatchNPDCalculator   ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class HPatchNPDCalculator : public protocols::pack_daemon::NPDPropCalculator
{
public:

	virtual
	core::Real
	calculate( core::pose::Pose const & p ) {
		return core::pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( p );
	}

};

class HPatchNPDCalculatorCreator : public protocols::pack_daemon::NPDPropCalculatorCreator
{
	virtual
	std::string
	calculator_name() const {return "hpatch"; }

	virtual
	protocols::pack_daemon::NPDPropCalculatorOP
	new_calculator() const { return new HPatchNPDCalculator; }
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////   HPatchByChainNPDCalculator   ////////////////////////
///////////////////////////////////////////////////////////////////////////////


class HPatchByChainNPDCalculator : public protocols::pack_daemon::NPDPropCalculator
{
public:
	virtual
	void
	setup(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & task
	){
		chains_ = pose.split_by_chain();
		task_ = task.clone();
		Size last_chain( 1 ), first_residue_for_chain( 1 );
		resid_2_chain_and_resid_.resize( pose.total_residue() );
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue( ii ).chain() != last_chain ) {
				last_chain = pose.residue( ii ).chain();
				first_residue_for_chain = ii;
			}
			resid_2_chain_and_resid_[ ii ] = std::make_pair( pose.residue( ii ).chain(), ii + 1 - first_residue_for_chain );
		}
	}

	virtual
	core::Real
	calculate( core::pose::Pose const & p ) {
		// MJO COMMENTING OUT BECAUSE IT IS UNUSED:
		// Size chain_offset = 0;
		for ( core::Size ii = 1; ii <= p.total_residue(); ++ii ) {
			if ( ! task_->being_packed( ii ) ) continue;
			chains_[ resid_2_chain_and_resid_[ ii ].first ]->replace_residue(
				resid_2_chain_and_resid_[ ii ].second, p.residue( ii ), false );
		}

		core::Real hpatch_sum = 0;
		for ( Size ii = 1; ii <= chains_.size(); ++ii ) {
			hpatch_sum += core::pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( *chains_[ ii ] );
		}
		return hpatch_sum;
	}
private:
	utility::vector1< core::pose::PoseOP > chains_;
	core::pack::task::PackerTaskOP task_;
	utility::vector1< std::pair< core::Size, core::Size > > resid_2_chain_and_resid_;
};

class HPatchByChainNPDCalculatorCreator : public protocols::pack_daemon::NPDPropCalculatorCreator
{
	virtual
	std::string
	calculator_name() const {return "hpatch_by_chain"; }

	virtual
	protocols::pack_daemon::NPDPropCalculatorOP
	new_calculator() const { return new HPatchByChainNPDCalculator; }
};



///////////////////////////////////////////////////////////////////////////////
///////////////////////      BunsCalculator     ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class NBuriedUnsatsCalcultor : public protocols::pack_daemon::NPDPropCalculator
{
public:
	virtual
	void
	setup(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & task
	){
		pose_ = new core::pose::Pose( pose );
		task_ = task.clone();
		sfxn_ = core::scoring::getScoreFunction();

		//register calculators

		//This comment brackets a snippet of code that is automatically stripped from the release, usually for semi-forbidden interactions with the unreleased devel library
		//If you want non-devel code stripped from the release, see the release machinery in tools/release and contact the release manager (Steven Lewis smlewi@gmail.com at this time)
		///DONOTRELEASE_TOP
		if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
			core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new devel::vardist_solaccess::VarSolDistSasaCalculator;
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
		}
		///DONOTRELEASE_BOTTOM
		//end automatic stripping comment

		if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "num_hbonds" ) ) {
			core::pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );
		}

		if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "unsat" ) ) {
			core::pose::metrics::PoseMetricCalculatorOP unsat_calculator = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds");
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );
		}

	}

	virtual
	core::Real
	calculate( core::pose::Pose const & p ) {
		for ( core::Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
			if ( ! task_->being_packed( ii ) ) continue;
			pose_->replace_residue( ii, p.residue( ii ), false );
		}
		(*sfxn_)(*pose_);
		basic::MetricValue< core::Size > nburied_unsats;
		pose_->metric( "unsat", "all_bur_unsat_polars", nburied_unsats );
		return nburied_unsats.value();
	}

private:
	core::pose::PoseOP             pose_;
	core::pack::task::PackerTaskOP task_;
	core::scoring::ScoreFunctionOP sfxn_;
};

class NBuriedUnsatsCalcultorCreator : public protocols::pack_daemon::NPDPropCalculatorCreator
{
	virtual
	std::string
	calculator_name() const {return "nbunsats"; }

	virtual
	protocols::pack_daemon::NPDPropCalculatorOP
	new_calculator() const { return new NBuriedUnsatsCalcultor; }
};





///////////////////////////////////////////////////////////////////////////////
///////////////////////           main          ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////



int main( int argc, char ** argv )
{
	try {

	using namespace utility;
	using namespace protocols::pack_daemon;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::genetic_algorithm;

	NEW_OPT( msd::entity_resfile, "Resfile for the entity elements which are shared between the multiple states", "" );
	NEW_OPT( msd::fitness_file,   "Fitness function file specifying the multiple states and the objective function to optimize", "" );
	NEW_OPT( msd::double_lazy_ig_mem_limit, "The amount of memory, in MB, that each double-lazy interaction graph should be allowed to allocate toward rotamer pair energies.", 0 );
	NEW_OPT( msd::dont_score_bbhbonds, "Disable the default activation of the decompose_bb_hb_into_pair_energies flag for hbonds", false );
	NEW_OPT( msd::exclude_background_energies, "Disable the default activation of the inclusion of background one-body and background/background two-body interaction energies in the state energies (which until now held only the packer energies)", false );
	NEW_OPT( msd::seed_sequences, "Seed the GA's population with the given input sequences", "" );
	NEW_OPT( msd::fill_gen1_from_seed_sequences, "Fill the entirety of the first generation from perturbations off the seed sequences that were provided", false );
	NEW_OPT( msd::seed_sequence_from_input_pdb, "Seed the GA's population with the given input sequences using the sequence already present in a given pdb file; requires the use of the msd::seed_sequence_using_correspondence_file flag ", "" );
	NEW_OPT( msd::seed_sequence_using_correspondence_file, "The name of the correspondence file to guide the seeding of the GA's population with the sequence from a particular pdb", "" );

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
	ds->add_npdpro_calculator_creator( new HPatchNPDCalculatorCreator );
	ds->add_npdpro_calculator_creator( new HPatchByChainNPDCalculatorCreator );
	ds->add_npdpro_calculator_creator( new NBuriedUnsatsCalcultorCreator );

	core::scoring::ScoreFunctionOP sfxn = core::scoring::getScoreFunction();

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
	//	ds->add_pack_daemon( 1, bound_pdb,   entity_correspondence_file, secondary_resfile );
	//} else {
	//	ds->add_pack_daemon( 2, unbound_pdb, entity_correspondence_file, secondary_resfile );
	//}


	if ( mpi_rank() == 0 ) {
		protocols::pack_daemon::MPIMultistateFitnessFunctionOP func = new protocols::pack_daemon::MPIMultistateFitnessFunction;
		protocols::pack_daemon::DynamicAggregateFunctionDriverOP daf = new DynamicAggregateFunctionDriver;
		daf->set_num_entity_elements( ds->entity_task()->total_residue() );
		daf->set_score_function( *sfxn ); // assume one score function for the entire
		utility::io::izstream daf_file( daf_filename );
		try {
			daf->initialize_from_input_file( ds, daf_file );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Caught exception" << std::endl;
			std::cerr << e.msg() << std::endl;
			exit(1);
		}
		func->daemon_set( ds );
		func->set_num_pack_daemons(   daf->num_states()         );
		func->set_num_npd_properties( daf->num_npd_properties() );
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
		protocols::genetic_algorithm::Mutate1RandomizerOP rand = new protocols::genetic_algorithm::Mutate1Randomizer;
		// reset the default value of 1.0, but the mutation-rate variable is not used by the Mutate1Randomizer!
		rand->set_mutation_rate( 0.0 /*option[ ms::mutate_rate ]()*/ );
		protocols::pack_daemon::initialize_ga_randomizer_from_entity_task( rand, ds->entity_task() );

		// done setting up randomizer
		ga.set_rand( rand );
		ga.set_func( func );
		// </stolen code>

		/// Now initialize the GA with a population from which to begin exploration.
		/// Initialize this population completely randomly so that the native sequence
		/// is not arrived at unfairly.

		if ( option[ msd::seed_sequences ].user() ) {
			utility::vector1< std::string > seedseqs = option[ msd::seed_sequences ];
			for ( Size ii = 1; ii <= seedseqs.size(); ++ii ) {
				if ( seedseqs[ ii ].size() != ds->entity_task()->total_residue() ) {
					utility_exit_with_message( "Input seed sequence " + seedseqs[ ii ] + " has " + utility::to_string( seedseqs[ ii ].size() ) + " elements; must have the same number of elements as the number specified in the entity resfile (" + utility::to_string( ds->entity_task()->total_residue()) + ")" );
				}
				ga.add_entity( protocols::multistate_design::entity_elements_from_1letterstring( seedseqs[ ii ] ) );
			}
		}
		if ( option[ msd::seed_sequence_from_input_pdb ].user() ) {
			std::string seq = read_native_sequence_for_entity_elements( ds->entity_task()->total_residue() );
			ga.add_entity( protocols::multistate_design::entity_elements_from_1letterstring( seq ) );
		}

		if ( option[ msd::fill_gen1_from_seed_sequences ] && option[ msd::seed_sequences ].user() ) {
			ga.fill_with_perturbations_of_existing_entities();
		} else {
			ga.fill_with_random_entities();
		}
		// clear parents for the next generation
		// ga.clear_parents(); // do I need this?
		// loop over generations
		while ( !ga.complete() ) {
			clock_t starttime = clock();
			if (ga.current_generation_complete()) ga.evolve_next_generation();
			ga.evaluate_fitnesses();
			if ( TR.visible( t_debug )) {
				TR(t_debug) << "Generation " << ga.current_generation() << ":" << std::endl;
				ga.print_population( TR(t_debug) );
			}
			clock_t stoptime = clock();
			TR << "Generation " << ga.current_generation() << " took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC
				<< " seconds; best fitness = " << ga.best_fitness_from_current_generation() << std::endl;
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
			func->print_load_balance_statistics( TR );
			func->reset_load_balance_statistics();
#endif
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
	}
	return 0;
}


