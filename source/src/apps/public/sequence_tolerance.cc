// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file sequence_tolerance.cc
/// @brief App for predicting the sequence tolerance of a structure or set of structures

#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>
#include <protocols/multistate_design/PackingState.hh>
#include <protocols/multistate_design/MetricCalculatorFitnessFunction.hh>
#include <protocols/multistate_design/MultiStateEntity.hh>
#include <protocols/toolbox/pose_metric_calculators/DecomposeAndReweightEnergiesCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/ResidueDecompositionByChainCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/SurfaceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
using namespace protocols::genetic_algorithm;
using namespace protocols::multistate_design;

//#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/PDBOutput.hh>
//#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
//#include <protocols/dna/RotamerDNAHBondFilter.hh>
#include <protocols/dna/util.hh> // add_constraints_from_file
//using namespace protocols::dna;

#include <protocols/viewer/viewers.hh>

#include <devel/init.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/dna/setup.hh>

// AUTO-REMOVED #include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static thread_local basic::Tracer TR( "app.sequence_tolerance" );

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> // lead_zero_string_of

// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>
#include <list> // PackerTask allowed_residue_types

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/ms.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <utility/excn/Exceptions.hh>

using namespace core;
	using namespace pose::metrics;
using namespace conformation;
using namespace chemical;
using namespace pack;
	using namespace task;
		using namespace operation;
using namespace scoring;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
using namespace protocols::toolbox::pose_metric_calculators;

#include <basic/options/option_macros.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <basic/MetricValue.hh>



OPT_1GRP_KEY(RealVector, seq_tol, fitness_master_weights)
OPT_1GRP_KEY(Boolean, seq_tol, unsat_polars)
OPT_1GRP_KEY(Boolean, seq_tol, surface)

void *
sequence_tolerance_main( void * );

int
main( int argc, char * argv[] )
{
	try {
	NEW_OPT(seq_tol::fitness_master_weights, "fitness function master weights", utility::vector1<core::Real>());
	NEW_OPT(seq_tol::unsat_polars, "calculate the number of buried unsatisfied polar atoms", false);
	NEW_OPT(seq_tol::surface, "calculate the the surface score", false);
	devel::init( argc, argv );
	protocols::viewer::viewer_main( sequence_tolerance_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

// for sorting std::pairs by the second element
template <typename A, typename B>
bool cmpscnd( std::pair<A,B> const & a, std::pair<A,B> const & b )
{
	return a.second < b.second;
}

void *
sequence_tolerance_main( void * )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	basic::prof_reset();

	// get filenames
	typedef utility::vector1< utility::file::FileName > Filenames;
	Filenames pdbnames;

	if ( !option[ in::file::s ].user() ) return 0;
	pdbnames = option[ in::file::s ]().vector();
	// load pdb
	pose::PoseOP pose( new pose::Pose );
	core::import_pose::pose_from_pdb( *pose, pdbnames.front() );
	protocols::dna::add_constraints_from_file( *pose ); // (if specified by options)

	// set up genetic algorithm
	GeneticAlgorithm ga;
	ga.set_max_generations( option[ ms::generations ] );
	ga.set_max_pop_size( option[ ms::pop_size ]() );
	ga.set_num_to_propagate( 1 );
	ga.set_frac_by_recomb( option[ ms::fraction_by_recombination ]() );

	// create an entity template object from which all entities will be cloned
	Entity::OP entity_template( new protocols::multistate_design::MultiStateEntity );
	ga.set_entity_template(entity_template);

	// enable checkpointing and load if the checkpoint files exist
	ga.set_checkpoint_prefix( option[ ms::checkpoint::prefix ]() );
	ga.set_checkpoint_write_interval( option[ ms::checkpoint::interval ]() );
	ga.set_checkpoint_gzip( option[ ms::checkpoint::gz ]() );
	ga.set_checkpoint_rename( option[ ms::checkpoint::rename ]() );
	ga.read_checkpoint();

	// set up the EntityRandomizer
	PositionSpecificRandomizer::OP rand( new PositionSpecificRandomizer );
	rand->set_mutation_rate( option[ ms::mutate_rate ]() );
	rand->set_entity_template(entity_template);

	TaskFactoryOP taskfactory( new TaskFactory );
	taskfactory->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
	if ( option[ packing::resfile ].user() ) {
		taskfactory->push_back( TaskOperationCOP( new operation::ReadResfile ) );
	}
	PackerTaskOP ptask = taskfactory->create_task_and_apply_taskoperations( *pose );

	// figure out design positions/choices from PackerTask
	utility::vector1< Size > design_positions;
	for ( Size i(1), end( ptask->total_residue() ); i <= end; ++i ) {
		// ignore DNA positions
		if ( !pose->residue_type(i).is_protein() ) continue;
		ResidueLevelTask const & rtask( ptask->residue_task(i) );
		if ( rtask.being_designed() ) {
			design_positions.push_back(i);
			// will be passed to randomizer
			EntityElements choices;
			// to avoid duplicate AA's (such as for multiple histidine ResidueTypes)
			std::set< core::chemical::AA > aaset;
			std::list< ResidueTypeCOP > const & allowed( rtask.allowed_residue_types() );
			for ( std::list< ResidueTypeCOP >::const_iterator t( allowed.begin() ), end( allowed.end() );
						t != end; ++t ) {
				core::chemical::AA aa( (*t)->aa() );
				// avoid duplicate AA's (such as for multiple histidine ResidueTypes)
				if ( aaset.find( aa ) != aaset.end() ) continue;
				aaset.insert(aa);
				TR(t_debug) << "adding choice " << aa << std::endl;
				choices.push_back( protocols::genetic_algorithm::EntityElementOP( new PosType( i, aa ) ) );
			}
			rand->append_choices( choices );
		}
	}

	// done setting up the EntityRandomizer
	ga.set_rand( rand );

	// set up the FitnessFunction
	MultiStatePackerOP func( new MultiStatePacker(option[ ms::num_packs ]()) );
	func->set_aggregate_function(MultiStateAggregateFunction::COP( new MultiStateAggregateFunction() ));
	ScoreFunctionOP score_function( core::scoring::get_score_function() );
	func->set_scorefxn( score_function );

	// set up the PoseMetricCalculators to be used by the SingleStateFitnessFunction
	core::pose::metrics::CalculatorFactory & calculator_factory(core::pose::metrics::CalculatorFactory::Instance());

	ResidueDecompositionByChainCalculatorOP res_decomp_calculator( new ResidueDecompositionByChainCalculator() );
	calculator_factory.register_calculator("residue_decomposition", res_decomp_calculator);

	DecomposeAndReweightEnergiesCalculatorOP decomp_reweight_calculator( new DecomposeAndReweightEnergiesCalculator("residue_decomposition") );
	if ( option[ seq_tol::fitness_master_weights ].user() ) {
		decomp_reweight_calculator->master_weight_vector(option[ seq_tol::fitness_master_weights ]());
	}
	calculator_factory.register_calculator("fitness", decomp_reweight_calculator);

	// done setting up the PoseMetricCalculators
	protocols::multistate_design::MetricCalculatorFitnessFunctionOP metric_fitness_function( new protocols::multistate_design::MetricCalculatorFitnessFunction("fitness", "weighted_total") );

	// evaluate and output the fitness of the input pose
	func->scorefxn()->score(*pose);
	TR(t_info) << "Residue Decomposition: " << pose->print_metric("residue_decomposition", "residue_decomposition") << std::endl;
	TR(t_info) << "Fitness Function Master Weights: " << pose->print_metric("fitness", "master_weight_vector") << std::endl;
	TR(t_info) << "Fitness of Starting Structure: " << '\n' << pose->print_metric("fitness", "summary") << std::endl;

	// tell the FitnessFunction to add the vector of components without master weighting to MultiStateEntity objects
	func->add_metric_value_getter(
		"fitness_comp",
		MetricValueGetter("fitness", "weighted_total_no_master_vector", basic::MetricValueBaseCOP( new basic::MetricValue<utility::vector1<Real> > ))
	);

	if ( option[ seq_tol::unsat_polars ] ) {
		// tell the FitnessFunction to add buried unsatisfied polars to MultiStateEntity objects
		calculator_factory.register_calculator("sasa", PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() ));
		calculator_factory.register_calculator("num_hbonds", PoseMetricCalculatorOP( new NumberHBondsCalculator() ));
		calculator_factory.register_calculator("unsat_polars", PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds") ));
		func->add_metric_value_getter(
			"unsat_polars",
			MetricValueGetter("unsat_polars", "all_bur_unsat_polars", basic::MetricValueBaseCOP( new basic::MetricValue<Size> ))
		);
	}

	if ( option[ seq_tol::surface ] ) {
		// tell the FitnessFunction to add the surface score to MultiStateEntity objects
		calculator_factory.register_calculator("surface", PoseMetricCalculatorOP( new SurfaceCalculator() ));
		func->add_metric_value_getter(
			"surface",
			MetricValueGetter("surface", "total_surface", basic::MetricValueBaseCOP( new basic::MetricValue<Real> ))
		);
	}

	// add target state to the FitnessFunction
	PackingStateOP target_state( new PackingState( *pose, true ) );
	target_state->fitness_function( metric_fitness_function );
	target_state->create_packer_data( func->scorefxn(), ptask );
	func->add_state( target_state );

	TR(t_info) << "There are " << func->num_positive_states() << " positive states and "
	           << func->num_negative_states() << " negative states" << std::endl;

	// do single-state designs to find best theoretical single-state energy
	func->single_state_design();

	// done setting up the FitnessFunction
	ga.set_func( func );

	// start the genetic algorithm from scratch if not resuming from a checkpoint
	if (ga.population(ga.current_generation()).size() == 0) {
		// add single-state design sequence(s) to genetic algorithm starting population and parents
		SingleStateCOPs states( func->positive_states() );
		TR(t_info) << "Adding single-state design entities:" << std::endl;
		for ( SingleStateCOPs::const_iterator s( states.begin() ), end( states.end() ); s != end; ++s ) {
			EntityElements traits;
			for ( utility::vector1<Size>::const_iterator i( design_positions.begin() ),
						end( design_positions.end() ); i != end; ++i ) {
				PosType pt( *i, (*s)->pose().residue_type(*i).aa() );
				traits.push_back( protocols::genetic_algorithm::EntityElementOP( new PosType( pt ) ) );
				TR(t_info) << pt.to_string() << " ";
			}
			ga.add_entity( traits );
			ga.add_parent_entity( traits );
			TR(t_info) << std::endl;
		}

		// make more entities by mutation of single-state seeds
		ga.fill_by_mutation( option[ ms::pop_from_ss ]() );
		// the rest are fully random
		ga.fill_with_random_entities();
		// clear parents for the next generation
		ga.clear_parents();
	}

	// loop over generations
	while ( !ga.complete() ) {
		if (ga.current_generation_complete()) ga.evolve_next_generation();
		TR(t_info) << "Generation " << ga.current_generation() << ":" << std::endl;
		ga.evaluate_fitnesses();
		ga.print_generation_statistics(TR(t_info), ga.current_generation());
	}

	// output resulting solution(s)
	protocols::dna::PDBOutputOP pdboutput( new protocols::dna::PDBOutput );
	pdboutput->reference_pose( *pose );
	pdboutput->score_function( *score_function );

	std::string prefix("result");
	if ( option[ OptionKeys::out::prefix ].user() ) prefix = option[ out::prefix ]();

	// sort local copy of sequence/fitness cache
	typedef GeneticAlgorithm::TraitEntityHashMap TraitEntityHashMap;
	TraitEntityHashMap const & cache( ga.entity_cache() );
	utility::vector1< Entity::OP > sortable;
// 	std::copy( cache.begin(), cache.end(), sortable.begin() ); // FAIL(?)
	for ( TraitEntityHashMap::const_iterator it( cache.begin() ), end( cache.end() ); it != end; ++it ) {
		sortable.push_back( it->second );
	}
	std::sort( sortable.begin(), sortable.end(), lt_OP_deref< Entity > );

	TR(t_info) << "Evaluated " << sortable.size() << " sequences out of " << rand->library_size()
	           << " possible.\nBest sequences:\n";
	// list and output top solutions
	int counter(1);
	for ( utility::vector1<Entity::OP>::const_iterator it( sortable.begin() ), end( sortable.end() );
				it != end && counter <= option[ ms::numresults ](); ++it, ++counter ) {
		protocols::genetic_algorithm::Entity & entity(**it);
		// apply sequence to existing positive state(s)
		func->evaluate_positive_states( entity );
		// copy pose
		pose::Pose solution_pose = func->positive_states().front()->pose();
		// output pdb with information
		std::string traitstring;
		for ( core::Size traitnum = 1; traitnum <= entity.traits().size(); ++traitnum ) {
			PosTypeCOP postype = utility::pointer::dynamic_pointer_cast< protocols::multistate_design::PosType const > ( entity.traits()[traitnum] );
			if ( ! postype ) {
				utility_exit_with_message( "PosType dynamic cast failed for entity element with name " + entity.traits()[traitnum]->name() );
			}
			traitstring += core::chemical::oneletter_code_from_aa( postype->type() );
		}
		std::string pdbname( prefix + "_" + lead_zero_string_of(counter,4) + "_" + traitstring + ".pdb" );
		std::list< std::string > extra_lines;
		std::ostringstream ms_info;
		ms_info << "MultiState Fitness: " << F(5,4,entity.fitness());
		extra_lines.push_back( ms_info.str() );
    ms_info.str(""); // funky way to 'empty' ostringstream
    ms_info << "MultiState Fitness Offset: " << F(5,4,entity.fitness() - (*(sortable.begin()))->fitness());
    extra_lines.push_back( ms_info.str() );
		ms_info.str(""); // funky way to 'empty' ostringstream
		ms_info << "MultiState Fitness Offset: " << F(5,4,entity.fitness() - (*(sortable.begin()))->fitness());
		extra_lines.push_back( ms_info.str() );
		ms_info.str(""); // funky way to 'empty' ostringstream
		ms_info << "MultiState Sequence:";
		for ( EntityElements::const_iterator
				pos( entity.traits().begin() ), end( entity.traits().end() );
				pos != end; ++pos ) {
			ms_info << " " << (*pos)->to_string();
			TR(t_info) << (*pos)->to_string() << " ";
		}
		TR(t_info) << "fitness " << F(5,4,entity.fitness()) << '\n';
		extra_lines.push_back( ms_info.str() );
		pdboutput->add_info( "multistate_design", extra_lines, false );
		(*pdboutput)( solution_pose, pdbname );
	}
	TR(t_info) << std::flush;

	basic::prof_show();
	return 0;
}
