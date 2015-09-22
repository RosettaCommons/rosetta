// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ga.cc
/// @brief ashworth

#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>
#include <protocols/multistate_design/PackingState.hh>
#include <protocols/multistate_design/PartitionAggregateFunction.hh>
using namespace protocols::genetic_algorithm;
using namespace protocols::multistate_design;

#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/RotamerDNAHBondFilter.hh>
#include <protocols/dna/util.hh>
using namespace protocols::dna;

#include <protocols/viewer/viewers.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/dna/setup.hh>

#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
using utility::string_split;

#include <basic/prof.hh>
#include <basic/Tracer.hh>
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "app.ga", t_info );

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> // lead_zero_string_of

#include <fstream>
#include <iostream>
#include <string>
#include <list> // PackerTask allowed_residue_types

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/ms.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace pack;
	using namespace task;
		using namespace operation;
using namespace scoring;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

void
add_dna_states(
	MultiStatePacker & func,
	pose::Pose const & pose,
	PackerTaskCOP ptask
)
{
	// temporary Pose to build states
	pose::Pose mutpose( pose );
	ResidueTypeSet const & rts( mutpose.residue(1).residue_type_set() );
	ResidueTypeCOPs dna_types(
		ResidueTypeSelector().set_property("DNA").exclude_variants().select( rts )
	);

	// get DNA chains info
	DnaChains dna_chains;
	find_basepairs( mutpose, dna_chains );

	TR(t_info) << "\nBuilding dna target state:\n";

	// construct and add the target state
	for ( DnaPositions::const_iterator itr( dna_chains.begin() ); itr != dna_chains.end(); ++itr ) {
		// limit to dna design positions
		DnaPosition const & pos( itr->second );
		Size const index( pos.top() );
		assert( index == itr->first );
		// resfile key "TARGET" indicates positions at which multistate design will be targeted
		if ( !ptask->has_behavior("TARGET",index) ) continue;

		ResidueTypeCOP target_type( ptask->residue_task( index ).target_type() );
		runtime_assert( pos.paired() );
		ResidueTypeCOP bot_type( ptask->residue_task( pos.bottom() ).target_type() );
		TR(t_info) << mutpose.pdb_info()->chain( index ) << '.'
		           << mutpose.pdb_info()->number( index ) << '.' << target_type->name() << '/'
		           << mutpose.pdb_info()->chain( pos.bottom() ) << '.'
		           << mutpose.pdb_info()->number( pos.bottom() ) << '.' << bot_type->name() << ", ";
		substitute_residue( mutpose, index, *target_type );
		substitute_residue( mutpose, pos.bottom(), *bot_type );
	}
	TR(t_info) << '\n';
	PackingStateOP target_state = new PackingState( mutpose, true );
	target_state->create_packer_data( func.scorefxn(), ptask );
	func.add_state( target_state );

	TR(t_info) << "Building dna sequence competitors:\n";

	// build competitor DNA negative states and add to the MultiStateDesign instance
	for ( DnaPositions::const_iterator itr( dna_chains.begin() ); itr != dna_chains.end(); ++itr ) {
		// limit to dna design positions
		DnaPosition const & pos( itr->second );
		Size const index( pos.top() );
		assert( index == itr->first );
		if ( !ptask->has_behavior("TARGET",index) ) continue;

		TR(t_info) << mutpose.pdb_info()->number( index ) << "/" << mutpose.pdb_info()->number( pos.bottom() ) << ":";

		// remember the starting type in order to restore later
		ResidueType const & orig_top( mutpose.residue_type( index ) );
		ResidueType const & orig_bot( mutpose.residue_type( pos.bottom() ) );

		// add a negative state for every single basepair substitution
		for ( ResidueTypeCOPs::const_iterator rt( dna_types.begin() ); rt != dna_types.end(); ++rt ) {
			std::string name( (*rt)->name() );
			if ( (*rt)->name3() == orig_top.name3() ) continue; // add mutants only
			ResidueType const & bot_type( rts.name_map( dna_comp_name_str( name ) ) );
			TR(t_info) << " " << name << "/" << bot_type.name() << ",";
			substitute_residue( mutpose, index, **rt );
			substitute_residue( mutpose, pos.bottom(), bot_type );
			PackingStateOP competitor_state = new PackingState( mutpose, false );
			competitor_state->share_packer_data_from( *target_state );
			func.add_state( competitor_state );
		}

		// restore the original basepair at this position
		substitute_residue( mutpose, index, orig_top );
		substitute_residue( mutpose, pos.bottom(), orig_bot );
	}
	TR(t_info) << std::endl;
}

// for sorting std::pairs by the second element
template <typename A, typename B>
bool cmpscnd( std::pair<A,B> const & a, std::pair<A,B> const & b )
{
	return a.second < b.second;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void *
ga_main( void * )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	basic::prof_reset();

	// get filenames
	typedef vector1< utility::file::FileName > Filenames;
	Filenames pdbnames;

	if ( !option[ in::file::s ].user() ) return 0;
	pdbnames = option[ in::file::s ]().vector();
	// load pdb
	pose::PoseOP pose = new pose::Pose;
	core::import_pose::pose_from_pdb( *pose, pdbnames.front() );
	add_constraints_from_file( *pose ); // (if specified by options)

	// initialization necessary for scoring of (with?) DNA
	scoring::dna::set_base_partner( *pose );

	// set up genetic algorithm
	GeneticAlgorithm ga;
	ga.set_max_generations( option[ ms::generations ] );
	ga.set_max_pop_size( option[ ms::pop_size ]() );
	ga.set_num_to_propagate( static_cast<core::Size>(0.5*option[ ms::pop_size ]) );
	ga.set_frac_by_recomb( option[ ms::fraction_by_recombination ]() );

	// set up sequence randomizer
	PositionSpecificRandomizer::OP rand = new PositionSpecificRandomizer;
	rand->set_mutation_rate( option[ ms::mutate_rate ]() );

	TaskFactoryOP taskfactory = new TaskFactory;
	taskfactory->push_back( new InitializeFromCommandline );
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		taskfactory->push_back( new ReadResfile );
	}
	// DNA-specific TaskFactory
	taskfactory->push_back( new RestrictDesignToProteinDNAInterface );
	taskfactory->push_back( new AppendRotamer( new RotamerDNAHBondFilter ) );
	PackerTaskOP ptask = taskfactory->create_task_and_apply_taskoperations( *pose );

	// figure out design positions/choices from PackerTask
	vector1< Size > design_positions;
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
				choices.push_back( new PosType( i, aa ) );
			}
			rand->append_choices( choices );
		}
	}
	// done setting up randomizer
	ga.set_rand( rand );

	// set up fitness function
	MultiStatePackerOP func = new MultiStatePacker(
		option[ ms::num_packs ]()
	);

	func->set_aggregate_function(
		new PartitionAggregateFunction(option[ ms::Boltz_temp ](), option[ ms::anchor_offset ]())
	);

	// get scorefunction, give to fitness function
	std::string weights_tag("dna");
	if ( option[ score::weights ].user() ) weights_tag = option[ score::weights ]();
	ScoreFunctionOP score_function( ScoreFunctionFactory::create_score_function( weights_tag ) );
	func->set_scorefxn( score_function );

	// add target and competitor states to fitness function
	add_dna_states( *func, *pose, ptask );

	TR(t_info) << "There are " << func->num_positive_states() << " positive states and "
	           << func->num_negative_states() << " negative states" << std::endl;

	// do single-state designs to find best theoretical single-state energy
	func->single_state_design();
	// done setting up fitness function
	ga.set_func( func );

	// enable checkpointing
	ga.set_checkpoint_prefix( option[ ms::checkpoint::prefix ]() );
	ga.set_checkpoint_write_interval( option[ ms::checkpoint::interval ]() );
	ga.set_checkpoint_gzip( option[ ms::checkpoint::gz ]() );
	ga.set_checkpoint_rename( option[ ms::checkpoint::rename ]() );
	ga.read_checkpoint();

	// start the genetic algorithm from scratch if not resuming from a checkpoint
	if (ga.population(ga.current_generation()).size() == 0) {
		// add single-state design sequence(s) to genetic algorithm starting population
		SingleStateCOPs states( func->positive_states() );
		TR(t_info) << "Adding single-state design entities:" << std::endl;
		for ( SingleStateCOPs::const_iterator s( states.begin() ), end( states.end() ); s != end; ++s ) {
			EntityElements traits;
			for ( vector1<Size>::const_iterator i( design_positions.begin() ),
						end( design_positions.end() ); i != end; ++i ) {
				PosType pt( *i, (*s)->pose().residue_type(*i).aa() );
				traits.push_back( new PosType( pt ) );
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
		ga.print_population( TR(t_info) );
	}

	// output resulting solution(s)
	protocols::dna::PDBOutputOP pdboutput = new protocols::dna::PDBOutput;
	pdboutput->reference_pose( *pose );
	pdboutput->score_function( *score_function );

	std::string prefix("result");
	if ( option[ OptionKeys::out::prefix ].user() ) prefix = option[ out::prefix ]();

	// sort local copy of sequence/fitness cache
	typedef GeneticAlgorithm::TraitEntityHashMap TraitEntityHashMap;
	TraitEntityHashMap const & cache( ga.entity_cache() );
	utility::vector1<Entity::OP> sortable;
// 	std::copy( cache.begin(), cache.end(), sortable.begin() ); // FAIL(?)
	for ( TraitEntityHashMap::const_iterator it( cache.begin() ), end( cache.end() ); it != end; ++it ) {
		sortable.push_back( it->second );
	}
	std::sort( sortable.begin(), sortable.end(), lt_OP_deref< Entity > );

	TR(t_info) << "Evaluated " << sortable.size() << " sequences out of " << rand->library_size()
	           << " possible.\nBest sequences:\n";
	// list and output top solutions
	int counter(0);
	for ( utility::vector1<Entity::OP>::const_iterator it( sortable.begin() ), end( sortable.end() );
				it != end && counter < option[ ms::numresults ](); ++it, ++counter ) {
		protocols::genetic_algorithm::Entity & entity(**it);
		// apply sequence to existing positive state(s)
		func->evaluate_positive_states( entity );
		// copy pose
		pose::Pose solution_pose = func->positive_states().front()->pose();
		// output pdb with information
		std::string pdbname( prefix + "_ms_" + lead_zero_string_of(counter,4) + ".pdb" );
		std::list< std::string > extra_lines;
		std::ostringstream ms_info;
		ms_info << "MultiState Fitness: " << F(5,4,entity.fitness());
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

////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv[] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		devel::init( argc, argv );
		protocols::viewer::viewer_main( ga_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
