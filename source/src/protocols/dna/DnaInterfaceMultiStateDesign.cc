// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DnaInterfaceMultiStateDesign.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/dna/DnaInterfaceMultiStateDesign.hh>
#include <protocols/dna/DnaInterfaceMultiStateDesignCreator.hh>

// Package headers
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/RotamerDNAHBondFilter.hh>
#include <protocols/dna/util.hh> // find_basepairs, substitute_residue
#include <protocols/multistate_design/MultiStatePacker.hh>
#include <protocols/multistate_design/PackingState.hh>
#include <protocols/multistate_design/PartitionAggregateFunction.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh> // string_split
#include <utility/tag/Tag.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/ms.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


namespace protocols {
namespace dna {

using utility::vector1;
using utility::string_split;
using namespace core;
	using namespace chemical;
	using namespace basic::options;
	using namespace pack;
		using namespace task;
			using namespace operation;
	using namespace scoring;

using namespace ObjexxFCL::format;

using namespace multistate_design;
using namespace genetic_algorithm;

using basic::t_error;
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static basic::Tracer TR("protocols.dna.DnaInterfaceMultiStateDesign",t_info);

std::string
DnaInterfaceMultiStateDesignCreator::keyname() const
{
	return DnaInterfaceMultiStateDesignCreator::mover_name();
}

protocols::moves::MoverOP
DnaInterfaceMultiStateDesignCreator::create_mover() const {
	return new DnaInterfaceMultiStateDesign;
}

std::string
DnaInterfaceMultiStateDesignCreator::mover_name()
{
	return "DnaInterfaceMultiStateDesign";
}

DnaInterfaceMultiStateDesign::DnaInterfaceMultiStateDesign()
	: protocols::simple_moves::PackRotamersMover( DnaInterfaceMultiStateDesignCreator::mover_name() ),
		gen_alg_(0),
		multistate_packer_(0),
		dna_chains_(0),
		// option flags/parameters: default to command line options
		// parse_my_tag method may change them
		generations_(               option[ OptionKeys::ms::generations ]() ),
		pop_size_(                  option[ OptionKeys::ms::pop_size ]() ),
		num_packs_(                 option[ OptionKeys::ms::num_packs ]() ),
		pop_from_ss_(               option[ OptionKeys::ms::pop_from_ss ]() ),
		numresults_(                option[ OptionKeys::ms::numresults ]() ),
		fraction_by_recombination_( option[ OptionKeys::ms::fraction_by_recombination ]() ),
		mutate_rate_(               option[ OptionKeys::ms::mutate_rate ]() ),
		boltz_temp_(                option[ OptionKeys::ms::Boltz_temp ]() ),
		anchor_offset_(             option[ OptionKeys::ms::anchor_offset ]() ),
		checkpoint_prefix_(         option[ OptionKeys::ms::checkpoint::prefix ]() ),
		checkpoint_interval_(       option[ OptionKeys::ms::checkpoint::interval ]() ),
		checkpoint_gz_(             option[ OptionKeys::ms::checkpoint::gz ]() ),
		checkpoint_rename_(         option[ OptionKeys::ms::checkpoint::rename ]() )
{}

DnaInterfaceMultiStateDesign::~DnaInterfaceMultiStateDesign(){}

void
DnaInterfaceMultiStateDesign::apply( Pose & pose )
{
	initialize( pose );
	run();
	output_results( pose );
}

void
DnaInterfaceMultiStateDesign::copy_dna_chains( DnaChainsOP dna_chains )
{
	runtime_assert( dna_chains );
	dna_chains_ = new DnaChains( *dna_chains );
}

void
DnaInterfaceMultiStateDesign::copy_targeted_dna( DnaDesignDefOPs const & defs )
{
	targeted_dna_ = defs;
}

void
DnaInterfaceMultiStateDesign::initialize( Pose & pose )
{
	// clear out (any) pre-existing info strings (from previous apply calls)
	info().clear();

	if ( ! dna_chains_ ) {
		dna_chains_ = new DnaChains;
		find_basepairs( pose, *dna_chains_ );
	}

	if ( ! score_function() ) {
		std::string weights_tag("dna");
		if ( option[ OptionKeys::score::weights ].user() )
			weights_tag = option[ OptionKeys::score::weights ]();
		score_function( ScoreFunctionFactory::create_score_function( weights_tag ) );
	}

	// always start with a fresh GeneticAlgorithm
	// important when reusing DnaInterfaceMultistateDesign mover
	gen_alg_ = new GeneticAlgorithm; /// APL does this require a PosType ctor?

	// set up genetic algorithm
	gen_alg_->set_max_generations( generations_ );
	gen_alg_->set_max_pop_size( pop_size_ );
	gen_alg_->set_num_to_propagate( static_cast<core::Size>( 0.5 * pop_size_ ) );
	gen_alg_->set_frac_by_recomb( fraction_by_recombination_ );

	// set up sequence randomizer
	PositionSpecificRandomizer::OP rand = new PositionSpecificRandomizer; /// APL does this require a PosType ctor?
	rand->set_mutation_rate( mutate_rate_ );

	TaskFactoryOP my_tf;
	// if PackRotamerMover base class has no initialized TaskFactory, create default one here
	if ( ! task_factory() ) {
		// DNA-specific TaskFactory -> PackerTask -> figure out positions to design
		my_tf = new TaskFactory;
		my_tf->push_back( new InitializeFromCommandline );
		if ( option[ OptionKeys::packing::resfile ].user() ) my_tf->push_back( new ReadResfile );
		RestrictDesignToProteinDNAInterfaceOP rest_to_dna_int = new RestrictDesignToProteinDNAInterface;
		rest_to_dna_int->copy_dna_chains( dna_chains_ );
		if ( ! targeted_dna_.empty() ) rest_to_dna_int->copy_targeted_dna( targeted_dna_ );
		my_tf->push_back( rest_to_dna_int );
	} else { // TaskFactory already exists, add to it
		// (temporary? parser has no support for RotamerOperations yet)
		my_tf = new TaskFactory( *task_factory() );
	}
	// a protein-DNA hbonding filter for ex rotamers that the PackerTask makes available to the rotamer set during rotamer building (formerly known as 'rotamer explosion')
	RotamerDNAHBondFilterOP rot_dna_hb_filter( new RotamerDNAHBondFilter );
	my_tf->push_back( new AppendRotamer( rot_dna_hb_filter ) );
	task_factory( my_tf ); // PackRotamersMover base class setter

	PackerTaskOP ptask = task_factory()->create_task_and_apply_taskoperations( pose );

	// figure out design positions/choices from PackerTask
	vector1< Size > design_positions;
	for ( Size i(1), end( ptask->total_residue() ); i <= end; ++i ) {
		// ignore DNA positions
		if ( !pose.residue_type(i).is_protein() ) continue;
		ResidueLevelTask const & rtask( ptask->residue_task(i) );
		if ( rtask.being_designed() ) {
			design_positions.push_back(i);
			// will be passed to randomizer
			vector1< EntityElementOP > choices;
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
	gen_alg_->set_rand( rand );
	TR(t_info) << "There will be " << rand->library_size() << " possible sequences." << std::endl;

	// set up fitness function
	multistate_packer_ = new MultiStatePacker( num_packs_ );

	multistate_packer_->set_aggregate_function(
		new PartitionAggregateFunction( boltz_temp_, anchor_offset_ ) );

	multistate_packer_->set_scorefxn( score_function() );

	// add target and competitor states to fitness function
	add_dna_states( pose, ptask );

	TR(t_info) << "There are " << multistate_packer_->num_positive_states() << " positive states and "
	           << multistate_packer_->num_negative_states() << " negative states" << std::endl;

	// do single-state designs to find best theoretical single-state energy
	multistate_packer_->single_state_design();
	// done setting up fitness function
	gen_alg_->set_func( multistate_packer_ );

	// enable checkpointing
	gen_alg_->set_checkpoint_prefix( checkpoint_prefix_ );
	gen_alg_->set_checkpoint_write_interval( checkpoint_interval_ );
	gen_alg_->set_checkpoint_gzip( checkpoint_gz_ );
	gen_alg_->set_checkpoint_rename( checkpoint_rename_ );
	gen_alg_->read_checkpoint();

	// start the genetic algorithm from scratch if not resuming from a checkpoint
	if ( gen_alg_->population( gen_alg_->current_generation() ).size() == 0 ) {
		// add single-state design sequence(s) to genetic algorithm starting population
		SingleStateCOPs states( multistate_packer_->positive_states() );
		TR(t_info) << "Adding single-state design entities:" << std::endl;
		for ( SingleStateCOPs::const_iterator s( states.begin() ), end( states.end() );
					s != end; ++s ) {
			EntityElements traits;
			for ( vector1< Size >::const_iterator i( design_positions.begin() ),
					end( design_positions.end() ); i != end; ++i ) {
				PosTypeOP pt = new PosType( *i, (*s)->pose().residue_type(*i).aa() );
				traits.push_back( pt );
				TR(t_info) << pt->to_string() << " ";
			}
			gen_alg_->add_entity( traits );
			gen_alg_->add_parent_entity( traits );
			TR(t_info) << std::endl;
		}

		// make more entities by mutation of single-state seeds
		gen_alg_->fill_by_mutation( pop_from_ss_ );
		// the rest are fully random
		gen_alg_->fill_with_random_entities();
		// clear parents for the next generation
		gen_alg_->clear_parents();
	}
}

void
DnaInterfaceMultiStateDesign::run()
{
	// loop over generations
	while ( !gen_alg_->complete() ) {
		if ( gen_alg_->current_generation_complete() ) gen_alg_->evolve_next_generation();
		TR(t_info) << "Generation " << gen_alg_->current_generation() << ":" << std::endl;
		gen_alg_->evaluate_fitnesses();
		gen_alg_->print_population( TR(t_info) );
	}
}

void
DnaInterfaceMultiStateDesign::output_results( Pose & pose )
{
	protocols::dna::PDBOutputOP pdboutput = new protocols::dna::PDBOutput;
	pdboutput->score_function( *multistate_packer_->scorefxn() );
	pdboutput->reference_pose( pose );

	std::string prefix("result");
	if ( option[ OptionKeys::out::prefix ].user() ) prefix = option[ OptionKeys::out::prefix ]();

	// sort local copy of sequence/fitness cache
	typedef GeneticAlgorithm::TraitEntityHashMap TraitEntityHashMap;
	TraitEntityHashMap const & cache( gen_alg_->entity_cache() );
	vector1< EntityOP > sortable;
// 	std::copy( cache.begin(), cache.end(), sortable.begin() ); // FAIL(?)
	for ( TraitEntityHashMap::const_iterator it( cache.begin() ), end( cache.end() );
				it != end; ++it ) {
		sortable.push_back( it->second );
	}
	std::sort( sortable.begin(), sortable.end(), lt_OP_deref< Entity > );

	TR(t_info) << "Evaluated " << sortable.size() << " sequences.\nBest sequences:\n";
	// list and output top solutions
	Size counter(0);

	for ( vector1< EntityOP >::const_iterator it( sortable.begin() ),
				end( sortable.end() ); it != end; ++it ) {
		Entity & entity(**it);
		// apply sequence to existing positive state(s)
		multistate_packer_->evaluate_positive_states( entity );
		// copy pose
		Pose solution_pose = multistate_packer_->positive_states().front()->pose();
		// output pdb with information
		std::string pdbname( prefix + "_ms_" + ObjexxFCL::lead_zero_string_of(counter,4) + ".pdb" );
		Strings extra_lines;
		std::ostringstream ms_info;
		ms_info << "REMARK MultiState Fitness: " << F(5,4,entity.fitness());
		extra_lines.push_back( ms_info.str() );
		ms_info.str(""); // funky way to 'empty' ostringstream
		ms_info << "REMARK MultiState Sequence:";
		for ( EntityElements::const_iterator
				pos( entity.traits().begin() ), end( entity.traits().end() );
				pos != end; ++pos ) {
			ms_info << " " << (*pos)->to_string();
			TR(t_info) << (*pos)->to_string() << " ";
		}
		TR(t_info) << "fitness " << F(5,4,entity.fitness()) << '\n';
		extra_lines.push_back( ms_info.str() );
		if ( counter == 0 ) {
			// copy top result to input pose
			pose = solution_pose;
			// save info for top result
			info().insert( info().end(), extra_lines.begin(), extra_lines.end() );
			// ensure that pose has up-to-date score information
			(*multistate_packer_->scorefxn())(pose);
		}
		// set numresults to 0 to suppress output
		if ( counter >= numresults_ ) break;
		pdboutput->add_info( "multistate_design", extra_lines, false );
		(*pdboutput)( solution_pose, pdbname );
		++counter;
	}
	TR(t_info) << std::endl;
}

std::string
DnaInterfaceMultiStateDesign::get_name() const {
	return DnaInterfaceMultiStateDesignCreator::mover_name();
}

///@brief parse "XML" Tag (specifically in the context of the parser/scripting scheme)
void DnaInterfaceMultiStateDesign::parse_my_tag(
	TagPtr const tag,
	moves::DataMap & datamap,
	protocols::filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose
)
{
	// flags/parameters (override options settings)
	if ( tag->hasOption("generations") ) generations_ = tag->getOption<Size>("generations");
	if ( tag->hasOption("pop_size") ) pop_size_ = tag->getOption<Size>("pop_size");
	if ( tag->hasOption("num_packs") ) num_packs_ = tag->getOption<Size>("num_packs");
	if ( tag->hasOption("pop_from_ss") ) pop_from_ss_ = tag->getOption<Size>("pop_from_ss");
	if ( tag->hasOption("numresults") ) numresults_ = tag->getOption<Size>("numresults");
	if ( tag->hasOption("fraction_by_recombination") )
		fraction_by_recombination_ = tag->getOption<Real>("fraction_by_recombination");
	if ( tag->hasOption("mutate_rate") ) mutate_rate_ = tag->getOption<Real>("mutate_rate");
	if ( tag->hasOption("boltz_temp") ) boltz_temp_ = tag->getOption<Real>("boltz_temp");
	if ( tag->hasOption("anchor_offset") ) anchor_offset_ = tag->getOption<Real>("anchor_offset");
		// checkpointing options
	if ( tag->hasOption("checkpoint_prefix") )
		checkpoint_prefix_ = tag->getOption<std::string>("checkpoint_prefix");
	if ( tag->hasOption("checkpoint_interval") )
		checkpoint_interval_ = tag->getOption<Size>("checkpoint_interval");
	if ( tag->hasOption("checkpoint_gz") ) checkpoint_gz_ = tag->getOption<bool>("checkpoint_gz");
	if ( tag->hasOption("checkpoint_rename") )
		checkpoint_rename_ = tag->getOption<bool>("checkpoint_rename");

	// calls to PackRotamersMover base class methods
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

///@brief required in the context of the parser/scripting scheme
moves::MoverOP
DnaInterfaceMultiStateDesign::fresh_instance() const
{
	return new DnaInterfaceMultiStateDesign;
}

///@brief required in the context of the parser/scripting scheme
moves::MoverOP
DnaInterfaceMultiStateDesign::clone() const
{
	return new DnaInterfaceMultiStateDesign( *this );
}

void
DnaInterfaceMultiStateDesign::add_dna_states(
	Pose const & pose,
	PackerTaskCOP ptask
)
{
	runtime_assert( dna_chains_ );
	runtime_assert( multistate_packer_ );

	// temporary copy of Pose used to build DNA target and competitor states
	Pose mutpose( pose );
	ResidueTypeSet const & rts( mutpose.residue(1).residue_type_set() );
	ResidueTypeCOPs dna_types(
		ResidueSelector().set_property("DNA").exclude_variants().select( rts )
	);

	TR(t_info) << "\nBuilding dna target state:\n";

	// construct and add the target state
	for ( DnaPositions::const_iterator itr( dna_chains_->begin() );
				itr != dna_chains_->end(); ++itr ) {
		// limit to dna design positions
		DnaPosition const & pos( itr->second );
		Size const index( pos.top() );
		runtime_assert( index == itr->first );
		// resfile key "TARGET" indicates positions at which multistate design will be targeted
		if ( !ptask->has_behavior("TARGET",index) ) continue;

		std::ostringstream pdbtag;
		if ( pose.pdb_info() ) {
			pdbtag << pose.pdb_info()->chain(index) << '.' << pose.pdb_info()->number(index);
		} else {
			pdbtag << pose.chain(index) << '.' << index;
		}
		ResidueTypeCOP target_type( ptask->residue_task( index ).target_type() );
		if ( ! target_type ) {
			TR(t_error) << "No target type found for " << pdbtag.str() << '\n'
				<< "(Did the DNA definition string indicate a target nucleotide type?" << std::endl;
			utility_exit();
		}
		runtime_assert( pos.paired() );
		std::ostringstream pdbtag_btm;
		if ( pose.pdb_info() ) {
			pdbtag_btm << pose.pdb_info()->chain( pos.bottom() )
				<< '.' << pose.pdb_info()->number( pos.bottom() );
		} else {
			pdbtag_btm << pose.chain( pos.bottom() ) << '.' << pos.bottom();
		}
		ResidueTypeCOP bot_type( ptask->residue_task( pos.bottom() ).target_type() );
		if ( ! target_type ) {
			TR(t_error) << "No target type found for " << pdbtag_btm.str() << '\n'
				<< "(Did the DNA definition string indicate a target nucleotide type?" << std::endl;
			utility_exit();
		}
		TR(t_info) << pdbtag.str() << '.' << target_type->name() << '/';
		TR(t_info) << pdbtag_btm.str() << '.' << bot_type->name() << ", ";
		substitute_residue( mutpose, index, *target_type );
		substitute_residue( mutpose, pos.bottom(), *bot_type );
	}
	TR(t_info) << '\n';
	PackingStateOP target_state = new PackingState( mutpose, true );
	target_state->create_packer_data( multistate_packer_->scorefxn(), ptask );
	multistate_packer_->add_state( target_state );

	TR(t_info) << "Building dna sequence competitors:\n";

	// build competitor DNA negative states and add to the MultiStateDesign instance
	for ( DnaPositions::const_iterator itr( dna_chains_->begin() );
	      itr != dna_chains_->end(); ++itr ) {
		// limit to dna design positions
		DnaPosition const & pos( itr->second );
		Size const index( pos.top() );
		assert( index == itr->first );
		if ( !ptask->has_behavior("TARGET",index) ) continue;
		if ( pose.pdb_info() ) {
			TR(t_info) << pose.pdb_info()->number( index ) << "/"
				<< pose.pdb_info()->number( pos.bottom() ) << ":";
		} else {
			TR(t_info) << index << "/" << pos.bottom() << ":";
		}

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
			multistate_packer_->add_state( competitor_state );
		}

		// restore the original basepair at this position
		substitute_residue( mutpose, index, orig_top );
		substitute_residue( mutpose, pos.bottom(), orig_bot );
	}
	TR(t_info) << std::endl;
}

} // namespace dna
} // namespace protocols
