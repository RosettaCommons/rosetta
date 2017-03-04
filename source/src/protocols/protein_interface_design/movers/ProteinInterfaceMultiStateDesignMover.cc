// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ProteinInterfaceMultiStateDesignMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@uw.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/ProteinInterfaceMultiStateDesignMover.hh>
#include <protocols/protein_interface_design/movers/ProteinInterfaceMultiStateDesignMoverCreator.hh>

// Project headers
#include <core/import_pose/import_pose.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>
#include <protocols/multistate_design/PackingState.hh>
#include <protocols/multistate_design/PartitionAggregateFunction.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.hh>

#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/ms.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using utility::vector1;
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
static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.ProteinInterfaceMultiStateDesignMover", t_info );

// XRW TEMP std::string
// XRW TEMP ProteinInterfaceMultiStateDesignMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ProteinInterfaceMultiStateDesignMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ProteinInterfaceMultiStateDesignMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ProteinInterfaceMultiStateDesignMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ProteinInterfaceMultiStateDesignMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ProteinInterfaceMS";
// XRW TEMP }

ProteinInterfaceMultiStateDesignMover::ProteinInterfaceMultiStateDesignMover() :
	protocols::simple_moves::PackRotamersMover( ProteinInterfaceMultiStateDesignMover::mover_name() ),
	gen_alg_(/* 0 */),
	multistate_packer_(/* 0 */),
	// option flags/parameters: default to command line options
	// parse_my_tag method may change them
	rb_jump_( 1 ),
	scorefxn_( /* 0 */ ),
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
	checkpoint_rename_(         option[ OptionKeys::ms::checkpoint::rename ]() ),
	unbound_( true ),
	unfolded_( true ),
	input_is_positive_( true ),
	use_unbound_for_sequence_profile_( true ),
	bump_threshold_( 1.0 ),
	compare_energy_to_ground_state_( false ),
	fname_prefix_("")
{
	state_poses_.clear();
	saved_state_poses_.clear();
	state_positive_.clear();
	state_unfolded_.clear();
	state_unbound_.clear();
}

ProteinInterfaceMultiStateDesignMover::~ProteinInterfaceMultiStateDesignMover(){}

void
part_complex( core::pose::PoseOP pose, core::Size const rb_jump ){
	using namespace protocols::moves;
	rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( *pose, rb_jump ) );
	translate->step_size( 1000.0 );
	translate->apply( *pose );
}

void
unfold_complex( core::pose::PoseOP pose ){
	for ( core::Size i( 1 ), end( pose->size() ); i<=end; ++i ) {
		// unfold the unfolded protein
		if ( !pose->residue_type(i).is_protein() ) continue;
		std::string const restype( pose->residue(i).type().name() );
		if ( restype == "PRO" || ( i>1 && pose->residue(i-1).type().name() == "PRO" ) ) continue; // don't unfold prolines and preprolines
		/// numbers for extended from Nobu
		pose->set_phi( i, -130.0 );
		pose->set_psi( i,  130.0 );
	}
}


void
ProteinInterfaceMultiStateDesignMover::apply( Pose & pose )
{
	initialize( pose );
	run();
	output_results( pose );
	output_alternative_states( pose );
}

// XRW TEMP std::string
// XRW TEMP ProteinInterfaceMultiStateDesignMover::get_name() const {
// XRW TEMP  return ProteinInterfaceMultiStateDesignMover::mover_name();
// XRW TEMP }

void
ProteinInterfaceMultiStateDesignMover::restrict_sequence_profile(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP const ptask
) const
{
	using namespace core::chemical;
	using namespace core::pack;
	using namespace core::pack::task;
	using core::pose::Pose;
	using core::pose::PoseOP;
	using utility::vector1;
	using std::list;
	using std::vector;
	using protocols::simple_filters::EnergyPerResidueFilter;
	using core::scoring::fa_rep;

	unsigned long const seq_space_before( sequence_space( ptask ) );
	TR<<"Total number of sequence possibilites: "<<seq_space_before<<std::endl;
	if ( !use_unbound_for_sequence_profile_ ) {
		return;
	}
	TR<<"Restricting the packer task to residues that would not clash in the unbound monomer..."<<std::endl;
	/// Turn all designable positions to ala. Freeze everything else
	PoseOP ala_pose( new Pose( pose ) );
	part_complex( ala_pose, rb_jump_ );
	PackerTaskOP ala_task = ptask->clone();
	vector1< bool > allow_ala( num_canonical_aas, false );
	allow_ala[ aa_ala ] = true;
	vector< Size > designable; //which residues are designable? used below
	designable.clear();
	for ( Size i( 1 ), end( ala_task->total_residue() ); i <= end; ++i ) {
		if ( !pose.residue_type( i ).is_protein() ) continue;
		ResidueLevelTask & rtask( ala_task->nonconst_residue_task(i) );
		if ( !rtask.being_designed() ) { // undesignable
			rtask.prevent_repacking();
		} else { //designable
			rtask.restrict_absent_canonical_aas( allow_ala );
			designable.push_back( i );
		}
	}
	pack_rotamers( *ala_pose, *scorefxn_, ala_task );

	// scorefxn for the bump check
	core::scoring::ScoreFunctionOP bump_scorefxn( new core::scoring::ScoreFunction );
	bump_scorefxn->reset();
	bump_scorefxn->set_weight( fa_rep, 1.0 );

	for ( Size const pos : designable ) {
		EnergyPerResidueFilter const eprf( pos, bump_scorefxn, fa_rep, 0 );
		core::Real const ref_bump_energy( eprf.compute( *ala_pose ) );
		PackerTaskOP template_substitution_task( ptask->clone() ); //prevent repacking at all but positions but pos
		for ( Size i( 1 ); i<=pose.size(); ++i ) {
			if ( i!=pos ) {
				template_substitution_task->nonconst_residue_task(i).prevent_repacking();
			}
		}

		ResidueLevelTask & rtask( ptask->nonconst_residue_task( pos ) );
		list< ResidueTypeCOP > const & allowed( rtask.allowed_residue_types() );
		Pose ala_pose_and_single_residue( *ala_pose );
		vector1< bool > allowed_aas_in_pos( num_canonical_aas, false );
		for ( ResidueTypeCOP const t : allowed ) {
			AA const aa( t->aa() );
			PackerTaskOP specific_substitution_task( template_substitution_task->clone() );
			utility::vector1< bool > allow_aa( num_canonical_aas, false );
			allow_aa[ aa ] = true;
			specific_substitution_task->nonconst_residue_task(pos).restrict_absent_canonical_aas( allow_aa );
			rotamer_trials( ala_pose_and_single_residue, *bump_scorefxn, specific_substitution_task );
			core::Real const bump_energy( eprf.compute( ala_pose_and_single_residue ) );
			if ( bump_energy - ref_bump_energy <= bump_threshold_ ) {
				allowed_aas_in_pos[ aa ] = true;
			}
		}///foreach ResidueTypeCOP const t
		rtask.restrict_absent_canonical_aas( allowed_aas_in_pos );
		AA const aa_in_pose( pose.residue( pos ).aa() );
		if ( !allowed_aas_in_pos[ aa_in_pose ] ) {
			TR.Fatal <<"Native identity "<<pose.residue( pos ).name3()<<" at position "<<pos<<" in input pdb is not allowed by bump_test! Increase the bump_cutoff from the current "<<bump_threshold_<<std::endl;
			utility_exit();
		}//fi
	}//foreach pos
	unsigned long const seq_space_after( sequence_space( ptask ) );
	TR<<"Finished restricting. Total number of sequences after restriction: "<<seq_space_after<<'\n';
	TR<<"Orders of magnitude change: "<<log10( (double)seq_space_after ) - log10( (double)seq_space_before )<<std::endl; // REQUIRED FOR WINDOWS
}

unsigned long
ProteinInterfaceMultiStateDesignMover::sequence_space( core::pack::task::PackerTaskCOP ptask ) const
{
	using namespace core::pack::task;

	unsigned long size( 1 );
	for ( core::Size i( 1 ); i<=ptask->total_residue(); ++i ) {
		ResidueLevelTask const & rtask( ptask->residue_task( i ) );
		if ( !rtask.being_designed() ) continue;
		core::Size const pos_allowed( rtask.allowed_residue_types().size() );
		size *= pos_allowed;
	}
	return( size );
}

void
ProteinInterfaceMultiStateDesignMover::initialize( Pose & pose )
{
	// clear out (any) pre-existing info strings (from previous apply calls)
	info().clear();

	// always start with a fresh GeneticAlgorithm
	// important when reusing ProteinInterfaceMultistateDesign mover
	gen_alg_ = GeneticAlgorithmOP( new GeneticAlgorithm );

	// set up genetic algorithm
	gen_alg_->set_max_generations( generations_ );
	gen_alg_->set_max_pop_size( pop_size_ );
	gen_alg_->set_num_to_propagate( static_cast<core::Size>( 0.5 * pop_size_ ) );
	gen_alg_->set_frac_by_recomb( fraction_by_recombination_ );

	// set up sequence randomizer
	PositionSpecificRandomizer::OP rand( new PositionSpecificRandomizer );
	rand->set_mutation_rate( mutate_rate_ );

	TaskFactoryOP my_tf;
	// if PackRotamerMover base class has no initialized TaskFactory, create default one here
	if ( ! task_factory() ) {
		// Protein-interface design-specific TaskFactory -> PackerTask -> figure out positions to design
		my_tf = TaskFactoryOP( new TaskFactory );
	} else { // TaskFactory already exists, add to it
		my_tf = TaskFactoryOP( new TaskFactory( *task_factory() ) );
	}

	my_tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	if ( option[ OptionKeys::packing::resfile ].user() ) my_tf->push_back( TaskOperationCOP( new ReadResfile ) );

	task_factory( my_tf ); // PackRotamersMover base class setter

	PackerTaskOP ptask = task_factory()->create_task_and_apply_taskoperations( pose );

	// figure out design positions/choices from PackerTask
	restrict_sequence_profile( pose, ptask );
	vector1< Size > design_positions;
	for ( Size i(1), end( ptask->total_residue() ); i <= end; ++i ) {
		if ( !pose.residue_type(i).is_protein() ) continue;
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
	// done setting up randomizer
	gen_alg_->set_rand( rand );
	TR(t_info) << "There will be " << rand->library_size() << " possible sequences." << std::endl;

	// set up fitness function
	multistate_packer_ = multistate_design::MultiStatePackerOP( new MultiStatePacker( num_packs_ ) );

	multistate_packer_->set_aggregate_function(
		MultiStateAggregateFunction::COP( MultiStateAggregateFunction::OP( new PartitionAggregateFunction( boltz_temp_, anchor_offset_, compare_energy_to_ground_state_ ) ) ) );

	multistate_packer_->set_scorefxn( scorefxn_ );

	// add target and competitor to fitness function
	add_states( pose );

	TR(t_info) << "There are " << multistate_packer_->num_positive_states() << " positive states and "
		<< multistate_packer_->num_negative_states() << " negative states" << std::endl;

	// do single-state designs to find best theoretical single-state energy
	multistate_packer_->single_state_design( );

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
			for ( vector1<Size>::const_iterator i( design_positions.begin() ),
					end( design_positions.end() ); i != end; ++i ) {
				PosType pt( *i, (*s)->pose().residue_type(*i).aa() );
				traits.push_back( protocols::genetic_algorithm::EntityElementOP( new PosType( pt ) ));
				TR(t_info) << pt.to_string() << " ";
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
	task( ptask ); //For InterfaceSequenceRecapitulation ot konw of the designable task
}

void
ProteinInterfaceMultiStateDesignMover::run()
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
ProteinInterfaceMultiStateDesignMover::output_alternative_states( core::pose::Pose const & output_pose ) const
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP tf( new TaskFactory( *task_factory() ) );// Allow all repackable residues to move, but not redesign
	tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	PackerTaskCOP ptask_output_pose = tf->create_task_and_apply_taskoperations( output_pose );
	std::string const output_pose_fname( fname_prefix_ + "_ms_pos_0000.pdb" );
	core::pose::Pose copy_pose( output_pose );
	pack_rotamers( copy_pose, *scorefxn_, ptask_output_pose );
	(*scorefxn_)(copy_pose );
	copy_pose.dump_scored_pdb( output_pose_fname, *scorefxn_ );

	if ( saved_state_poses_.size() == 0 ) return;
	runtime_assert( saved_state_poses_.size() == state_positive_.size() );
	for ( core::Size count( 1 ); count<=saved_state_poses_.size(); ++count ) {
		core::pose::Pose state_i( *saved_state_poses_[ count ] );
		(*scorefxn_)(state_i); // to set up the energy graph or else interface task operation will work well. arggg!
		PackerTaskCOP unmodifed_ptask( task_factory()->create_task_and_apply_taskoperations( state_i ));
		for ( core::Size resi( 1 ); resi<=state_i.size(); ++resi ) {//thread output-pose's sequence on alternative state
			if ( !unmodifed_ptask->residue_task( resi ).being_designed() ) continue;
			state_i.replace_residue( resi, output_pose.residue( resi ), true );
		}
		state_i.update_residue_neighbors();
		PackerTaskCOP ptask_statei = tf->create_task_and_apply_taskoperations( state_i );
		pack_rotamers( state_i, *scorefxn_, ptask_statei );
		(*scorefxn_)( state_i );

		std::string const neg_pos( state_positive_[ count ] ? "_pos_" : "_neg_" );

		std::string const pdbname( fname_prefix_ + "_ms" + neg_pos + ObjexxFCL::lead_zero_string_of( count, 4 ) + ".pdb" );
		state_i.dump_scored_pdb( pdbname, *scorefxn_ );
		TR<<"\nDumped "<<pdbname<<'\n';
	}
	TR.flush();
}

void
ProteinInterfaceMultiStateDesignMover::output_results( Pose & pose )
{
	protocols::dna::PDBOutputOP pdboutput( new protocols::dna::PDBOutput );
	pdboutput->score_function( *scorefxn_ );
	pdboutput->reference_pose( pose );

	std::string prefix("result");
	if ( option[ OptionKeys::out::prefix ].user() ) prefix = option[ OptionKeys::out::prefix ]();

	typedef GeneticAlgorithm::TraitEntityHashMap TraitEntityHashMap;
	TraitEntityHashMap const & cache( gen_alg_->entity_cache() );
	utility::vector1<Entity::OP> sortable;
	for ( TraitEntityHashMap::const_iterator it( cache.begin() ), end( cache.end() ); it != end; ++it ) {
		sortable.push_back( it->second );
	}
	std::sort( sortable.begin(), sortable.end(), lt_OP_deref< Entity > );

	TR(t_info) << "Evaluated " << sortable.size() << " sequences.\nBest sequences:\n";
	// list and output top solutions
	Size counter(0);

	for ( vector1< Entity::OP >::const_iterator it( sortable.begin() ),
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
		for ( EntityElements::const_iterator pos( entity.traits().begin() ), end( entity.traits().end() );
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
			(*scorefxn_)(pose);
		}
		// set numresults to 0 to suppress output
		if ( counter >= numresults_ ) break;
		pdboutput->add_info( "multistate_design", extra_lines, false );
		(*pdboutput)( solution_pose, pdbname );
		++counter;
	}
	TR(t_info) << std::endl;
}

void ProteinInterfaceMultiStateDesignMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	// flags/parameters (override options settings)
	if ( tag->hasOption("generations") ) generations_ = tag->getOption<Size>("generations");
	if ( tag->hasOption("pop_size") ) pop_size_ = tag->getOption<Size>("pop_size");
	if ( tag->hasOption("num_packs") ) num_packs_ = tag->getOption<Size>("num_packs");
	if ( tag->hasOption("pop_from_ss") ) pop_from_ss_ = tag->getOption<Size>("pop_from_ss");
	if ( tag->hasOption("numresults") ) numresults_ = tag->getOption<Size>("numresults");
	if ( tag->hasOption("fraction_by_recombination") ) {
		fraction_by_recombination_ = tag->getOption<Real>("fraction_by_recombination");
	}
	if ( tag->hasOption("mutate_rate") ) mutate_rate_ = tag->getOption<Real>("mutate_rate");
	if ( tag->hasOption("boltz_temp") ) boltz_temp_ = tag->getOption<Real>("boltz_temp");
	if ( tag->hasOption("anchor_offset") ) anchor_offset_ = tag->getOption<Real>("anchor_offset");
	// checkpointing options
	if ( tag->hasOption("checkpoint_prefix") ) {
		checkpoint_prefix_ = tag->getOption<std::string>("checkpoint_prefix");
	}
	if ( tag->hasOption("checkpoint_interval") ) {
		checkpoint_interval_ = tag->getOption<Size>("checkpoint_interval");
	}
	if ( tag->hasOption("checkpoint_gz") ) checkpoint_gz_ = tag->getOption<bool>("checkpoint_gz");
	if ( tag->hasOption("checkpoint_rename") ) {
		checkpoint_rename_ = tag->getOption<bool>("checkpoint_rename");
	}

	fname_prefix_ = tag->getOption< std::string >( "output_fname_prefix", "" );
	// calls to PackRotamersMover base class methods
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );

	TaskFactoryCOP tf = protocols::rosetta_scripts::parse_task_operations( tag, datamap );
	if ( tf ) task_factory( tf );
	rb_jump_ = tag->getOption< core::Size >( "rb_jump", 1 );

	unfolded_ = tag->getOption< bool >( "unfolded", 1 );
	unbound_ = tag->getOption< bool >( "unbound", 1 );
	input_is_positive_ = tag->getOption< bool >( "input_is_positive", 1 );
	use_unbound_for_sequence_profile_ = tag->getOption< bool >( "unbound_for_sequence_profile", unbound_ );
	bump_threshold_ = tag->getOption< core::Real >( "profile_bump_threshold", 1.0 );
	/// Read additional positive and negative states
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	bool at_least_one_negative_state( unfolded_ || unbound_ );
	compare_energy_to_ground_state_ = tag->getOption< bool >( "compare_to_ground_state", branch_tags.size() );
	TR<<"Compare energy to ground state set to: "<<compare_energy_to_ground_state_<<std::endl;
	for ( TagCOP const btag : branch_tags ) {
		if ( unfolded_ || unbound_ ) {
			TR<<"ERROR: If you specify additional pdb files as states, it is assumed that those would have different energies than the starting pdb. As such, comparison of energies across different states is automatically done by grounding each pdb file to its starting 'best-score design' and comparing energy differences from that state. The energies of unbound and unfolded states then become tricky to interpret. You can use anchor_offset to get much of the effect of these additional states. Or, ask Sarel."<<std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("");
		}
		std::string const fname( btag->getOption< std::string >( "pdb" ) );
		bool const unbound( btag->getOption< bool >( "unbound", 0 ) );
		bool const unfolded( btag->getOption< bool >( "unfolded", 0 ) );

		core::pose::PoseOP new_pose( new core::pose::Pose );
		state_poses_.push_back( new_pose );
		core::import_pose::pose_from_file( *new_pose, fname , core::import_pose::PDB_file);
		saved_state_poses_.push_back( core::pose::PoseOP( new core::pose::Pose( *new_pose ) ) ); //deep copying new pose so that its saved throughout the run
		state_unbound_.push_back( unbound );
		state_unfolded_.push_back( unfolded );

		TaskFactoryCOP state_tf = protocols::rosetta_scripts::parse_task_operations( btag, datamap );
		state_task_factory_.push_back( state_tf );

		if ( btag->getName() == "Positive" ) {
			state_positive_.push_back( true );
		} else if ( btag->getName() == "Negative" ) {
			state_positive_.push_back( false );
			at_least_one_negative_state = true;
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption( "Name "+btag->getName()+" is not recognized in ProteinInterfaceMultistateDesign::parse_my_tag." );
		}
		compare_energy_to_ground_state_ = true;
	}
	runtime_assert( at_least_one_negative_state );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
ProteinInterfaceMultiStateDesignMover::fresh_instance() const
{
	return moves::MoverOP( new ProteinInterfaceMultiStateDesignMover );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
ProteinInterfaceMultiStateDesignMover::clone() const
{
	return moves::MoverOP( new ProteinInterfaceMultiStateDesignMover( *this ) );
}

/// @details we build one target (bound) and two competitor (unbound and unfolded) states.
void
ProteinInterfaceMultiStateDesignMover::add_states(
	Pose const & pose
)
{
	using namespace core::pose;
	using namespace protocols::multistate_design;

	runtime_assert( multistate_packer_ != 0 );

	runtime_assert( task_factory() != 0 );
	PackerTaskCOP ptask = task_factory()->create_task_and_apply_taskoperations( pose );

	Pose const bound( pose );
	PoseOP unbound( new core::pose::Pose( pose ) );
	PoseOP unfolded( new core::pose::Pose( pose ) );

	if ( unbound_ ) {
		part_complex( unbound, rb_jump_ );
		part_complex( unfolded, rb_jump_ );
	}
	if ( unfolded_ ) {
		unfold_complex( unfolded );
	}
	PackingStateOP bound_state( new PackingState( bound, input_is_positive_ ) );
	PackingStateOP unbound_state( new PackingState( *unbound, false ) );
	PackingStateOP unfolded_state( new PackingState( *unfolded, false ) );

	bound_state->create_packer_data( scorefxn_, ptask );
	/// Sharing data is only appropriate if the interaction graph is identical
	/// between the states (e.g., in specificity calculations), but that is
	/// clearly not the case here.
	// unfolded_state->share_packer_data_from( *bound_state );
	// unbound_state->share_packer_data_from( *bound_state );
	// unfolded_state->create_packer_data( scorefxn_, ptask );
	unbound_state->create_packer_data( scorefxn_, ptask );
	unfolded_state->create_packer_data( scorefxn_, ptask );

	multistate_packer_->add_state( bound_state );
	if ( unbound_ ) {
		multistate_packer_->add_state( unbound_state );
	}
	if ( unfolded_ ) {
		multistate_packer_->add_state( unfolded_state );
	}

	runtime_assert( state_unbound_.size() == state_unfolded_.size() );
	runtime_assert( state_unbound_.size() == state_positive_.size() );
	runtime_assert( state_unbound_.size() == state_poses_.size() );
	runtime_assert( state_unbound_.size() == state_task_factory_.size() );

	for ( core::Size i( 1 ); i<=state_poses_.size(); ++i ) {
		if ( state_unbound_[ i ] || state_unfolded_[ i ] ) {
			part_complex( state_poses_[ i ], rb_jump_ );
			part_complex( saved_state_poses_[ i ], rb_jump_ );
		}
		if ( state_unfolded_[ i ] ) {
			unfold_complex( state_poses_[ i ] );
			unfold_complex( saved_state_poses_[ i ] );
		}
		PackingStateOP state( new PackingState( *state_poses_[ i ], state_positive_[ i ] ) );
		PackerTaskOP state_ptask = task_factory()->create_task_and_apply_taskoperations( *state_poses_[ i ] );
		if ( state_task_factory_[ i ] ) {
			state_task_factory_[ i ]->modify_task( pose, state_ptask );
		}
		state->create_packer_data( scorefxn_, state_ptask );
		multistate_packer_->add_state( state );
	}
}

std::string ProteinInterfaceMultiStateDesignMover::get_name() const {
	return mover_name();
}

std::string ProteinInterfaceMultiStateDesignMover::mover_name() {
	return "ProteinInterfaceMS";
}

std::string subtag_for_msd( std::string const & name ) {
	return "stfmsd_" + name + "_type";
}

void ProteinInterfaceMultiStateDesignMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "generations", xsct_non_negative_integer, "Number of population generations" )
		+ XMLSchemaAttribute( "pop_size", xsct_non_negative_integer, "Size of population" )
		+ XMLSchemaAttribute( "num_packs", xsct_non_negative_integer, "Number of packing rounds" )
		+ XMLSchemaAttribute( "pop_from_ss", xsct_non_negative_integer, "Number of population generated from mutating single-state seeds" )
		+ XMLSchemaAttribute( "numresults", xsct_non_negative_integer, "Number of results to preserve" )
		+ XMLSchemaAttribute( "fraction_by_recombination", xsct_real, "Fraction of the population generated by recombination of the parent generation" )
		+ XMLSchemaAttribute( "mutate_rate", xsct_real, "Mutation rate for simulation" )
		+ XMLSchemaAttribute( "boltz_temp", xsct_real, "Temperature of Monte Carlo simulation" )
		+ XMLSchemaAttribute( "anchor_offset", xsct_real, "Anchoring offset for faux alteranative states" )
		+ XMLSchemaAttribute( "checkpoint_prefix", xs_string, "Prefix for checkpointing files" )
		+ XMLSchemaAttribute( "checkpoint_interval", xsct_non_negative_integer, "Frequency of checkpointing" )
		+ XMLSchemaAttribute( "checkpoint_gz", xsct_rosetta_bool, "Checkpoint as gzipped files" )
		+ XMLSchemaAttribute( "checkpoint_rename", xsct_rosetta_bool, "Rename checkpoint files" )
		+ XMLSchemaAttribute( "output_fname_prefix", xs_string, "Prefix for output filename" );

	rosetta_scripts::attributes_for_parse_score_function( attlist );
	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "rb_jump", xsct_non_negative_integer, "Rigid-body jump that defines the interface, numbered from 1" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "unfolded", xsct_rosetta_bool, "Are we providing an unfolded alternative state?", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "unbound", xsct_rosetta_bool, "Are we providing an unbound alternative state?", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "input_is_positive", xsct_rosetta_bool, "Ought the input state to be subject to positive selection?", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "profile_bump_threshold", xsct_real, "The energy above which a bump check fails", "1.0" )
		+ XMLSchemaAttribute( "compare_to_ground_state", xsct_rosetta_bool, "Do we have a ground state comparison? Should be true, by default, if there are subelements" );

	AttributeList subtag_attributes;

	subtag_attributes + XMLSchemaAttribute::required_attribute( "pdb", xs_string, "File name containing this alternative state" )
		+ XMLSchemaAttribute::attribute_w_default( "unbound", xsct_rosetta_bool, "This alternative state is unbound", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "unfolded", xsct_rosetta_bool, "This alternative state is unfolded", "0" );

	rosetta_scripts::attributes_for_parse_task_operations( subtag_attributes );

	XMLSchemaSimpleSubelementList ssl;

	ssl.add_simple_subelement( "Positive", subtag_attributes, "Tags describing individual alternative states in the calculation"/*, 0 minoccurs*/ )
		.add_simple_subelement( "Negative", subtag_attributes, "Tags describing individual alternative states in the calculation"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_msd );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string ProteinInterfaceMultiStateDesignMoverCreator::keyname() const {
	return ProteinInterfaceMultiStateDesignMover::mover_name();
}

protocols::moves::MoverOP
ProteinInterfaceMultiStateDesignMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ProteinInterfaceMultiStateDesignMover );
}

void ProteinInterfaceMultiStateDesignMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinInterfaceMultiStateDesignMover::provide_xml_schema( xsd );
}



} // namespace movers
} // namespace protein_interface_design
} // namespace protocols
