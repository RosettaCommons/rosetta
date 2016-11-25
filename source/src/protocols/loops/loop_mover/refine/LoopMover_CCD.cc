// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/LoopMover_Perturb_CCD.cc
/// @brief kinematic loop closure main protocols
/// @author Chu Wang
/// @author Mike Tyka

// Unit headers
#include <protocols/loops/loops_main.hh>

// Package headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleFactory.hh>

#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCDCreator.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/loops_definers/util.hh>
#include <protocols/moves/MonteCarlo.hh>

// Rosetta Headers

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>


#include <core/pose/symmetry/util.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <protocols/loops/util.hh>
#include <basic/Tracer.hh> // tracer output

//Utility and numeric Headers
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// option key includes

#include <basic/options/keys/MonteCarlo.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using namespace ObjexxFCL::format;

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////
using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_mover.refine.LoopMover_Refine_CCD" );

Real const CHAINBREAK_SCORE_RAMP_FACTOR = 10./3.;

//constructors
LoopMover_Refine_CCD::LoopMover_Refine_CCD()
: LoopMover(),
	set_fold_tree_from_loops_(false),
	user_defined_move_map_(false),
	debug_(false),
	outer_cycles_(3),
	max_inner_cycles_(200),
	repack_period_(20),
	temp_initial_(1.5),
	temp_final_(0.5)
{
	read_options();
	set_scorefxn( get_fa_scorefxn() );
	protocols::moves::Mover::type("LoopMover_Refine_CCD");
	set_default_settings();
}

LoopMover_Refine_CCD::LoopMover_Refine_CCD(
	protocols::loops::LoopsOP  loops_in
) : LoopMover( loops_in ),
	set_fold_tree_from_loops_(false),
	user_defined_move_map_(false),
	debug_(false),
	outer_cycles_(3),
	max_inner_cycles_(200),
	repack_period_(20),
	temp_initial_(1.5),
	temp_final_(0.5)
{
	read_options();
	set_scorefxn( get_fa_scorefxn() );
	protocols::moves::Mover::type("LoopMover_Refine_CCD");
	set_default_settings();
}


LoopMover_Refine_CCD::LoopMover_Refine_CCD(
	protocols::loops::LoopsOP  loops_in,
	core::scoring::ScoreFunctionOP  scorefxn
) : LoopMover( loops_in ),
	set_fold_tree_from_loops_(false),
	user_defined_move_map_(false),
	debug_(false),
	outer_cycles_(3),
	max_inner_cycles_(200),
	repack_period_(20),
	temp_initial_(1.5),
	temp_final_(0.5)
{
	read_options();
	set_scorefxn( scorefxn );
	protocols::moves::Mover::type("LoopMover_Refine_CCD");
	set_default_settings();
}

//destructor
LoopMover_Refine_CCD::~LoopMover_Refine_CCD() {}

// XRW TEMP std::string
// XRW TEMP LoopMover_Refine_CCD::get_name() const {
// XRW TEMP  return "LoopMover_Refine_CCD";
// XRW TEMP }

void
LoopMover_Refine_CCD::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Loops:\n" << get_loops();
	output <<   "Outer cycles:        " << outer_cycles_ << "\nMax inner cycles:    " << max_inner_cycles_ <<
		"\nRepack period:       " << repack_period_ << "\nInitial temperature: " << temp_initial_ <<
		"\nFinal temperature:   " << temp_final_ <<  "\nSet fold tree from loop?: " <<
		(set_fold_tree_from_loops_ ? "True" : "False") << "\nMovemap:  ";
	if ( move_map() != 0 ) { output << std::endl; move_map()->show(output);}
	else { output << "none" << std::endl; }
}

//clone
protocols::moves::MoverOP LoopMover_Refine_CCD::clone() const {
	return protocols::moves::MoverOP( new LoopMover_Refine_CCD(*this) );
}


void LoopMover_Refine_CCD::set_default_settings()
{
	redesign_loop_ = false;
	packing_isolated_to_active_loops_ = false;
	flank_residue_min_ = false; // added by JQX
	move_map_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	inner_cycles_ = max_inner_cycles_;
	current_cycle_number_ = 0;
	original_fold_tree_ = NULL;
}
loop_mover::LoopMover::MoveMapOP LoopMover_Refine_CCD::move_map() const
{
	return move_map_;
}
void LoopMover_Refine_CCD::move_map( LoopMover::MoveMapOP mm )
{
	if ( mm ) {
		move_map_ = mm;
		user_defined_move_map_ = true;
	} else {
		user_defined_move_map_ = false;
	}
}

protocols::loops::LoopsCOP LoopMover_Refine_CCD::get_loops() const {
	return loops();
}

std::ostream &operator<< ( std::ostream &os, LoopMover_Refine_CCD const &mover )
{
	mover.show(os);
	return os;
}

void
LoopMover_Refine_CCD::read_options()
{
	using namespace basic::options;
	outer_cycles_ = option[ OptionKeys::loops::refine_outer_cycles ]();
	if ( option[ OptionKeys::loops::max_inner_cycles ].user() ) {
		max_inner_cycles_ = option[ OptionKeys::loops::max_inner_cycles ]();
	}
	if ( option[ OptionKeys::loops::repack_period ].user() ) {
		repack_period_ = option[ OptionKeys::loops::repack_period ]();
	}
	if ( option[ OptionKeys::MonteCarlo::temp_initial ].user() ) {
		temp_initial_ = option[ OptionKeys::MonteCarlo::temp_initial ]();
	}
	if ( option[ OptionKeys::MonteCarlo::temp_final ].user() ) {
		temp_final_ = option[ OptionKeys::MonteCarlo::temp_final ]();
	}

	debug_ = option[ OptionKeys::loops::debug ].user();

	repack_neighbors_ = (! basic::options::option[ basic::options::OptionKeys::loops::fix_natsc ]);
}

void LoopMover_Refine_CCD::set_task_factory(
	core::pack::task::TaskFactoryCOP task_factory_in
)
{
	// make local, non-const copy from const input
	runtime_assert( task_factory_in != 0 );
	task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory( *task_factory_in ) );
}

core::pack::task::TaskFactoryCOP LoopMover_Refine_CCD::get_task_factory() const { return task_factory_; }


void
LoopMover_Refine_CCD::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose ){
	using namespace basic::options;
	packing_isolated_to_active_loops_ = false;
	//using parser implies that the fold tree probably isn't set correctly
	set_fold_tree_from_loops( tag->getOption< bool >( "set_fold_tree_from_loops", true ) );
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	bool specified_movemap( false );
	for ( utility::tag::TagCOP tag : branch_tags ) {
		if ( tag->getName() == "MoveMap" ) specified_movemap = true;
		break;
	}
	if ( specified_movemap ) {
		move_map( LoopMover::MoveMapOP( new core::kinematics::MoveMap ) );
		move_map_->set_bb( false );
		move_map_->set_chi( false );
		move_map_->set_jump( false );
		protocols::rosetta_scripts::parse_movemap( tag, pose, move_map_, data, false/*don't reset movemap, keep falses, unless stated otherwise*/ );
	}
	if ( tag->hasOption( "loops" ) ) {
		loops( loops_definers::load_loop_definitions(tag, data, pose) );
	}
	if ( tag->hasOption( "scorefxn" ) ) this->set_scorefxn( data.get< core::scoring::ScoreFunction * >( "scorefxns", tag->getOption<std::string>( "scorefxn" ) )->clone() );

	if ( tag->hasOption("task_operations") ) {
		core::pack::task::TaskFactoryOP task_factory = protocols::rosetta_scripts::parse_task_operations( tag, data );
		this->set_task_factory( task_factory );
	} else task_factory_ = NULL;

	if ( tag->hasOption( "loops_from_cache" ) ) set_use_loops_from_observer_cache( tag->getOption<bool>( "loops_from_cache", 1 ) );

	if ( tag->hasOption( "outer_cycles" ) ) outer_cycles_ = tag->getOption<core::Size>( "outer_cycles", option[ OptionKeys::loops::refine_outer_cycles ]() );
	if ( tag->hasOption( "max_inner_cycles" ) ) max_inner_cycles_ = tag->getOption<core::Size>( "max_inner_cycles", 250 );
	temp_initial( tag->getOption< core::Real >( "temp_initial", 1.5 ) );
	temp_final( tag->getOption< core::Real >( "temp_final", 0.5 ) );
}

void LoopMover_Refine_CCD::apply( core::pose::Pose & pose )
{
	using namespace scoring;
	using namespace basic::options;

	// scheduler
	int const fast( option[ OptionKeys::loops::fast ] ); // why is this an int?
	inner_cycles_ = std::min( max_inner_cycles_, fast ? loops()->loop_size() : Size(10)*loops()->loop_size() );

	/// must be called once the Pose has become available.
	resolve_loop_indices( pose );

	if ( ! get_native_pose() ) set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ) );
	if ( use_loops_from_observer_cache() ) this->set_loops_from_pose_observer_cache( pose );

	setup_foldtree_and_add_cutpoint_variants( pose );

	// Use 'get_new_ramping_scorefxn' in case the result of scorefunction() has changed since the last apply. Set weights.
	get_new_ramping_scorefxn()->set_weight( chainbreak, (1. * CHAINBREAK_SCORE_RAMP_FACTOR) );

	Real const temperature_ramp_factor = std::pow( (temp_final_ / temp_initial_),
		Real(1.0f / (outer_cycles_ * inner_cycles_)) );

	// make sure we have scored before we instantiate monte carlo and before we ask for a tenA neighbor graph
	(*ramping_scorefxn())(pose);

	protocols::moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *ramping_scorefxn(), temp_initial_ ) );
	pack::task::PackerTaskOP pack_task = get_packer_task( pose );
	pack::pack_rotamers( pose, *ramping_scorefxn(), pack_task );
	std::string move_type = "repack";
	mc->boltzmann( pose, move_type );
	mc->show_scores();

	debugging_output( pose );

	// set minimization degrees of freedom for all loops
	setup_movemap( pose, *loops(), pack_task->repacking_residues(), move_map_ );

	LoopRefineInnerCycleOP inner_cycle = LoopRefineInnerCycleFactory::get_instance()->create_inner_cycle(
		IC_RefineCCDStandard,
		LoopMover_Refine_CCDAP( utility::pointer::static_pointer_cast< LoopMover_Refine_CCD >(get_self_ptr()) ),
		mc,
		ramping_scorefxn(),
		task_factory_
	);

	inner_cycle->set_native_pose( get_native_pose() ); // Native pose may be used in debugging steps that use loop RMSD

	for ( Size i = 1; i <= outer_cycles_; ++i ) {
		increase_chainbreak_weight_and_update_monte_carlo( i, ramping_scorefxn(), *mc, pose );
		for ( current_cycle_number_ = 1; current_cycle_number_ <= inner_cycles_; ++current_cycle_number_ ) {
			mc->set_temperature( mc->temperature() * temperature_ramp_factor );
			tr().Info << "refinement cycle (outer/inner): " << i << "/" << outer_cycles_ << " ";
			tr().Info << current_cycle_number_ << "/" << inner_cycles_ << " " << std::endl;

			inner_cycle->apply( pose );

			tr().Info << std::flush;
		} //inner_cycle
	} //outer_cycle

	mc->show_counters();
	pose = mc->lowest_score_pose();

	if ( set_fold_tree_from_loops_ ) { //if requested, put back old foldtree
		loops::remove_cutpoint_variants( pose );
		pose.fold_tree( *original_fold_tree_ );
		(*ramping_scorefxn())(pose);
	}
}

/// @brief setup an appropriate movemap for the given loops
/// @param[in] loops The loops to model.
/// @param[in] allow_repack Indicates whether or not to allow a position to
///  repack.
/// @param[out] movemap Output movemap, all settings added here.
/// @remarks will enforce the false movemap
void LoopMover_Refine_CCD::setup_movemap(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	utility::vector1< bool > const & allow_repack,
	core::kinematics::MoveMapOP & movemap
)
{
	if ( user_defined_move_map_ ) {
		movemap = move_map_;
		return;
	}
	loops_set_move_map( loops, allow_repack, *movemap );
	enforce_false_movemap( movemap );
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	}
}

basic::Tracer & LoopMover_Refine_CCD::tr() const
{
	return TR;
}

// XRW TEMP LoopMover_Refine_CCDCreator::~LoopMover_Refine_CCDCreator() {}

// XRW TEMP moves::MoverOP LoopMover_Refine_CCDCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new LoopMover_Refine_CCD() );
// XRW TEMP }

// XRW TEMP std::string LoopMover_Refine_CCDCreator::keyname() const {
// XRW TEMP  return "LoopMover_Refine_CCD";
// XRW TEMP }

core::scoring::ScoreFunctionOP LoopMover_Refine_CCD::get_new_ramping_scorefxn()
{
	if ( scorefxn() != 0 ) {
		ramping_scorefxn_ = scorefxn()->clone();
	} else {
		ramping_scorefxn_ = get_fa_scorefxn();
	}
	return ramping_scorefxn_;
}

core::scoring::ScoreFunctionOP LoopMover_Refine_CCD::ramping_scorefxn()
{
	if ( ! ramping_scorefxn_ ) {
		return get_new_ramping_scorefxn();
	}
	return ramping_scorefxn_;
}

void LoopMover_Refine_CCD::setup_foldtree_and_add_cutpoint_variants( core::pose::Pose & pose )
{
	if ( set_fold_tree_from_loops_ ) {
		core::kinematics::FoldTree f_new;
		original_fold_tree_ = core::kinematics::FoldTreeOP( new core::kinematics::FoldTree( pose.fold_tree() ) );
		loops::fold_tree_from_loops( pose, *( this->loops() ), f_new);
		pose.fold_tree( f_new );
	}
	loops::add_cutpoint_variants( pose );
}

core::pack::task::PackerTaskOP LoopMover_Refine_CCD::get_packer_task( core::pose::Pose const & pose )
{
	pack::task::PackerTaskOP base_packer_task;
	// create default Packer behavior if none has been set via TaskFactory
	if ( task_factory_ == 0 ) {
		// the default Packer behavior is defined here
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using toolbox::task_operations::RestrictToLoopsAndNeighbors;
		using toolbox::task_operations::RestrictToLoopsAndNeighborsOP;
		task_factory_ = core::pack::task::TaskFactoryOP( new TaskFactory );
		task_factory_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
		task_factory_->push_back( TaskOperationCOP( new IncludeCurrent ) );
		task_factory_->push_back( TaskOperationCOP( new NoRepackDisulfides ) );

		RestrictToLoopsAndNeighborsOP restrict_to_loops_and_neighbors( new RestrictToLoopsAndNeighbors() );

		// This can be simplified by making a constructor that takes these settings as arguments
		restrict_to_loops_and_neighbors->set_cutoff_distance( 10.0 );
		restrict_to_loops_and_neighbors->set_design_loop( redesign_loop_ );
		restrict_to_loops_and_neighbors->set_include_neighbors( repack_neighbors_ );
		restrict_to_loops_and_neighbors->set_loops( loops() );

		task_factory_->push_back( restrict_to_loops_and_neighbors );

		// additional default behavior: packing restricted to active loops
		packing_isolated_to_active_loops_ = true;
	}
	base_packer_task = task_factory_->create_task_and_apply_taskoperations( pose );
	base_packer_task->set_bump_check( true );
	return base_packer_task;
}

void LoopMover_Refine_CCD::increase_chainbreak_weight_and_update_monte_carlo(
	Size iteration_number,
	scoring::ScoreFunctionOP local_scorefxn,
	protocols::moves::MonteCarlo & mc,
	pose::Pose & pose
) {
	// increase CHAINBREAK weight and update monte carlo
	local_scorefxn->set_weight( scoring::chainbreak, Real(iteration_number)* CHAINBREAK_SCORE_RAMP_FACTOR );
	mc.score_function( *local_scorefxn );
	// recover low
	mc.recover_low( pose );
	// score info

	tr().Info << "cycle: " << iteration_number << "  " << (*local_scorefxn)(pose) << std::endl;
}

void LoopMover_Refine_CCD::debugging_output( core::pose::Pose & pose )
{
	if ( debug_ ) {
		pose.dump_pdb("tmp_fa_repack.pdb");
		std::ofstream out("score.tmp_repack_fa");
		out << "scoring for repack_fa " << (*ramping_scorefxn())(pose) << std::endl;
		ramping_scorefxn()->show( out );
		out << pose.energies();
	}
}

std::string LoopMover_Refine_CCD::get_name() const {
	return mover_name();
}

std::string LoopMover_Refine_CCD::mover_name() {
	return "LoopMover_Refine_CCD";
}

void LoopMover_Refine_CCD::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("set_fold_tree_from_loops", xsct_rosetta_bool,
		"Set fold tree using loop info. Default = true","true" )
		+ XMLSchemaAttribute("scorefxn", xs_string, "Set score function to be used." )
		+ XMLSchemaAttribute::attribute_w_default("loops_from_cache", xsct_rosetta_bool,
		"Allow the loops to be set from the segments stored in the poses observer cache."
		"Makes it possible to have LoopMovers be part of parser protocols where the loops"
		"were determined by some previous on the fly step . Default = true","1" )
		+ XMLSchemaAttribute("outer_cycles", xsct_non_negative_integer, "Set number of outer cycles." )
		+ XMLSchemaAttribute::attribute_w_default("max_inner_cycles", xsct_non_negative_integer,
		"Finish inner cycles latest after n cycles. Default = 250", "250" )
		+ XMLSchemaAttribute::attribute_w_default("temp_initial", xsct_real, "Initial Boltzman temperature. Default = 1.5", "1.5" )
		+ XMLSchemaAttribute::attribute_w_default("temp_final", xsct_real, "Final Boltzman temperature. Default = 0.5", "0.5" );
	loops_definers::attributes_for_load_loop_definitions( attlist );
	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	AttributeList subelement_attributes;
	XMLSchemaSimpleSubelementList subelement_list;
	rosetta_scripts::append_subelement_for_parse_movemap_w_datamap( xsd, subelement_list );
	subelement_list.add_simple_subelement("MoveMap", subelement_attributes, "Define the move map.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Perform loop refinement using CCD for loop closure.", attlist, subelement_list );
}

std::string LoopMover_Refine_CCDCreator::keyname() const {
	return LoopMover_Refine_CCD::mover_name();
}

protocols::moves::MoverOP
LoopMover_Refine_CCDCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopMover_Refine_CCD );
}

void LoopMover_Refine_CCDCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopMover_Refine_CCD::provide_xml_schema( xsd );
}


} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
