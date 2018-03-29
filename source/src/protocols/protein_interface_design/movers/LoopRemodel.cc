// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/LoopRemodel.cc
/// @brief Parseable class to run loop perturbation or refinement between a given loop between start/end (inclusive). Intended to sample loop conformations at an interface.
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopRemodel.hh>
#include <protocols/protein_interface_design/movers/LoopRemodelCreator.hh>

// Package headers

// Project headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loops_main.hh> // for various loop utility fxns
#include <protocols/loops/Loops.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/methods/fragment_functions.hh> // smallmer_from_largemer
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <protocols/minimization_packing/MinMover.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/dssp/Dssp.hh> // for getting secondary structure for fragments

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

#include <core/fragment/Frame.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/Jump.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <protocols/calc_taskop_movers/DesignRepackMover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.LoopRemodel" );

// XRW TEMP std::string
// XRW TEMP LoopRemodelCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LoopRemodel::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoopRemodelCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoopRemodel );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LoopRemodel::mover_name()
// XRW TEMP {
// XRW TEMP  return "LoopRemodel";
// XRW TEMP }

LoopRemodel::~LoopRemodel() = default;

protocols::moves::MoverOP
LoopRemodel::clone() const {
	return( protocols::moves::MoverOP( new LoopRemodel( *this ) ) );
}

LoopRemodel::LoopRemodel() :
	calc_taskop_movers::DesignRepackMover( LoopRemodel::mover_name() )
{}

LoopRemodel::LoopRemodel(
	std::string const & protocol,
	std::string const & loop_start,
	std::string const & loop_end,
	core::Size const cycles,
	bool const perturb,
	bool const refine,
	bool const hurry,
	core::scoring::ScoreFunctionOP hires_score,
	core::scoring::ScoreFunctionOP lores_score,
	protocols::loops::LoopsCOP loops,
	core::fragment::FragSetOP frag1,
	core::fragment::FragSetOP frag3,
	core::fragment::FragSetOP frag9
) :
	calc_taskop_movers::DesignRepackMover( LoopRemodel::mover_name() ),
	protocol_( protocol ),
	loop_start_( loop_start ),
	loop_end_( loop_end ),
	cycles_( cycles ),
	perturb_( perturb),
	refine_( refine),
	hurry_( hurry )
{
	hires_score_ = hires_score;
	lores_score_ = lores_score->clone();
	if ( loops ) {
		loops_ = protocols::loops::LoopsOP( new protocols::loops::Loops( *loops ) );
	}
	frag1_ = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( *frag1 ) );
	frag3_ = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( *frag3 ) );
	frag9_ = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( *frag9 ) );
}

void
LoopRemodel::apply( core::pose::Pose & pose )
{
	using namespace protocols::loops;
	using core::pack::task::operation::TaskOperationCOP;

	core::pose::Pose native_pose = pose;

	LoopsOP loops;
	if ( loops_ ) {
		loops = LoopsOP( new protocols::loops::Loops( *loops_ ) ); // Make a copy
	} else if ( !loop_start_.empty() && !loop_end_.empty() ) {
		core::Size loop_start = core::pose::parse_resnum( loop_start_, pose );
		core::Size loop_end   = core::pose::parse_resnum( loop_end_, pose );
		runtime_assert( loop_end > loop_start );
		runtime_assert( (loop_end - loop_start) >= 3 );
		runtime_assert( loop_start > 1 );
		runtime_assert( loop_end < pose.size() );
		core::Size const cutpt = (loop_start+loop_end)/2; // put cutpoint in the middle of the loop
		loops = protocols::loops::LoopsOP( new protocols::loops::Loops );
		loops->add_loop( protocols::loops::Loop( loop_start, loop_end, cutpt ) );
	} else {
		utility_exit_with_message("In LoopRemodel mover - must specify either a Loops definition or a loop start and end.");
	}

	if ( loops->size() == 0 )  {
		TR << "No loops found!" << std::endl;
		return; // bounce out if we didn't define any loops
	} else TR << *loops << std::endl;

	if ( loops->size() > 0 ) {
		// set up temporary fold tree for loop closure
		TR.Debug << "Original FoldTree " << pose.fold_tree() << std::endl;
		core::kinematics::FoldTree old_ft( pose.fold_tree() );
		//core::kinematics::FoldTree new_ft;
		core::kinematics::FoldTree new_ft = protocols::forge::methods::fold_tree_from_loops( pose, *loops );
		//fold_tree_from_loops( pose, *loops, new_ft, false /*terminal_cutpoint*/ );
		pose.fold_tree( new_ft );
		add_cutpoint_variants( pose );
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent ) );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::NoRepackDisulfides ) );
		// perturbation and closure
		if ( hurry_ ) {
			core::Real mc_kt( 1.0 );
			if ( basic::options::option[ basic::options::OptionKeys::loops::remodel_init_temp ].user() ) mc_kt = basic::options::option[ basic::options::OptionKeys::loops::remodel_init_temp ]();
			core::Size const outer_cycles( cycles_ );
			protocols::moves::MonteCarlo outer_mc( pose, *hires_score_, mc_kt );
			outer_mc.set_autotemp( true, mc_kt );
			for ( core::Size i = 1; i <= outer_cycles; ++i ) {
				TR.Debug << "outer_cycle " << i << " kt=" << outer_mc.temperature() << std::endl;
				core::Size const inner_cycles( 20 );
				protocols::moves::MonteCarlo inner_mc( pose, *hires_score_, mc_kt );
				inner_mc.set_autotemp( true, mc_kt );
				for ( core::Size j = 1; j <= inner_cycles; ++j ) {
					TR.Debug << "inner_cycle " << j << " kt=" << inner_mc.temperature() << std::endl;
					for ( auto it = loops->v_begin(); it != loops->v_end(); ++it ) {
						// make a temporary loop/loops set to use in this scope
						Loop loop( *it );

						loops::loop_closure::kinematic_closure::KinematicMoverOP kinmover( new loops::loop_closure::kinematic_closure::KinematicMover );
						if ( perturb_ ) kinmover->set_idealize_loop_first( true );
						core::Size const cycles = 100;
						kinmover->set_temperature( mc_kt );
						kinmover->set_sfxn(hires_score_);
						kinmover->set_vary_bondangles( true );
						kinmover->set_sample_nonpivot_torsions( true );
						kinmover->set_rama_check( true );

						protocols::loops::loop_closure::kinematic_closure::KinematicWrapper kinwrapper( kinmover, loop, cycles );
						kinwrapper.apply( pose );
					} // for all loops
					inner_mc.boltzmann( pose );
				} // for inner_cycles
				inner_mc.show_counters();
				if ( basic::options::option[ basic::options::OptionKeys::loops::kic_recover_last ].value() ) {
					pose = inner_mc.last_accepted_pose();
				} else inner_mc.recover_low( pose );

				//} // if perturb

				(*hires_score_)(pose); // score the pose (for safety & to get a good graph state)
				// repack to relieve clashes

				if ( prevent_repacking().size() ) {
					using namespace core::pack::task::operation;
					OperateOnCertainResiduesOP prevent_repacking_on_certain_res( new OperateOnCertainResidues );
					prevent_repacking_on_certain_res->residue_indices( prevent_repacking() );
					prevent_repacking_on_certain_res->op( ResLvlTaskOperationCOP( new PreventRepackingRLT ) );
					task_factory->push_back( prevent_repacking_on_certain_res );
				}

				if ( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ) {
					task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );
				}
				if ( !design() ) task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
				core::pack::task::PackerTaskOP task = task_factory->create_task_and_apply_taskoperations( pose );

				if ( design() ) {
					for ( auto it = loops->v_begin(); it != loops->v_end(); ++it ) {
						Loop loop( *it );
						for ( core::Size i = 1; i <= pose.size(); ++i ) {
							if ( i >= loop.start() && i <= loop.stop() ) continue; // design
							else task->nonconst_residue_task( i ).restrict_to_repacking(); // repack only
						}
					}
				}
				protocols::minimization_packing::PackRotamersMover pack( hires_score_, task );
				pack.apply( pose );

				if ( refine_ ) {
					loops_set_move_map( pose, *loops, refine_, *movemap ); // bb, except for omega and all sidechains
					core::scoring::ScoreFunctionOP copy_score( hires_score_->clone() );
					copy_score->set_weight( core::scoring::chainbreak, 10.0 ); // upweight chainbreak, to strongly disfavor breaks
					copy_score->set_weight( core::scoring::omega, 0.5 ); // omega term to keep backbone healthy
					protocols::minimization_packing::MinMoverOP minmover( new protocols::minimization_packing::MinMover( movemap, copy_score, "lbfgs_armijo_nonmonotone", 1e-5, true) ); // DJM has reported better results with dfpmin
					minmover->apply( pose );
				} // if refine
				outer_mc.boltzmann( pose );
			} // for outer_cycles
			outer_mc.recover_low( pose );
			outer_mc.show_counters();
		} else { // if hurry // !hurry
			// pose will always start full atom
			//protocols::moves::MonteCarlo mc( pose, *scorefxn_repack_, mc_kt );
			SaveAndRetrieveSidechains retrieve_sc( pose );
			retrieve_sc.allsc( true );
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
			if ( protocol_ == "kinematic" ) {
				if ( perturb_ ) {
					for ( auto it = loops->v_begin(); it != loops->v_end(); ++it ) {
						it->set_extended( true ); // set all loops to extended (needed for kinematic mover to really perturb)
					}
					protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC perturb( loops, lores_score_ );
					perturb.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					perturb.apply( pose );
				}
				core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
				retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose
				if ( refine_ ) {
					protocols::loops::loop_mover::refine::LoopMover_Refine_KIC refine( loops, hires_score_ );
					refine.set_redesign_loop( design() ); // design?
					//if( task_factory() ) refine.set_task_factory( task_factory() ); // if we have a task factory set, then we should pass it to the loop mover
					refine.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					pose.update_residue_neighbors();
					refine.apply( pose );
				}
			} else if ( protocol_ == "ccd" ) { // protocol == kinematic
				pose.update_residue_neighbors();
				core::scoring::dssp::Dssp dssp( pose );
				dssp.insert_ss_into_pose( pose );
				std::string const full_ss = pose.secstruct();
				std::string const full_sequence = pose.sequence();

				bool const pick_status = pick_loop_frags( loops, full_sequence, full_ss );
				if ( !pick_status ) {
					set_last_move_status( protocols::moves::FAIL_RETRY );
					return;
				}
				if ( perturb_ ) {
					//protocols::loops::LoopMover_Perturb_QuickCCD perturb(*loops, lores_score_ );
					for ( auto it = loops->v_begin(); it != loops->v_end(); ++it ) {
						it->set_extended( true ); // set all loops to extended (needed for CCD mover to really perturb)
					}
					protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCD perturb( loops, lores_score_ );
					perturb.add_fragments( frag1_ );
					perturb.add_fragments( frag3_ );
					perturb.add_fragments( frag9_ );
					//perturb.randomize_loop( true );
					perturb.set_strict_loops( true );
					perturb.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					perturb.apply( pose );
				}
				core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
				retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose
				if ( refine_ ) {
					protocols::loops::loop_mover::refine::LoopMover_Refine_CCD refine( loops, hires_score_ );
					refine.add_fragments( frag1_ );
					refine.add_fragments( frag3_ );
					refine.add_fragments( frag9_ );
					refine.set_redesign_loop( design() );
					//if( task_factory() ) refine.set_task_factory( task_factory() ); // if we have a task factory set, then we should pass it to the loop mover
					refine.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					refine.apply( pose );
				}
			} else if ( protocol_ == "remodel" ) { // protocol == ccd
				// remodel starts as fa
				core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
				retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose

				core::scoring::dssp::Dssp dssp( pose );
				dssp.insert_ss_into_pose( pose );
				std::string const full_ss = pose.secstruct();
				std::string const full_sequence = pose.sequence();
				bool const pick_status = pick_loop_frags( loops, full_sequence, full_ss );
				if ( !pick_status ) {
					set_last_move_status( protocols::moves::FAIL_RETRY );
					return;
				}
				protocols::forge::remodel::RemodelLoopMover remodel( loops );
				remodel.scorefunction( *hires_score_ );
				remodel.add_fragments( frag1_ );
				remodel.add_fragments( frag3_ );
				remodel.add_fragments( frag9_ );
				remodel.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
				for ( core::Size i = 1; i <= cycles_; ++i ) {
					remodel.apply( pose );
					if ( (remodel.get_last_move_status() == protocols::moves::MS_SUCCESS) || (remodel.get_last_move_status() == protocols::moves::FAIL_DO_NOT_RETRY) ) break;
				}
			}
		}

		// revert to the original FT and take out cutpoints
		remove_cutpoint_variants( pose, true );
		pose.fold_tree( old_ft );
		TR.Debug << "Reverted FoldTree " << pose.fold_tree() << std::endl;
	} else TR << "No loops found!" << std::endl;
}

// XRW TEMP std::string
// XRW TEMP LoopRemodel::get_name() const {
// XRW TEMP  return LoopRemodel::mover_name();
// XRW TEMP }

// true if all fragments picked.
// false if something went wrong
bool
LoopRemodel::pick_loop_frags( protocols::loops::LoopsCOP loops_in, std::string const & full_sequence, std::string const & full_ss )
{
	using namespace core;
	using namespace core::fragment;
	using namespace protocols::loops;

	LoopsCOP loops( LoopsOP( new Loops( *loops_in ) ) );
	for ( core::Size frag_length = 3; frag_length <= 9; frag_length+=6 ) { // frag3 and frag9
		TR << "Finding " << frag_length <<"mer loop fragments..." << std::endl;
		if ( frag_length == 3 ) frag3_ = core::fragment::FragSetOP( new ConstantLengthFragSet( frag_length ) );
		else if ( frag_length == 9 ) frag9_ = core::fragment::FragSetOP( new ConstantLengthFragSet( frag_length ) );

		for ( auto const & it : *loops ) {
			// make a temporary loop/loops set to use in this scope
			LoopCOP loop( LoopOP( new Loop(it ) ) );
			if ( loop->size() < frag_length ) continue; // fragment extends past loop

			for ( core::Size i=loop->start(); i <= loop->stop() - frag_length; ++i ) {

				// figure out ss and sequence that we're currently working on
				std::string ss = full_ss.substr( i - 1, frag_length); // subtract 1 to get string indexing -> pose indexing
				if ( ss.length() < frag_length ) {
					ss.append( frag_length - ss.length(), 'D' );
				}
				TR.Debug << "Window ss: " << ss << "\n";
				std::string aa = full_sequence.substr( i - 1, frag_length);
				if ( aa.length() < frag_length ) {
					aa.append( frag_length - aa.length(), '.' );
				}
				TR.Debug << "Window aa: " << aa << std::endl;

				// pick fragments
				FragDataOPs list;
				if ( design() ) list =  picking_old::vall::pick_fragments_by_ss( ss, 4000, true /*add random noise*/ ); //magic number: 4000 fragments
				else list = picking_old::vall::pick_fragments_by_ss_plus_aa( ss, aa, 4000, true );
				for ( FragDataOPs::const_iterator it = list.begin(); it != list.end(); ++it ) {
					TR.Debug << (*it)->size() << " " << (*it)->sequence() << std::endl;
				}

				// add frames to fragset
				TR.Debug << "Adding frame: "<< i << "-" << i+frag_length << ": " << ss << " " << aa << std::endl;
				core::fragment::FrameOP frame( new core::fragment::Frame( i, frag_length ) );
				frame->add_fragment( list );
				if ( frag_length == 3 ) frag3_->add( frame );
				else if ( frag_length == 9 ) frag9_->add( frame );
			} // for all residues in each loop
		} // for all loops
	} // frag3 and frag9
	if ( frag3_->size() ) {
		frag1_ = core::fragment::FragSetOP( new ConstantLengthFragSet( 1 ) );
		frag1_->add( *protocols::forge::methods::smallmer_from_largemer( frag3_->begin(), frag3_->end(), 1 ) );
	}

	// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
	picking_old::FragmentLibraryManager::get_instance()->clear_Vall();
	if ( (frag1_->size() > 0) || (frag3_->size() > 0) || (frag9_->size() > 0) ) return true;
	else return false;
}


void
LoopRemodel::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	protocol_ = tag->getOption<std::string>( "protocol", "ccd" );
	perturb_ = tag->getOption<bool>( "perturb", false );
	refine_ = tag->getOption<bool>( "refine", true );

	hires_score_ = protocols::rosetta_scripts::parse_score_function( tag, "refine_score", data )->clone();
	lores_score_ = protocols::rosetta_scripts::parse_score_function( tag, "perturb_score", data, "score4L" )->clone();

	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );

	bool auto_loops = tag->getOption<bool>( "auto_loops", false );
	bool const des = tag->getOption<bool>( "design", false );
	design( des ); // set baseclass design flag

	hurry_ = tag->getOption<bool>( "hurry", false );
	cycles_ = tag->getOption<Size>( "cycles", 10 );
	runtime_assert( cycles_ > 0 );

	// populate loops
	if ( auto_loops ) {
		if ( !data.has( "loops", "found_loops" ) ) {
			TR << "Loops not present in basic::datacache::DataMap! Be sure to add LoopFinder before LoopRemodel!" << std::endl;
			return;
		}
		loops_ = data.get_ptr<protocols::loops::Loops>( "loops", "found_loops" ); // from LoopFinder
	} else {
		loop_start_ = core::pose::get_resnum_string( tag, "loop_start_", "" );
		loop_end_   = core::pose::get_resnum_string( tag, "loop_end_", "" );
	}

	if ( !auto_loops && ( loop_start_.empty() || loop_end_.empty() ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "LoopRemodel must have either auto_loops set, or both loop_start and loop_end designations");
	}

	TR << "LoopRemodel mover: auto_loops="<<auto_loops<<" loop_start="<<loop_start_<<" loop_end="<<loop_end_<<" design="<<design()<<
		" perturb="<<perturb_<<" refine="<<refine_<<" hurry=" <<hurry_<< " cycles=" << cycles_ <<" hires_score="<< rosetta_scripts::get_score_function_name(tag, "refine_score") << std::endl;
}

std::string LoopRemodel::get_name() const {
	return mover_name();
}

std::string LoopRemodel::mover_name() {
	return "LoopRemodel";
}

void LoopRemodel::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction closuremethod;
	closuremethod.name( "closuremethod" );
	closuremethod.base_type( xs_string );
	closuremethod.add_restriction( xsr_enumeration, "ccd" );
	closuremethod.add_restriction( xsr_enumeration, "automatic" );
	closuremethod.add_restriction( xsr_enumeration, "kinematic" );
	xsd.add_top_level_element( closuremethod );

	attlist + XMLSchemaAttribute::attribute_w_default( "protocol", "closuremethod", "Method for loop closure", "ccd" )
		+ XMLSchemaAttribute::attribute_w_default( "perturb", xsct_rosetta_bool, "Perturb before closure?", "1" )
		// Unlike in LoopMoverFromCommandLine, refine has a default value independent of the
		// value of protocol. HUH?
		+ XMLSchemaAttribute::attribute_w_default( "refine", xsct_rosetta_bool, "Use a refinement protocol?", "1" );

	rosetta_scripts::attributes_for_parse_score_function( attlist, "refine_score" );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "perturb_score" );

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "auto_loops", xsct_rosetta_bool, "Automatically detect/assign loops", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool, "Do design!", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "hurry", xsct_rosetta_bool, "Fast mode", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "cycles", xsct_non_negative_integer, "cycles", "10" );


	attlist + XMLSchemaAttribute::attribute_w_default( "loop_file", xs_string, "File from which loops definitions might be read", "loops.loops" )
		+ XMLSchemaAttribute( "loop_start_pdb_num", xsct_refpose_enabled_residue_number, "Loop start, in PDB/seqpos/reference pose numbering" )
		+ XMLSchemaAttribute( "loop_start_resnum", xsct_refpose_enabled_residue_number, "Loop start, in PDB/seqpos/reference pose numbering" )
		+ XMLSchemaAttribute( "loop_end_pdb_num", xsct_refpose_enabled_residue_number, "Loop end, in PDB/seqpos/reference pose numbering" )

		+ XMLSchemaAttribute( "loop_end_pdbnum", xsct_refpose_enabled_residue_number, "Loop end, in PDB/seqpos/reference pose numbering" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string LoopRemodelCreator::keyname() const {
	return LoopRemodel::mover_name();
}

protocols::moves::MoverOP
LoopRemodelCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopRemodel );
}

void LoopRemodelCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopRemodel::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
