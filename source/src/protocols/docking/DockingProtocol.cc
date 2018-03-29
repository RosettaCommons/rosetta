// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   DockingProtocol.cc
///
/// @brief
/// @author Monica Berrondo
/// @author Sid Chaudhury
/// @author Jeff Gray
/// @author Modified by Daisuke Kuroda

// Unit Headers
#include <protocols/docking/DockingProtocol.hh>

// Package Headers
#include <protocols/docking/util.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/docking/DockingEnsemble.hh>
#include <protocols/docking/DockFilters.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingLowResEnsemble.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockingHighResFactory.hh>
#include <protocols/docking/DockingHighResLegacy.hh> // there is some design thing that only this flavor of DHR supports
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockingProtocolCreator.hh> // Support the scriptor

// Project Headers
#include <core/chemical/ChemicalManager.fwd.hh>


#include <core/pose/Pose.hh>

#include <protocols/simple_task_operations/InterfaceTaskOperation.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/RepackSidechainsMover.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <basic/datacache/DataMap.hh>

#include <core/io/raw_data/ScoreMap.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS

#include <core/types.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <algorithm>

#include <core/import_pose/import_pose.hh>
//#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.docking.DockingProtocol" );
using namespace core;

namespace protocols {
namespace docking {

// Constructors
DockingProtocol::DockingProtocol()
{
	user_defined_ = false;
	if_ensemble_ = false; // valgrind complains about uninitialized variable w/o this.
	init(utility::tools::make_vector1<core::SSize>(1), false, false, true, nullptr, nullptr);
}

DockingProtocol::DockingProtocol(
	Size const rb_jump,
	bool const low_res_protocol_only,
	bool const docking_local_refine,
	bool const autofoldtree,
	core::scoring::ScoreFunctionOP docking_score_low,
	core::scoring::ScoreFunctionOP docking_score_high
) : Mover()
{
	user_defined_ = true;
	init(utility::tools::make_vector1<core::SSize>(rb_jump), low_res_protocol_only, docking_local_refine, autofoldtree, docking_score_low, docking_score_high);
}

DockingProtocol::DockingProtocol(
	DockJumps const movable_jumps,
	bool const low_res_protocol_only,
	bool const docking_local_refine,
	bool const autofoldtree,
	core::scoring::ScoreFunctionOP docking_score_low,
	core::scoring::ScoreFunctionOP docking_score_high
) : Mover()
{
	user_defined_ = true;
	init(movable_jumps, low_res_protocol_only, docking_local_refine, autofoldtree, docking_score_low, docking_score_high);
}
// End of constructors


void DockingProtocol::init(
	DockJumps const movable_jumps,
	bool const low_res_protocol_only,
	bool const docking_local_refine,
	bool const autofoldtree,
	core::scoring::ScoreFunctionOP docking_score_low,
	core::scoring::ScoreFunctionOP docking_score_high
)
{

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	Mover::type( "DockingProtocol" );
	movable_jumps_ = movable_jumps;

	// setup all the booleans with default values
	// they will get overwritten by the options and/or passed values
	set_default();
	init_from_options();
	if ( user_defined_ ) {
		low_res_protocol_only_ = low_res_protocol_only;
		docking_local_refine_ = docking_local_refine;
		autofoldtree_ = autofoldtree;
	}

	// correctly set up the score functions from either passed in values or defaults
	if ( docking_score_low == nullptr ) {
		//docking_scorefxn_low_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
		docking_scorefxn_low_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen",
			basic::options::option[ basic::options::OptionKeys::score::patch ]()
		);

	} else {
		docking_scorefxn_low_ = docking_score_low;
	}

	if ( docking_score_high == nullptr ) {
		docking_scorefxn_high_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" ) ;
		docking_scorefxn_pack_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) ;
		docking_scorefxn_output_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking" );
	} else {
		docking_scorefxn_high_ = docking_score_high;
		docking_scorefxn_pack_ = docking_score_high;
		docking_scorefxn_output_ = docking_score_high;
	}

	if ( option[ score::weights ].user() ) {
		docking_scorefxn_high_ = core::scoring::get_score_function();
		docking_scorefxn_output_ = docking_scorefxn_high_;
	}

	if ( option[ score::pack_weights ].user() ) {
		docking_scorefxn_pack_ = core::scoring::ScoreFunctionFactory::create_score_function( option[ score::pack_weights ]() );
	}

	// set up objects based on the boolean values defined above
	setup_objects();
}

void
DockingProtocol::set_default()
{
	//using namespace basic::options;
	using namespace core::scoring;

	reporting_ = true;
	autofoldtree_ = true;

	low_res_protocol_only_ = false;
	docking_local_refine_ = false;
	dock_min_ = false;
	rt_min_ = false;
	sc_min_ = false;
	partners_ = "_";
	cst_weight_ = 0.0;
	cst_fa_weight_ = 0.0;
	use_csts_ = false;
	score_cutoff_ = 1000000.0;
	no_filters_ = false;
	use_legacy_protocol_ = false;
	ignore_default_docking_task_ = false;
	design_ = false;
	if_ensemble_ = false;

	lowres_inner_cycles_ = 50;
	lowres_outer_cycles_ = 10;

	// initialize the ensemble movers and file strings
	ensemble1_ = nullptr;
	ensemble2_ = nullptr;
	ensemble1_filename_ = "";
	ensemble2_filename_ = "";
	recover_sidechains_filename_ = "";

}

void DockingProtocol::setup_objects()
{
	// initialize to null
	lowres_filter_ = nullptr;
	highres_filter_ = nullptr;

	// initialize constraint set mover
	docking_constraint_ = nullptr;

	// stores the sequence of the previous pose, so that the DockingProtocol can re setup the fold tree
	previous_sequence_ = "";

	fold_tree_ = core::kinematics::FoldTree(); // apl NOTE: was NULL, but fold_tree_ is not a pointer.
	init_task_factory_ = nullptr;
	perturber_ = nullptr;
	docking_lowres_mover_ = nullptr;
	docking_highres_mover_ = nullptr;

	// Residue movers
	to_centroid_ = protocols::simple_moves::SwitchResidueTypeSetMoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	//generate to_all_atom mover:  to_all_atom_ =
	using protocols::moves::MoverOP;
	protocols::moves::SequenceMoverOP to_all_atom_and_repack( new protocols::moves::SequenceMover );
	to_all_atom_and_repack->add_mover( MoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) ) );
	to_all_atom_and_repack->add_mover( MoverOP( new protocols::minimization_packing::RepackSidechainsMover( docking_scorefxn_pack_ ) ) );
	to_all_atom_=to_all_atom_and_repack;

	sync_objects_with_flags();
}

void DockingProtocol::sync_objects_with_flags()
{
	if ( !docking_local_refine_ ) {
		if ( ensemble2_filename_ != "" ) {
			if_ensemble_ = true;
		}

		if ( !perturber_ ) {
			perturber_ = protocols::docking::DockingInitialPerturbationOP( new DockingInitialPerturbation( movable_jumps_, true /*slide into contact*/ ) );
		}

		if ( !docking_lowres_mover_ ) {
			// Modified by DK
			if ( if_ensemble_ ) {
				docking_lowres_mover_ = protocols::docking::DockingLowResOP( new DockingLowResEnsemble( docking_scorefxn_low_, movable_jumps_ ) );
			} else {
				docking_lowres_mover_ = protocols::docking::DockingLowResOP( new DockingLowRes( docking_scorefxn_low_, movable_jumps_ ) );
			}
		}
		if ( !no_filters_ && !lowres_filter_ ) {
			lowres_filter_ = protocols::docking::DockingLowResFilterOP( new protocols::docking::DockingLowResFilter() );
		}
	} else {
		perturber_ = nullptr;
		docking_lowres_mover_ = nullptr;
	}

	if ( !low_res_protocol_only_ ) {
		DHR_Type t = DHR_Type::DockMCMProtocol; // default
		if ( dock_min_ ) {
			t = DHR_Type::DockMinMover;
			if ( docking_highres_mover_ && docking_highres_mover_->get_name() != "DockMinMover" ) {
				docking_highres_mover_ = nullptr;
			}
		} else if ( use_legacy_protocol_ ) {
			t = DHR_Type::DockingHighResLegacy;
			if ( docking_highres_mover_ && docking_highres_mover_->get_name() != "DockingHighResLegacy" ) {
				docking_highres_mover_ = nullptr;
			}
		} else {
			if ( docking_highres_mover_ && docking_highres_mover_->get_name() != "DockMCMProtocol" ) {
				docking_highres_mover_ = nullptr;
			}
		}
		// TODO: Add support for SnugDock here.
		// create and configure the appropriate flavor of DockingHighRes
		if ( ! docking_highres_mover_ ) {
			docking_highres_mover_ = DockingHighResFactory::get_instance()->create_docking_high_res_mover( t );
			docking_highres_mover_->set_movable_jumps( movable_jumps_ );
			docking_highres_mover_->set_scorefxn( docking_scorefxn_high_ );
			docking_highres_mover_->set_scorefxn_pack( docking_scorefxn_pack_ );
			docking_highres_mover_->set_rt_min( rt_min_ );
			docking_highres_mover_->set_sc_min( sc_min_ );
			docking_highres_mover_->set_partners( partners_ );
		}
		if ( !no_filters_ && !highres_filter_ ) {
			highres_filter_ = protocols::docking::DockingHighResFilterOP( new protocols::docking::DockingHighResFilter() );
		}
	} else {
		docking_highres_mover_ = nullptr;
	}
	if ( no_filters_ ) {
		highres_filter_ = nullptr;
		lowres_filter_ = nullptr;
	}

	// set up constraint set mover
	if ( !use_csts_ ) {
		docking_constraint_ = nullptr;
	} else {
		if ( !docking_constraint_ ) {
			docking_constraint_ = protocols::constraint_movers::ConstraintSetMoverOP( new protocols::constraint_movers::ConstraintSetMover() );
		}
	}

	flags_and_objects_are_in_sync_ = true;
	first_apply_with_current_setup_ = true;
}

void
DockingProtocol::init_from_options()
{
	using namespace basic::options;

	// check for movable jumps option
	if ( option[ OptionKeys::docking::multibody ].user() ) {
		set_movable_jumps( option[ OptionKeys::docking::multibody ]() );
	}

	// These options default to false
	if ( option[ OptionKeys::docking::low_res_protocol_only ].user() ) {
		set_low_res_protocol_only(option[ OptionKeys::docking::low_res_protocol_only ]());
	}

	if ( option[ OptionKeys::docking::docking_local_refine ].user() ) {
		set_docking_local_refine(option[ OptionKeys::docking::docking_local_refine ]());
	}

	if ( option[ OptionKeys::docking::dock_min ].user() ) {
		set_dock_min(option[ OptionKeys::docking::dock_min ]());
	}

	// Sets the member directly so sync_objects_with_flags won't be triggered prematurely.
	// A public setter exists for this member.
	if ( option[ OptionKeys::docking::dock_rtmin ].user() ) {
		rt_min_ = option[ OptionKeys::docking::dock_rtmin ]();
	}

	// Sets the member directly so sync_objects_with_flags won't be triggered prematurely.
	// A public setter exists for this member.
	if ( option[ OptionKeys::docking::sc_min ].user() ) {
		sc_min_ = option[ OptionKeys::docking::sc_min ]();
	}

	// This defaults to "_"
	if ( option[ OptionKeys::docking::partners ].user() ) {
		set_partners(option[ OptionKeys::docking::partners ]());
	}

	// Defaults to 0
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		cst_weight_ = option[ OptionKeys::constraints::cst_weight ]();
	}

	// Defaults to 0
	if ( option[ OptionKeys::constraints::cst_fa_weight ].user() ) {
		cst_fa_weight_ = option[ OptionKeys::constraints::cst_fa_weight ]();
	}

	// Defaults to false
	set_use_constraints( option[ OptionKeys::constraints::cst_file ].user() ||
		option[ OptionKeys::constraints::cst_fa_file ].user() );

	//set native pose if asked for
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_file( *native_pose, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		set_native_pose( native_pose );
	} else {
		set_native_pose(nullptr);
	}

	if ( option[ OptionKeys::cluster::output_score_filter ].user() ) {
		score_cutoff_ = option[ OptionKeys::cluster::output_score_filter ]();
	}

	if ( option[ OptionKeys::docking::no_filters ].user() ) {
		set_no_filters(option[ OptionKeys::docking::no_filters ]);
	}

	if ( option[ OptionKeys::docking::use_legacy_protocol ].user() ) {
		set_use_legacy_protocol( option[ OptionKeys::docking::use_legacy_protocol ] );
	}

	if ( option[ OptionKeys::docking::ignore_default_docking_task].user() ) {
		set_ignore_default_docking_task( option[ OptionKeys::docking::ignore_default_docking_task ]() );
	}
	// docking low res options
	if ( option[ OptionKeys::docking::docking_centroid_inner_cycles ].user() ) {
		set_inner_cycles(option[ OptionKeys::docking::docking_centroid_inner_cycles ]() );
	}

	if ( option[ OptionKeys::docking::docking_centroid_outer_cycles ].user() ) {
		set_outer_cycles(option[ OptionKeys::docking::docking_centroid_outer_cycles ]());
	}

	// above cycles must be called before ensemble since the cycle number is overwritten if in ensemble mode
	if ( option[ OptionKeys::docking::ensemble1 ].user() ) {
		set_ensemble1(option[ OptionKeys::docking::ensemble1 ]());
	}

	if ( option[ OptionKeys::docking::ensemble2 ].user() ) {
		set_ensemble2(option[ OptionKeys::docking::ensemble2 ]());
	}

	if ( option[ OptionKeys::docking::recover_sidechains ].user() ) {
		set_recover_sidechains_filename(option[ OptionKeys::docking::recover_sidechains ]());
	}
}

void
DockingProtocol::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::multibody );
	option.add_relevant( OptionKeys::docking::low_res_protocol_only );
	option.add_relevant( OptionKeys::docking::docking_local_refine );
	option.add_relevant( OptionKeys::docking::dock_min );
	option.add_relevant( OptionKeys::docking::dock_rtmin );
	option.add_relevant( OptionKeys::docking::sc_min );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::docking::no_filters );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::constraints::cst_fa_weight );

	option.add_relevant( OptionKeys::run::score_only );
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::docking::use_legacy_protocol );
	option.add_relevant( OptionKeys::docking::ignore_default_docking_task );
	option.add_relevant( OptionKeys::docking::recover_sidechains );
	// low res protocol options
	option.add_relevant( OptionKeys::docking::docking_centroid_inner_cycles );
	option.add_relevant( OptionKeys::docking::docking_centroid_outer_cycles );
	option.add_relevant( OptionKeys::docking::ensemble1 );
	option.add_relevant( OptionKeys::docking::ensemble2 );

	// This takes care of the cst_file flags.
	constraint_movers::ConstraintSetMover::register_options();
}

void
DockingProtocol::finalize_setup( pose::Pose & pose ) //setup objects requiring pose during apply
{
	if ( autofoldtree_ ) {
		TR << "Setting docking foldtree" << std::endl;
		TR << "old fold tree: " << pose.fold_tree() << std::endl;
		setup_foldtree( pose, partners_, movable_jumps_ );
		TR << "new fold tree: " << pose.fold_tree() << std::endl;
	}

	fold_tree_ = pose.fold_tree();

	// check for native and input pose
	// input pose is used further down for the recover sidehchains mover
	if ( !get_input_pose() ) {
		core::pose::PoseOP input_pose( new core::pose::Pose( pose ) );
		set_input_pose( input_pose );
	}

	// if an ensemble file is defined for either partner, both partners must be defined
	if ( ensemble1_filename_ != "" || ensemble2_filename_ != "" ) {
		if ( ensemble1_filename_ == "" || ensemble2_filename_ == "" ) utility_exit_with_message( "Must define ensemble file for both partners");
		runtime_assert( movable_jumps_.size() == 1 ); // ensemble mode is only allowed for traditional docking
		core::Size const rb_jump = movable_jumps_[1];
		Size start_res(1), end_res(1), cutpoint(pose.fold_tree().cutpoint_by_jump( rb_jump ));

		//lowres_inner_cycles_ = 25; // Should be 50 (default value for traditional docking); modified by DK

		TR << "Ensemble 1: " << ensemble1_filename_ << std::endl;
		start_res = 1;
		end_res = cutpoint;

		ensemble1_ = protocols::docking::DockingEnsembleOP( new DockingEnsemble( start_res, end_res, rb_jump, ensemble1_filename_, "dock_ens_conf1", docking_scorefxn_low_, docking_scorefxn_high_ ) );

		TR << "Ensemble 2: " << ensemble2_filename_ << std::endl;
		start_res = cutpoint + 1;
		end_res = pose.size();

		ensemble2_ = protocols::docking::DockingEnsembleOP( new DockingEnsemble( start_res, end_res, rb_jump, ensemble2_filename_, "dock_ens_conf2", docking_scorefxn_low_, docking_scorefxn_high_ ) );

		// recover sidechains mover is not needed with ensemble docking since the sidechains are recovered from the partners in the ensemble file
		recover_sidechains_ = nullptr;

		// add ensemble conformer energies to low res scorefunction
		core::scoring::ScoreFunctionOP docking_scorefxn_ens, docking_highres_ens;
		docking_scorefxn_ens = docking_scorefxn_low_->clone();
		docking_scorefxn_ens->set_weight( core::scoring::dock_ens_conf, 1.0 );
		set_lowres_scorefxn( docking_scorefxn_ens );
		// pass the scorefunction to the low res mover
		docking_lowres_mover_->set_scorefxn( docking_scorefxn_low_ );

		//docking_highres_ens = new core::scoring::ScoreFunction( *docking_scorefxn_high_ );
		//docking_highres_ens->set_weight( core::scoring::dock_ens_conf, 1.0 );
		//set_highres_scorefxn( docking_highres_ens, docking_scorefxn_pack_ ); // sets csts for mc and minimization, but not packing
		// pass the score function to the high res mover
		//docking_highres_mover_->set_scorefxn( docking_scorefxn_high_ );
	} else { //if ensemble docking
		if ( recover_sidechains_filename_ != "" ) {
			if ( !recover_sidechains_ ) {
				core::pose::Pose a_pose;
				core::import_pose::pose_from_file( a_pose, recover_sidechains_filename_ , core::import_pose::PDB_file);
				recover_sidechains_ = protocols::simple_moves::ReturnSidechainMoverOP( new protocols::simple_moves::ReturnSidechainMover( a_pose ) );
			} //first initialization ?
		} else if ( get_input_pose() && get_input_pose()->is_fullatom() ) {
			recover_sidechains_ = protocols::simple_moves::ReturnSidechainMoverOP( new protocols::simple_moves::ReturnSidechainMover( *get_input_pose() ) );
		} else {
			// recover sidechains mover is not needed with ensemble docking since the sidechains are recovered from the partners in the ensemble fi
			recover_sidechains_ = nullptr;
		}
	} //if ensemble docking

	if ( docking_lowres_mover_ ) {
		// pass the ensemble movers to the lowres protocol
		if ( if_ensemble_ && docking_lowres_mover_->get_name() == "DockingLowResEnsemble" ) {
			auto* ensemble_mover = dynamic_cast< DockingLowResEnsemble* >(docking_lowres_mover_.get());

			ensemble_mover->set_ensemble1( ensemble1_ );
			ensemble_mover->set_ensemble2( ensemble2_ );
		}

		docking_lowres_mover_->set_inner_cycles( lowres_inner_cycles_ );
		docking_lowres_mover_->set_outer_cycles( lowres_outer_cycles_ );
	}

	// set relevant information to legacy high res mover
	if ( docking_highres_mover_ ) {
		if ( docking_highres_mover_->get_name() == "DockingHighResLegacy" && design_ ) {
			// TODO: Why are we still relying on DockingHighResLegacy for ANYTHING let alone for design?
			// We gotta find a way to move this functionality into the DockingHighRes base class (which would require some testing)

			// FIXME: This one line is the reason we have to #include the DockingHighResLegacy header >:O
			DockingHighResLegacyOP legacy_mover = utility::pointer::dynamic_pointer_cast< protocols::docking::DockingHighResLegacy > ( docking_highres_mover_ );
			legacy_mover->design( design_ );
		}
		// passes the task factory down the chain and allows setting of the default docking task
		if ( init_task_factory_  ) {
			docking_highres_mover_->set_task_factory( init_task_factory_ );
		}
		docking_highres_mover_->set_ignore_default_task( ignore_default_docking_task_ );
	}

	//initialize docking filters
	if ( highres_filter_ ) highres_filter_->set_score_cutoff( score_cutoff_ );

	if ( docking_constraint_ ) {
		TR << "setting up the constraint set mover" << std::endl;
		docking_constraint_->apply( pose );
		if ( cst_weight_ == 0.0 ) {
			cst_weight_ = 1.0;
		}

		if ( cst_fa_weight_ == 0.0 ) {
			cst_fa_weight_ = 1.0;
		}
	}
	// finish setting up constraints
	setup_constraints( pose );
}

//destructor
DockingProtocol::~DockingProtocol() = default;

/// @brief clone operator, calls the copy constructor
protocols::moves::MoverOP
DockingProtocol::clone() const {
	//return( new DockingProtocol( movable_jumps_, low_res_protocol_only_, docking_local_refine_, autofoldtree_, docking_scorefxn_low_, docking_scorefxn_high_ ) ); This is bad do not clone this way.
	return protocols::moves::MoverOP( new DockingProtocol(*this) );
}

/// @brief fresh_instance returns a default-constructed object for the JD2
protocols::moves::MoverOP
DockingProtocol::fresh_instance() const {
	return protocols::moves::MoverOP( new DockingProtocol() );
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool DockingProtocol::reinitialize_for_new_input() const {
	return true;
}

/// @brief copy ctor
DockingProtocol::DockingProtocol( DockingProtocol const & rhs ) :
	Mover(rhs),
	if_ensemble_( false ) // valgrind complains about uninitialized variable w/o this.
{
	initForEqualOperatorAndCopyConstructor(*this, rhs);
}

/// @brief assignment operator
DockingProtocol & DockingProtocol::operator=( DockingProtocol const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	Mover::operator=(rhs);
	initForEqualOperatorAndCopyConstructor(*this, rhs);
	return *this;
}

void DockingProtocol::initForEqualOperatorAndCopyConstructor(DockingProtocol & lhs, DockingProtocol const & rhs){
	//going through all of member data and assigning it
	lhs.user_defined_ = rhs.user_defined_;
	lhs.low_res_protocol_only_ = rhs.low_res_protocol_only_;
	lhs.reporting_ = rhs.reporting_;
	lhs.autofoldtree_ = rhs.autofoldtree_;
	lhs.flags_and_objects_are_in_sync_ = rhs.flags_and_objects_are_in_sync_;
	lhs.first_apply_with_current_setup_ = rhs.first_apply_with_current_setup_;
	lhs.sc_min_ = rhs.sc_min_;
	lhs.rt_min_ = rhs.rt_min_;
	lhs.dock_min_ = rhs.dock_min_;
	lhs.no_filters_ = rhs.no_filters_;
	lhs.use_legacy_protocol_ =  rhs.use_legacy_protocol_;
	lhs.docking_local_refine_ = rhs.docking_local_refine_;
	lhs.use_csts_ = rhs.use_csts_;
	lhs.cst_weight_ = rhs.cst_weight_;
	lhs.cst_fa_weight_ = rhs.cst_fa_weight_;
	lhs.score_cutoff_ = rhs.score_cutoff_;
	lhs.fold_tree_ = rhs.fold_tree_;
	lhs.partners_ = rhs.partners_;
	lhs.previous_sequence_ = rhs.previous_sequence_;
	lhs.lowres_inner_cycles_ = rhs.lowres_inner_cycles_;
	lhs.lowres_outer_cycles_ = rhs.lowres_outer_cycles_;
	lhs.movable_jumps_ = rhs.movable_jumps_;
	if ( rhs.docking_scorefxn_low_ ) lhs.docking_scorefxn_low_ = rhs.docking_scorefxn_low_->clone();
	if ( rhs.docking_scorefxn_high_ ) lhs.docking_scorefxn_high_ = rhs.docking_scorefxn_high_->clone();
	if ( rhs.docking_scorefxn_pack_ ) lhs.docking_scorefxn_pack_ = rhs.docking_scorefxn_pack_->clone();
	if ( rhs.docking_scorefxn_output_ ) lhs.docking_scorefxn_output_ = rhs.docking_scorefxn_output_->clone();
	if ( rhs.mc_ ) { //not used currently but might be needed later
		lhs.mc_ = protocols::moves::MonteCarloOP( new moves::MonteCarlo( *(rhs.mc_) ) );
	}
	if ( rhs.lowres_filter_ ) lhs.lowres_filter_ = utility::pointer::static_pointer_cast< protocols::docking::DockingLowResFilter > ( rhs.lowres_filter_->clone() );
	if ( rhs.highres_filter_ ) lhs.highres_filter_ = utility::pointer::static_pointer_cast< protocols::docking::DockingHighResFilter > ( rhs.highres_filter_->clone() );
	if ( rhs.docking_lowres_mover_ ) {
		lhs.docking_lowres_mover_ = utility::pointer::static_pointer_cast< protocols::docking::DockingLowRes > ( rhs.docking_lowres_mover_->clone() );
	}
	if ( rhs.docking_highres_mover_ ) {
		lhs.docking_highres_mover_ = utility::pointer::static_pointer_cast< protocols::docking::DockingHighRes > ( rhs.docking_highres_mover_->clone() );
	}
	if ( rhs.to_centroid_ ) lhs.to_centroid_ = utility::pointer::static_pointer_cast< protocols::simple_moves::SwitchResidueTypeSetMover > ( rhs.to_centroid_->clone() );
	if ( rhs.to_all_atom_ ) lhs.to_all_atom_ = rhs.to_all_atom_->clone();
	if ( rhs.ensemble1_ ) {
		lhs.ensemble1_ = protocols::docking::DockingEnsembleOP( new protocols::docking::DockingEnsemble( *(rhs.ensemble1_) ) );
		lhs.ensemble1_filename_ = rhs.ensemble1_filename_ ;
	}
	if ( rhs.ensemble2_ ) {
		lhs.ensemble2_ = protocols::docking::DockingEnsembleOP( new protocols::docking::DockingEnsemble( *(rhs.ensemble2_) ) );
		lhs.ensemble2_filename_ = rhs.ensemble2_filename_;
	}
	if ( rhs.docking_constraint_ ) {
		lhs.docking_constraint_ = utility::pointer::static_pointer_cast< protocols::constraint_movers::ConstraintSetMover > ( rhs.docking_constraint_->clone() );
	}
	if ( rhs.recover_sidechains_ ) {
		lhs.recover_sidechains_ = utility::pointer::static_pointer_cast< protocols::simple_moves::ReturnSidechainMover > ( rhs.recover_sidechains_->clone() );
	}
	if ( rhs.init_task_factory_ ) {
		lhs.init_task_factory_ = core::pack::task::TaskFactoryOP( new  core::pack::task::TaskFactory( *(rhs.init_task_factory_) ) );
	}
	lhs.design_ = rhs.design_;
	lhs.ignore_default_docking_task_ = rhs.ignore_default_docking_task_;
}

void DockingProtocol::set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_low )
{
	docking_scorefxn_low_ = docking_scorefxn_low;

}

void DockingProtocol::set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_high )
{
	docking_scorefxn_high_ = docking_scorefxn_high;
	docking_scorefxn_pack_ = docking_scorefxn_high;
	docking_scorefxn_output_ = docking_scorefxn_high;
}

void DockingProtocol::set_highres_scorefxn( // delete
	core::scoring::ScoreFunctionOP docking_scorefxn_high,
	core::scoring::ScoreFunctionOP docking_scorefxn_pack )
{
	docking_scorefxn_high_ = docking_scorefxn_high;
	docking_scorefxn_pack_ = docking_scorefxn_pack;
	docking_scorefxn_output_ = docking_scorefxn_high;
}

void DockingProtocol::set_highres_scorefxn(
	core::scoring::ScoreFunctionOP docking_scorefxn_high,
	core::scoring::ScoreFunctionOP docking_scorefxn_pack,
	core::scoring::ScoreFunctionOP docking_scorefxn_output )
{
	docking_scorefxn_high_ = docking_scorefxn_high;
	docking_scorefxn_pack_ = docking_scorefxn_pack;
	docking_scorefxn_output_ = docking_scorefxn_output;
}

void DockingProtocol::set_sc_min( bool sc_min )
{
	check_high_res_protocol();
	sc_min_ = sc_min;
	docking_highres_mover_->set_sc_min( sc_min );
}

void DockingProtocol::set_rt_min( bool rt_min )
{
	check_high_res_protocol();
	rt_min_ = rt_min;
	docking_highres_mover_->set_rt_min( rt_min );
}

void DockingProtocol::set_dock_min( bool const dock_min )
{
	//debug_assert ( !low_res_protocol_only_ ); // SJF this assert doesn't make sense b/c set_dock_min is called even if low_res_protocol_only_ is set through parse_my_tag
	dock_min_ = dock_min;
}

void DockingProtocol::set_no_filters( bool no_filters )
{
	if ( no_filters != no_filters_ ) {
		no_filters_ = no_filters;
		flags_and_objects_are_in_sync_ = false;
	}
}

void DockingProtocol::set_low_res_protocol_only( bool const low_res_protocol_only )
{
	if ( low_res_protocol_only != low_res_protocol_only_ ) {
		low_res_protocol_only_ = low_res_protocol_only;
		flags_and_objects_are_in_sync_ = false;
	}
}

void DockingProtocol::set_docking_local_refine( bool const docking_local_refine )
{
	if ( docking_local_refine != docking_local_refine_ ) {
		docking_local_refine_ = docking_local_refine;
		flags_and_objects_are_in_sync_ = false;
	}
}

void DockingProtocol::set_use_legacy_protocol( bool const use_legacy_protocol )
{
	if ( use_legacy_protocol != use_legacy_protocol_ ) {
		use_legacy_protocol_ = use_legacy_protocol;
		flags_and_objects_are_in_sync_ = false;
	}

}

/// @details Set the constraint weight for both centroid and full atom constraints to the same value
void DockingProtocol::set_cst_weight( core::Real const cst_weight )
{
	// Set the centroid constraint weight
	if ( cst_weight != cst_weight_ ) {
		cst_weight_ = cst_weight;
		runtime_assert( cst_weight_ != 0.0 );
		flags_and_objects_are_in_sync_ = false;
	}

	// Set the full atom constraint weight
	if ( cst_weight != cst_fa_weight_ ) {
		cst_fa_weight_ = cst_weight;
		runtime_assert( cst_fa_weight_ != 0.0 );
		flags_and_objects_are_in_sync_ = false;
	}
}

void DockingProtocol::set_use_constraints( bool const use_csts )
{
	if ( use_csts != use_csts_ ) {
		use_csts_ = use_csts;
		flags_and_objects_are_in_sync_ = false;
	}
}

void DockingProtocol::set_interface_definition_task_operation( protocols::simple_task_operations::InterfaceTaskOperationOP interface_definition )
{
	check_high_res_protocol();
	docking_highres_mover_->set_interface_definition_task_operation( interface_definition );
}

void DockingProtocol::set_additional_task_operarations( utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations )
{
	check_high_res_protocol();
	docking_highres_mover_->set_additional_task_operarations( additional_task_operations );
}

void DockingProtocol::add_additional_task_operaration( core::pack::task::operation::TaskOperationOP task_operation )
{
	check_high_res_protocol();
	docking_highres_mover_->add_additional_task_operaration( task_operation );
}

utility::vector1< core::pack::task::operation::TaskOperationOP > DockingProtocol::get_additional_task_operarations()
{
	check_high_res_protocol();
	return docking_highres_mover_->get_additional_task_operarations();
}


void DockingProtocol::check_high_res_protocol()
{
	debug_assert ( !low_res_protocol_only_ );
	if ( !docking_highres_mover_ ) {
		sync_objects_with_flags();
	}
}

/// @brief setup the constrainta for each application of the docking protocol
void DockingProtocol::setup_constraints( pose::Pose & pose )
{
	bool filter_use_csts = false;

	// Only checking cst_weight becuase this will only affect the lowres_filter
	if ( cst_weight_ != 0.0 ) {
		if ( ! pose.constraint_set() ) {
			TR << "No ConstraintSet has been created, skipping cst_weight. "
				<< "Use cst_file instead and set up desired ConstraintSet" << std::endl;
		} else {
			filter_use_csts = true;
		}
	}

	if ( lowres_filter_ ) { lowres_filter_->set_use_constraints( filter_use_csts ); }

	add_constraints_to_scorefunction();
}

/// @brief add distance constraint to the docking score functions
void DockingProtocol::add_constraints_to_scorefunction()
{
	// Add constraints to the low resolution docking mover if applicable
	if ( docking_lowres_mover_ && cst_weight_ != 0.0 ) {
		core::scoring::ScoreFunctionOP docking_scorefxn_cst = docking_scorefxn_low_->clone();
		docking_scorefxn_cst->set_weight( core::scoring::atom_pair_constraint, cst_weight_ );
		set_lowres_scorefxn( docking_scorefxn_cst );
		// pass the scorefunction to the low res mover
		docking_lowres_mover_->set_scorefxn( docking_scorefxn_low_ );
	}

	// Add constraints to the high resolution docking mover if applicable
	if ( docking_highres_mover_ && cst_fa_weight_ != 0.0 ) {
		core::scoring::ScoreFunctionOP docking_highres_cst = docking_scorefxn_output_->clone(); // needs to use the non-min score function to match assemble_domains test with constraints
		docking_highres_cst->set_weight( core::scoring::atom_pair_constraint, cst_fa_weight_ );
		set_highres_scorefxn( docking_highres_cst, docking_scorefxn_pack_ ); // sets csts for mc and minimization, but not packing
		// pass the score function to the high res mover
		docking_highres_mover_->set_scorefxn( docking_scorefxn_high_ );
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief main docking protocol
///
/// @details
/// Main function for docking. Includes the following steps:
///      0) prepack mode: prepare a starting structure for later runs
///   OR:
///      1) perturbation of decoy (see docking_perturb_decoy): changes
///         orientation of docking partners
///      2) low-resolution search:
///         refine structure with rigid-body, centroid-mode MC cycles
///      3) high-resolution search:
///         further refinement of structure in fullatom mode
void
DockingProtocol::apply( pose::Pose & pose )
{
	bool passed_lowres_filter = true;
	bool passed_highres_filter = true;
	std::map < std::string, core::Real > lowres_scores;

	if ( !get_native_pose() ) {
		TR << "Danger Will Robinson! Native is an impostor!" << std::endl;
		core::pose::PoseOP native_pose( new core::pose::Pose(pose) );
		set_native_pose( native_pose );
	}

	if ( !flags_and_objects_are_in_sync_ ) {
		sync_objects_with_flags();
	}

	if ( previous_sequence_.compare( pose.sequence() ) != 0 ) {
		first_apply_with_current_setup_ = true;
		previous_sequence_ = pose.sequence();
	}

	if ( first_apply_with_current_setup_ ) {
		finalize_setup( pose );
		first_apply_with_current_setup_ = false;
	}

	pose.fold_tree( fold_tree_ );
	show(TR);

	basic::prof_reset();

	// Low resolution docking
	if ( docking_lowres_mover_ ) {
		// convert to centroid mode
		to_centroid_->apply( pose );
		if ( docking_constraint_ ) {
			TR << "setting up the constraint set mover" << std::endl;
			docking_constraint_->apply( pose );
			if ( cst_weight_ == 0.0 ) {
				cst_weight_ = 1.0;
			}
			// finish setting up constraints
			setup_constraints( pose );
		}

		TR << pose.fold_tree() << std::endl;
		// make starting perturbations based on command-line flags over each movable jump
		if ( perturber_ ) {
			perturber_->apply( pose );
			TR << "finished initial perturbation" << std::endl;
		}

		// add scores to jd2 output
		if ( reporting_ && get_native_pose() ) protocols::jd2::add_string_real_pair_to_current_job("st_rmsd", calc_Lrmsd( pose, *get_native_pose(), movable_jumps_ ));

		docking_lowres_mover_->apply( pose );

		// Perform additional low resolution steps if they have been specified
		if ( additional_low_resolution_steps_ ) {
			additional_low_resolution_steps_->apply( pose );
		}

		// add scores to jd2 output
		if ( reporting_ ) {
			if ( ensemble1_ ) protocols::jd2::add_string_real_pair_to_current_job("conf_num1", ensemble1_->get_current_confnum() );
			if ( ensemble2_ ) protocols::jd2::add_string_real_pair_to_current_job("conf_num2", ensemble2_->get_current_confnum() );
			if ( get_native_pose() ) {
				// calculate and store the rms no matter which mode was used
				protocols::jd2::add_string_real_pair_to_current_job("rms", calc_Lrmsd( pose, *get_native_pose(), movable_jumps_ ));
				protocols::jd2::add_string_real_pair_to_current_job("cen_rms", calc_Lrmsd( pose, *get_native_pose(), movable_jumps_ ));
			}
			// store the low res scores for output
			core::io::raw_data::ScoreMap::score_map_from_scored_pose( lowres_scores, pose );
		}
		if ( lowres_filter_ ) {
			passed_lowres_filter = lowres_filter_->apply( pose );
			if ( !passed_lowres_filter ) lowres_filter_->report( TR, pose );
		}
	}

	// High resolution docking
	if ( passed_lowres_filter && docking_highres_mover_ ) {
		if ( docking_constraint_ ) {
			TR << "setting up the constraint set mover" << std::endl;
			docking_constraint_->apply( pose );
			if ( cst_fa_weight_ == 0.0 ) {
				cst_fa_weight_ = 1.0;
			}
			// finish setting up constraints
			setup_constraints( pose );
		}

		if ( ! pose.is_fullatom() ) {
			// Convert pose to high resolution and recover sidechains
			to_all_atom_->apply( pose );

			( * docking_highres_mover_->scorefxn() )( pose );
			jd2::write_score_tracer( pose, "Docking_to_all_atom" );
			if ( recover_sidechains_ ) {
				recover_sidechains_->apply( pose );
				( * docking_highres_mover_->scorefxn() )( pose );
				jd2::write_score_tracer( pose, "Docking_recovered_sidechains" );
			} else {
				if ( ensemble1_ ) {
					ensemble1_->recover_conformer_sidechains( pose );
				}
				if ( ensemble2_ ) {
					ensemble2_->recover_conformer_sidechains( pose );
				}
			}
		}
		docking_highres_mover_->apply( pose );
		( * docking_scorefxn_output_ )( pose );
		if ( highres_filter_ ) passed_highres_filter = highres_filter_->apply( pose );

		if ( reporting_ ) {
			// calculate and store the rms no matter which mode was used
			// this will overwrite the "rms" column from low res if it was calculated
			if ( get_native_pose() ) protocols::jd2::add_string_real_pair_to_current_job("rms", calc_Lrmsd( pose, *get_native_pose(), movable_jumps_ ) );

			// output low res scores if low res was run
			if ( !lowres_scores.empty() ) { //size() > 0 ) {
				for ( std::map< std::string, core::Real >::const_iterator pair=lowres_scores.begin(); pair!=lowres_scores.end(); ++pair ) {
					if ( pair->first == "dock_ens_conf" ) protocols::jd2::add_string_real_pair_to_current_job( "cen_dock_ens_conf", pair->second );
					else if ( pair->first != "total_score" ) protocols::jd2::add_string_real_pair_to_current_job( pair->first, pair->second );
				}
			}

			// TODO: metrics doesn't need to always take scorefunctions (for example, not necessary for irms or fnat)
			protocols::jd2::add_string_real_pair_to_current_job("I_sc", calc_interaction_energy( pose, docking_scorefxn_output_, movable_jumps_ ) );
			if ( get_native_pose() ) {
				protocols::jd2::add_string_real_pair_to_current_job("Irms", calc_Irmsd(pose, *get_native_pose(), docking_scorefxn_output_, movable_jumps_ ));
				protocols::jd2::add_string_real_pair_to_current_job("Fnat", calc_Fnat( pose, *get_native_pose(), docking_scorefxn_output_, movable_jumps_ ));
			}

			// pose.energies().show_total_headers( std::cout );
			// pose.energies().show_totals( std::cout );
			// TR << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
			// TR << pose.energies().total_energies()[ core::scoring::dock_ens_conf ] << std::endl;

			if ( ensemble1_ && ensemble2_ ) {
				protocols::jd2::add_string_real_pair_to_current_job("conf_score",
					pose.energies().total_energies()[ core::scoring::total_score ]
					- ensemble1_->highres_reference_energy()
					- ensemble2_->highres_reference_energy() );
			}
		}
	}

	if ( ! passed_lowres_filter || ! passed_highres_filter ) {
		set_last_move_status( protocols::moves::FAIL_RETRY );
	} else {
		set_last_move_status( protocols::moves::MS_SUCCESS );
	}

	// Check if pose is NOT full atom AND if the user wants full atom output, then do the right thing.
	if ( ! pose.is_fullatom() && basic::options::option[ basic::options::OptionKeys::out::file::fullatom ]() ) {
		// Convert pose to high resolution and recover sidechains
		// This is actually a sequnce mover that repacks after adding side chains -- do we always want that?
		to_all_atom_->apply( pose );
	}

	basic::prof_show();
}

//getters for const access to movers and data of docking protocol
protocols::simple_moves::SwitchResidueTypeSetMoverCOP DockingProtocol::to_centroid() const {
	return to_centroid_;
}
protocols::moves::MoverCOP DockingProtocol::to_all_atom() const{
	return to_all_atom_ ;
}
protocols::docking::DockingLowResCOP DockingProtocol::docking_lowres_mover() const{
	return docking_lowres_mover_;
}
protocols::docking::DockingHighResCOP DockingProtocol::docking_highres_mover() const{
	return docking_highres_mover_;
}
protocols::docking::DockingInitialPerturbationCOP DockingProtocol::perturber() const{
	return perturber_;
}

// Allow a developer to set a custom high resolution mover
void DockingProtocol::set_docking_highres_mover( protocols::docking::DockingHighResOP docking_highres_mover )
{
	docking_highres_mover_ = docking_highres_mover;

	// If docking_highres_mover_ has been configured with custom scorefxns, we will tell DockingProtocol.
	if ( docking_highres_mover_->scorefxn() && docking_highres_mover_->scorefxn_pack() ) {
		set_highres_scorefxn( docking_highres_mover_->scorefxn(), docking_highres_mover_->scorefxn_pack() );
	} else if ( docking_highres_mover_->scorefxn() && ! docking_highres_mover_->scorefxn_pack() ) {
		set_highres_scorefxn( docking_highres_mover_->scorefxn() );
	}
}

void DockingProtocol::add_additional_low_resolution_step( protocols::moves::MoverOP additional_low_resolution_mover )
{
	if ( ! additional_low_resolution_steps_ ) {
		additional_low_resolution_steps_ = protocols::moves::SequenceMoverOP( new moves::SequenceMover );
	}
	additional_low_resolution_steps_->add_mover( additional_low_resolution_mover );
}

/// @details  Show the complete setup of the docking protocol
void
DockingProtocol::show( std::ostream & out ) const {

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << ObjexxFCL::format::A( 47, "Rosetta 3 Docking Protocol" ) << ObjexxFCL::format::space( 27 ) << line_marker << std::endl;
	out << line_marker << ObjexxFCL::format::space( 74 ) << line_marker << std::endl;
	// Display the movable jumps that will be used in docking
	out << line_marker << " Dockable Jumps: ";

	int spaces_so_far = 23;

	bool first = true;
	for ( int movable_jump : movable_jumps_ ) {
		if ( !first ) {
			out << ", ";
			spaces_so_far += 2;
		} else first = false;

		out << ObjexxFCL::format::I( 1, movable_jump );
		spaces_so_far += 1;
	}

	int remaining_spaces = 80 - spaces_so_far;

	if ( remaining_spaces > 0 ) out << ObjexxFCL::format::space( 80 - spaces_so_far );
	out << line_marker << std::endl;

	// Display the state of the low resolution docking protocol that will be used
	out << line_marker << " Low Resolution Docking Protocol:  " << ( ( docking_lowres_mover_ ) ? ( " on " ) : ( "off " ) );
	out << ObjexxFCL::format::space( 35 ) << line_marker << std::endl;

	// Display the state of the low resolution docking protocol that will be used
	out << line_marker << " High Resolution Docking Protocol: " << ( ( docking_highres_mover_ ) ? ( " on " ) : ( "off " ) );
	out << ObjexxFCL::format::space( 35 ) << line_marker << std::endl;

	// Display the state of the filters (on or off)
	out << line_marker << " Low Resolution Filter:  " << ( ( lowres_filter_ ) ? ( " on " ) : ( "off " ) );
	out << ObjexxFCL::format::space( 45 ) << line_marker << std::endl;
	out << line_marker << " High Resolution Filter: " << ( ( highres_filter_ ) ? ( " on " ) : ( "off " ) );
	out << ObjexxFCL::format::space( 45 ) << line_marker << std::endl;

	// Close the box I have drawn
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
}

void
DockingProtocol::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	if ( tag->hasOption("docking_score_low" ) ) {
		std::string const score_low( tag->getOption<std::string>( "docking_score_low" ) );
		ScoreFunctionOP scorelo = data.get< ScoreFunction * >( "scorefxns", score_low )->clone();
		set_lowres_scorefxn(scorelo);
	}
	if ( tag->hasOption("docking_score_high" ) ) {
		std::string const score_high( tag->getOption<std::string>( "docking_score_high" ) );
		ScoreFunctionOP scorehi = data.get< ScoreFunction * >( "scorefxns", score_high )->clone();
		set_highres_scorefxn(scorehi);
	}
	//get through partners
	if ( tag->hasOption( "partners" ) ) {
		std::string const partners( tag->getOption<std::string>( "partners") );
		set_partners(partners);
	}
	//initialize other flags to control behavior
	//do high res step or not
	set_low_res_protocol_only( tag->getOption<bool>( "low_res_protocol_only", false ) );
	//skip the low res step if true
	set_docking_local_refine( tag->getOption<bool>( "docking_local_refine", false ) );
	//minimze final full atom structure?
	set_dock_min( tag->getOption<bool>( "dock_min", false ) );
	//ignore the default docking task and define your own
	set_ignore_default_docking_task( tag->getOption<bool>( "ignore_default_docking_task", false ) );
	if ( tag->hasOption( "task_operations" ) ) {
		set_task_factory(protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}

}//end parse_my_tag


// XRW TEMP std::string
// XRW TEMP DockingProtocolCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DockingProtocol::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DockingProtocolCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DockingProtocol() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DockingProtocol::mover_name()
// XRW TEMP {
// XRW TEMP  return "DockingProtocol";
// XRW TEMP }

std::string DockingProtocol::get_name() const {
	return mover_name();
}

std::string DockingProtocol::mover_name() {
	return "DockingProtocol";
}

void DockingProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "docking_score_low", xs_string, "Low-resolution docking score function" )
		+ XMLSchemaAttribute( "docking_score_high", xs_string, "High-resolution docking score function" )
		+ XMLSchemaAttribute( "partners", xs_string, "String specifying docking patners; underscore should separate the partners e.g. AB_C" )
		+ XMLSchemaAttribute::attribute_w_default( "low_res_protocol_only", xsct_rosetta_bool, "Only perform low-resolution docking", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "docking_local_refine", xsct_rosetta_bool, "Only perform high-resolution docking", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dock_min", xsct_rosetta_bool, "Use the DockMinMover", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_default_docking_task", xsct_rosetta_bool, "Ignore the default docking task and define your own. Unless this is specified, task operations will be appended to the default docking task.", "false" );
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Mover used for protein-protein docking.", attlist );
}

std::string DockingProtocolCreator::keyname() const {
	return DockingProtocol::mover_name();
}

protocols::moves::MoverOP
DockingProtocolCreator::create_mover() const {
	return protocols::moves::MoverOP( new DockingProtocol );
}

void DockingProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DockingProtocol::provide_xml_schema( xsd );
}


std::ostream & operator<<(std::ostream& out, const DockingProtocol & dp )
{
	dp.show( out );
	return out;
}

} //docking
} //protocols
