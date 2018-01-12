// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/symmetric_docking/SymDockProtocol.cc
///
/// @brief    Symmetric Protein-Protein Docking Protocol
/// @details  Dock symmetrical assemblies together - works with the symmetry
///           framework in Rosetta 3. Also includes some options for working with
///           membranes
///           Last Modified: 10/25/14
///
/// @author Ingemar Andre
/// @author Rebecca Alford (adding membranes & comments)

// Unit Headers
#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/symmetric_docking/SymSidechainMinMover.hh>
#include <protocols/symmetric_docking/SymDockProtocolCreator.hh>


#include <core/io/raw_data/ScoreMap.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>


#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/datacache/DiagnosticData.hh>

#include <core/types.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/symmetric_docking/SymDockingHiRes.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/util.hh>

#include <protocols/simple_filters/SAXSScoreFilter.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh> // PreventRepackingRLT
#include <core/pack/task/operation/ResFilters.hh> // ResidueLacksProperty
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>

#include <protocols/viewer/viewers.hh>
//for resfile reading
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/file/FileName.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <core/pose/symmetry/util.hh>


#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/scoring/electron_density/util.hh>


#include <core/import_pose/import_pose.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/prof.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace symmetric_docking {

using namespace core;
using namespace ObjexxFCL;

static basic::Tracer TR( "protocols.symmetric_docking.SymDockProtocol" );

void SymDock_main() {
	using namespace protocols::simple_moves::symmetry;
	SetupForSymmetryMoverOP setup_mover( new SetupForSymmetryMover );
	SymDockProtocolOP dock_mover( new SymDockProtocol );
	protocols::moves::SequenceMoverOP seq_mover( new protocols::moves::SequenceMover );
	seq_mover->add_mover( setup_mover );
	seq_mover->add_mover( dock_mover );
	protocols::jd2::JobDistributor::get_instance()->go( seq_mover );
}


SymDockProtocol::SymDockProtocol():
	Mover()
{
	Mover::type( "SymDockProtocol" );

	set_default();
	register_options();

}

SymDockProtocol::SymDockProtocol(
	bool const fullatom,
	bool const local_refine,
	bool const view
) : Mover()
{
	Mover::type( "SymDockProtocol" );
	set_default();
	set_fullatom(fullatom);
	set_local_refine(local_refine);
	set_view(view);
	register_options();
}

SymDockProtocol::SymDockProtocol(
	bool const fullatom,
	bool const local_refine,
	bool const view,
	core::scoring::ScoreFunctionOP docking_score_low,
	core::scoring::ScoreFunctionOP docking_score_high
) : Mover()
{
	Mover::type( "SymDockProtocol" );
	set_default();
	register_options();
	set_fullatom(fullatom);
	set_local_refine(local_refine);
	set_view(view);

	docking_score_low_ = docking_score_low;
	docking_score_high_ = docking_score_high;
	docking_score_high_min_ = docking_score_high;
	docking_score_pack_ = docking_score_high;
}

SymDockProtocol::~SymDockProtocol() = default;

//clone
protocols::moves::MoverOP
SymDockProtocol::clone() const {
	return( protocols::moves::MoverOP( new SymDockProtocol(  fullatom_, local_refine_, view_, docking_score_low_, docking_score_high_ ) ) );
}

//set functions
void SymDockProtocol::set_dock_rtmin( bool dock_rtmin_in ) { rtmin_=dock_rtmin_in; }
void SymDockProtocol::set_sc_min( bool sc_min_in ) { sc_min_=sc_min_in; }
void SymDockProtocol::set_max_repeats( Size const max_repeats_in ) { max_repeats_=max_repeats_in; }
void SymDockProtocol::set_dock_ppk( bool dock_ppk_in ) { dock_ppk_=dock_ppk_in; }
void SymDockProtocol::set_view( bool view_in ) { view_=view_in; }

void SymDockProtocol::set_fullatom( bool const fullatom_in )
{
	fullatom_ = fullatom_in;
	if ( !fullatom_ ) passed_highres_filter_ = true;
}

void SymDockProtocol::set_local_refine( bool const local_refine_in )
{
	local_refine_ = local_refine_in;
	if ( local_refine_ ) {
		set_fullatom(true);
		passed_lowres_filter_ = true;
	}
}

void SymDockProtocol::set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_score_low_in )
{
	docking_score_low_ = docking_score_low_in;
	docking_low_.reset();
}

void SymDockProtocol::set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_score_high_in )
{
	docking_score_high_ = docking_score_high_in;
	docking_score_high_min_ = docking_score_high_in;
	docking_high_.reset();
}

void SymDockProtocol::set_highres_scorefxn(
	core::scoring::ScoreFunctionOP docking_score_high_in,
	core::scoring::ScoreFunctionOP docking_score_pack_in )
{
	docking_score_high_ = docking_score_high_in;
	docking_score_high_min_ = docking_score_high_in;
	docking_score_pack_ = docking_score_pack_in;
}

void
SymDockProtocol::set_default()
{

	using namespace basic::options;
	using namespace core::scoring;

	passed_lowres_filter_ = true; // both default true. filter methods can set false
	passed_highres_filter_ = true;

	//how to turn on  atom_pair_constraint without patches by default to weight 1
	docking_score_low_  = ScoreFunctionFactory::create_score_function(  "interchain_cen" );
	docking_score_high_  = ScoreFunctionFactory::create_score_function( "docking" );
	docking_score_high_min_ = ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	docking_score_pack_ = get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );

	if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		docking_score_low_->set_weight(core::scoring::atom_pair_constraint, 1.0);
		docking_score_high_->set_weight(core::scoring::atom_pair_constraint, 1.0);
		docking_score_high_min_->set_weight(core::scoring::atom_pair_constraint, 1.0);
	}

	if ( option[ OptionKeys::docking::low_patch ].user() ) {
		//docking_score_low_  = ScoreFunctionFactory::create_score_function(  "interchain_cen", option[ OptionKeys::docking::low_patch ] );
		docking_score_low_->apply_patch_from_file(option[ OptionKeys::docking::low_patch ]);
	}

	if ( option[ OptionKeys::docking::high_patch ].user() ) {
		//docking_score_high_  = ScoreFunctionFactory::create_score_function( "docking", option[ OptionKeys::docking::high_patch ] );
		docking_score_high_->apply_patch_from_file(option[ OptionKeys::docking::high_patch ]);
	}

	if ( option[ OptionKeys::docking::high_min_patch ].user() ) {
		//docking_score_high_min_ = ScoreFunctionFactory::create_score_function( "docking",  option[ OptionKeys::docking::high_min_patch ] );
		docking_score_high_min_->apply_patch_from_file(option[ OptionKeys::docking::high_min_patch ]);
	}

	if ( option[ OptionKeys::docking::pack_patch ].user() ) {
		//docking_score_pack_ = ScoreFunctionFactory::create_score_function(  "standard", option[ OptionKeys::docking::pack_patch ] );
		docking_score_pack_->apply_patch_from_file(option[ OptionKeys::docking::pack_patch ]);
	}

	// score function setup
	design(false); //??

	dock_ppk_ = option[ OptionKeys::docking::dock_ppk ]();
	if ( dock_ppk_ ) set_local_refine(true);
	hurry( false );

	// add density if necessary
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *docking_score_low_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *docking_score_high_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *docking_score_high_min_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *docking_score_pack_ );
	}


}

void
SymDockProtocol::register_options()
{
	using namespace basic::options;

	set_fullatom(option[ OptionKeys::out::file::fullatom ]());
	set_local_refine(option[ OptionKeys::docking::docking_local_refine ]());

	set_dock_rtmin(option[ OptionKeys::docking::dock_rtmin ]());

	set_sc_min(option[ OptionKeys::docking::sc_min ]());
	set_max_repeats(option[ OptionKeys::docking::max_repeats ]());
	set_dock_ppk(option[ OptionKeys::docking::dock_ppk ]());

	//set native pose
	if ( basic::options::option[basic::options::OptionKeys::in::file::native].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_file( *native_pose, basic::options::option[basic::options::OptionKeys::in::file::native].value() , core::import_pose::PDB_file);
		this->set_native_pose( native_pose );
	}

}

// this is the docking protocol. Includes a low and high resolution stage. Very similar to heterodimeric docking protocol.
void
SymDockProtocol::apply( pose::Pose & pose )
{
	using namespace scoring;
	using namespace basic::options;
	using namespace moves;
	using utility::file::FileName;

	using namespace viewer;
	add_conformation_viewer( pose.conformation(), "start_pose", 450, 450 );

	//initialize docking protocol movers
	if ( !docking_low_ ) docking_low_ = protocols::symmetric_docking::SymDockingLowResOP( new SymDockingLowRes( docking_score_low_ ) );
	if ( !docking_high_ ) docking_high_ = protocols::symmetric_docking::SymDockingHiResOP( new SymDockingHiRes( docking_score_high_min_, docking_score_pack_ ) );

	// make sure the input pose has the right size
	core::pose::PoseOP input_pose( new core::pose::Pose() );
	*input_pose = pose;
	docking_low_->set_input_pose( input_pose );
	docking_high_->set_input_pose( input_pose );
	set_input_pose( input_pose );

	if ( init_task_factory_ ) {
		TR << "Setting non-default TaskFactory." << std::endl;
		docking_high_->task_factory( init_task_factory_ );
	}
	if ( design_ ) {
		TR << "Setting design to " << design_ << std::endl;
		docking_high_->design( design_ );
	}

	// Residue movers
	simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
	simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );
	simple_moves::ReturnSidechainMover recover_sidechains( *get_input_pose());

	core::scoring::ScoreFunctionOP docking_scorefxn;

	basic::prof_reset();

	core::pose::Pose starting_pose = pose;
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	if ( option[ OptionKeys::run::score_only ]() ) {
		score_only( pose );
		return;
	}

	core::Size const max_repeats( option[ OptionKeys::docking::max_repeats ]() );

	//start loop of decoy creation until filters are all passed
	for ( Size r = 1; r <= max_repeats; r++ ) {
		pose = starting_pose;

		//MonteCarloOP mc;
		if ( !local_refine_ ) {
			docking_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_low_ ) ;

			// first convert to centroid mode
			to_centroid.apply( pose );
			if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				protocols::simple_moves::ConstraintSetMover().apply(pose);
			}

			// make starting perturbations based on command-line flags over each movable jump (takes care of dock_pert and randomize)
			simple_moves::symmetry::SymDockingInitialPerturbation initial( true /*slide into contact*/ );
			initial.apply( pose );

			TR << "finished initial perturbation" << std::endl;
			if ( !hurry_ && get_native_pose() ) {
				Real st_rmsd = calc_rms( pose );
				score_map_["st_rmsd"] = st_rmsd;//jd1
				job->add_string_real_pair("st_rmsd", st_rmsd);
			}

			// low resolution docking

			TR.Debug << "DockingLowRes object created" << std::endl;
			docking_low_->apply( pose );

			//check low res filter to see if it should repeat low res or not
			if ( !hurry_ ) passed_lowres_filter_ = docking_lowres_filter( pose );

			// add scores to map for output
			if ( !hurry_ ) core::io::raw_data::ScoreMap::nonzero_energies( score_map_, docking_scorefxn, pose );
			if ( !hurry_ && get_native_pose() ) {
				Real cen_rms = calc_rms( pose );
				score_map_["cen_rms"] = cen_rms; //jd1
				job->add_string_real_pair("cen_rms", cen_rms);
			}
		}

		// only do this is full atom is true
		if ( fullatom_ && passed_lowres_filter_ ) {

			docking_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_high_ ) ;

			if ( !local_refine_ || !pose.is_fullatom() ) {
				to_all_atom.apply( pose );
				if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ) {
					protocols::simple_moves::ConstraintSetMover().apply(pose);
				}
				recover_sidechains.apply( pose );
				// recover_sidechains( pose, *get_input_pose());
			}

			// run high resolution docking
			TR.Debug << "DockingHighRes object created" << std::endl;
			docking_high_->apply( pose );

			if ( !hurry_ ) {
				Real interface_score = calc_interaction_energy( pose );
				// add scores to map for output
				core::io::raw_data::ScoreMap::nonzero_energies( score_map_, docking_scorefxn, pose );//jd1
				score_map_["I_sc"] = interface_score; //jd1
				job->add_string_real_pair("I_sc", interface_score);
				// check highres docking filter
				if ( !option[ OptionKeys::docking::dock_ppk ]() ) {
					passed_highres_filter_ = docking_highres_filter( pose );
				}
			}

			if ( passed_highres_filter_ && option[ OptionKeys::docking::kick_relax ].user() ) {
				protocols::relax::relax_pose( pose, core::scoring::get_score_function(), "s" );
			}
		}

		if ( passed_lowres_filter_ && passed_highres_filter_ ) break; // defaults true. filter methods can set false
		else  TR<<"REPEAT STRUCTURE "<< r <<std::endl;
	}//end of repeat decoy creation

	if ( !hurry_ ) {
		// calculate and store the rms no matter which mode was used
		if ( get_native_pose() ) {
			Real rms = calc_rms( pose );
			score_map_["rms"] = rms; //jd1
			job->add_string_real_pair("rms", rms);
		}

		if ( pose.is_fullatom() ) {
			Real interface_score_ = calc_interaction_energy( pose );
			score_map_["I_sc"] = interface_score_; //jd1
			job->add_string_real_pair("I_sc", interface_score_);
			if ( get_native_pose() ) {
				//Real Irms = calc_rms(pose, *get_native_pose(), docking_score_high_, movable_jumps_ );
				//score_map_["Irms"] = Irms; //jd1
				//job->add_string_real_pair("Irms", Irms);
				//job->add_string_real_pair("Fnat", Fnat);
			}
		} else {
			score_map_["I_sc"] = calc_interaction_energy( pose );
		}
	}

	if ( option[ OptionKeys::run::debug ]() ) basic::prof_show();

	// cache the score map to the pose
	using namespace basic::datacache;
	if ( !hurry_ ) pose.data().set(core::pose::datacache::CacheableDataType::SCORE_MAP, DataCache_CacheableData::DataOP( new basic::datacache::DiagnosticData(score_map_) ) );

}

bool
SymDockProtocol::docking_lowres_filter( core::pose::Pose & pose){

	using namespace core;
	using namespace core::scoring;
	using namespace basic::options;

	bool passed_filter = true;

	//criterion for failure
	Real interchain_contact_cutoff = 10.0;
	Real interchain_vdw_cutoff = 1.0;
	Real distance_constraint_cutoff = 1.0; //distance constraint is soft for now

	if ( option[ OptionKeys::docking::dock_lowres_filter ].user() ) {
		utility::vector1< Real > dock_filters = option[ OptionKeys::docking::dock_lowres_filter ]();
		interchain_contact_cutoff = dock_filters[1];
		interchain_vdw_cutoff = dock_filters[2];
		if ( dock_filters.size() > 2 ) {
			distance_constraint_cutoff = dock_filters[3];
		}
	}

	// fpd
	core::scoring::symmetry::SymmetricScoreFunction icc_sf, icvdw_sf, cst_sf;
	cst_sf.set_weight( core::scoring::atom_pair_constraint , 1.0 );
	icc_sf.set_weight( core::scoring::interchain_contact , 1.0 );
	icvdw_sf.set_weight( core::scoring::interchain_vdw , 1.0 );

	core::Real icc_score, icvdw_score, cst_score;
	icc_score = icc_sf(pose);
	icvdw_score = icvdw_sf(pose);
	cst_score = cst_sf(pose);

	if ( icc_score >= interchain_contact_cutoff ) passed_filter = false;
	if ( icvdw_score >= interchain_vdw_cutoff ) passed_filter = false;
	if ( cst_score >= distance_constraint_cutoff ) passed_filter = false;

	if ( ( option[basic::options::OptionKeys::filters::set_saxs_filter ].user() ) &&
			( option[basic::options::OptionKeys::score::saxs::ref_spectrum ].user() ) ) {
		protocols::simple_filters::SAXSScoreFilterOP saxs_filter( new protocols::simple_filters::SAXSScoreFilter() );
		if ( ! saxs_filter->apply(pose) ) {
			passed_filter = false;
		}
		core::pose::setPoseExtraScore( pose, "saxs_score", saxs_filter->recent_score());
	}


	if ( !passed_filter ) {
		TR << "STRUCTURE FAILED LOW-RES FILTER" << std::endl;
		//TR << " scores: " << pose.energies().total_energies()[ interchain_contact ] << "  " << pose.energies().total_energies()[ interchain_vdw ]
		//   << "  " << pose.energies().total_energies()[ atom_pair_constraint ] << std::endl;
		//TR << " cutoffs: " << interchain_contact_cutoff << "  " << interchain_vdw_cutoff << "  " << distance_constraint_cutoff << std::endl;
		TR << " scores: " << icc_score << "  " << icvdw_score << "  " << cst_score << std::endl;
		TR << " cutoffs: " << interchain_contact_cutoff << "  " << interchain_vdw_cutoff << "  " << distance_constraint_cutoff << std::endl;
	}

	return passed_filter;
}

bool
SymDockProtocol::docking_highres_filter( core::pose::Pose & pose){

	using namespace core;
	using namespace core::scoring;
	using namespace basic::options;

	bool passed_filter = true;
	if ( option[ OptionKeys::docking::dock_ppk ]() ) return passed_filter;

	Real score_cutoff = option[ OptionKeys::cluster::output_score_filter ]();
	//criterion for failure
	if ( pose.energies().total_energy() >= score_cutoff ) passed_filter = false;
	if ( score_map_["I_sc"] >= 0.0 ) passed_filter = false;

	if ( !passed_filter ) TR << "STRUCTURE FAILED HIGH-RES FILTER " << pose.energies().total_energy() << " " << score_map_["I_sc"] << std::endl;

	return passed_filter;
}

void
SymDockProtocol::recover_sidechains( core::pose::Pose & pose, const core::pose::Pose & native_pose )
{
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( native_pose );
	recover_sidechains.apply( pose );

	TR << "Doing initial repack of sidechains" << std::endl;

	//  Do initial pack over all residues within 1000A of the interface.

	using namespace moves;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	using core::pack::task::operation::TaskOperationOP;

	// pack over each movable interface
	TaskFactoryOP tf( new TaskFactory );
	tf->push_back( TaskOperationCOP( new OperateOnCertainResidues( ResLvlTaskOperationOP( new PreventRepackingRLT ), ResFilterOP( new ResidueLacksProperty("PROTEIN") ) ) ) );
	tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf->push_back( TaskOperationCOP( new IncludeCurrent ) );
	tf->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	tf->push_back( TaskOperationCOP( new NoRepackDisulfides ) );
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) tf->push_back( TaskOperationCOP( new ReadResfile ) );
	//tf->push_back( new SymRestrictTaskForDocking( docking_score_pack_, true, 1000 ) );
	tf->push_back( TaskOperationCOP( new RestrictToInterface( 1 ) ) );

	protocols::minimization_packing::PackRotamersMoverOP dock_pack( new protocols::minimization_packing::symmetry::SymPackRotamersMover(docking_score_pack_) );
	dock_pack->task_factory( tf );
	dock_pack->apply( pose );

	if ( rtmin_ ) {
		//  protocols::minimization_packing::RotamerTrialsMinMoverOP rtmin_trial = new protocols::minimization_packing::RotamerTrialsMinMover( docking_score_pack_, tf);
		//  rtmin_trial->apply( pose );
	}
	if ( basic::options::option[ basic::options::OptionKeys::docking::sc_min ]() ) {
		SymInterfaceSidechainMinMoverOP scmin_mover( new SymInterfaceSidechainMinMover(docking_score_pack_, 1000) );
		scmin_mover->apply(pose);
	}
}

// XRW TEMP std::string
// XRW TEMP SymDockProtocol::get_name() const {
// XRW TEMP  return "SymDockProtocol";
// XRW TEMP }


core::Real
SymDockProtocol::calc_interaction_energy( core::pose::Pose & pose ){
	using namespace scoring;
	using namespace core::conformation::symmetry;

	core::scoring::ScoreFunctionOP docking_scorefxn;
	core::pose::Pose complex_pose = pose;

	Real trans_magnitude = 1000;

	debug_assert( core::pose::symmetry::is_symmetric( pose ) );
	auto & symm_conf (
		dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	rigid::RigidBodyDofSeqTransMoverOP translate_away( new rigid::RigidBodyDofSeqTransMover( dofs ) );
	translate_away->step_size( trans_magnitude );

	//Don't use patches for computer Isc, problematic with constraints for unbound
	if ( fullatom_ ) {
		docking_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_pack_ ) ;
		docking_scorefxn->set_weight(core::scoring::atom_pair_constraint, 0.0);
		//docking_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "docking" );
	} else {
		docking_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_low_ ) ;
		docking_scorefxn->set_weight(core::scoring::atom_pair_constraint, 0.0);
		//docking_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
	}

	Real bound_energy = (*docking_scorefxn)( complex_pose );
	translate_away->apply( complex_pose );

	Real unbound_energy = (*docking_scorefxn)( complex_pose );
	return bound_energy - unbound_energy;

}

core::Real
SymDockProtocol::calc_rms( core::pose::Pose & pose ){

	using namespace core::conformation::symmetry;
	using namespace basic::options;

	debug_assert( core::pose::symmetry::is_symmetric( pose ) );
	auto & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	FArray1D_bool superpos ( pose.size(), false );
	for ( Size res=1; res <= symm_info->num_total_residues_without_pseudo(); ++res ) {
		superpos(res) = true;
	}
	if ( get_native_pose() ) {
		if ( option[ OptionKeys::symmetry::symmetric_rmsd ]() ) {
			return core::scoring::CA_rmsd_symmetric( *get_native_pose(), pose );
		} else {
			return core::scoring::rmsd_with_super_subset( *get_native_pose(), pose, superpos, core::scoring::is_protein_CA );
		}
	}
	return -1;
}

// turn on design of partner2 during docking. Not thoroughly tested!
void SymDockProtocol::design( bool const des ) { design_ = des; }
bool SymDockProtocol::design() const { return design_; }


void SymDockProtocol::hurry( bool const hurry ) { hurry_ = hurry; }

void
SymDockProtocol::task_factory( core::pack::task::TaskFactoryOP task )
{
	init_task_factory_ = task;
}

core::pack::task::TaskFactoryOP
SymDockProtocol::task_factory() const
{
	return( init_task_factory_ );
}

core::pack::task::TaskFactoryOP &
SymDockProtocol::task_factory()
{
	return( init_task_factory_ );
}


void
SymDockProtocol::score_only( core::pose::Pose & pose )
{
	using namespace scoring;
	using namespace moves;
	if ( fullatom_ ) {
		core::scoring::ScoreFunctionOP high_score = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_high_ );
		simple_moves::ScoreMover score_and_exit( high_score ) ;
		score_and_exit.insert_rms( calc_rms( pose) );
		score_and_exit.apply( pose );
	} else {
		simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( pose );
		core::scoring::ScoreFunctionOP low_score = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_low_ );
		simple_moves::ScoreMover score_and_exit( low_score );
		score_and_exit.insert_rms( calc_rms( pose ) );
		score_and_exit.apply( pose );
	}
}

void
SymDockProtocol::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	if ( tag->hasOption( "task_operations" ) ) {
		task_factory(protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}

	if ( tag->hasOption("docking_score_low" ) ) {
		std::string const score_low( tag->getOption<std::string>( "docking_score_low" ) );
		set_lowres_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data, score_low ) );
	}
	if ( tag->hasOption("docking_score_high" ) ) {
		std::string const score_high( tag->getOption<std::string>( "docking_score_high" ) );
		set_highres_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data, score_high ) );
	}


	//initialize other flags to control behavior

	//void set_dock_rtmin( bool dock_rtmin_in );
	if ( tag->hasOption("dock_rtmin" ) ) {
		bool const dock_rtmin( tag->getOption<bool>( "dock_rtmin" ) );
		set_dock_rtmin(dock_rtmin);
	}


	//void set_sc_min( bool sc_min_in );
	if ( tag->hasOption("sc_min" ) ) {
		bool const sc_min( tag->getOption<bool>( "sc_min" ) );
		set_sc_min(sc_min);
	}


	//void set_max_repeats( Size const max_repeats_in );
	if ( tag->hasOption("max_repeats" ) ) {
		bool const max_repeats( tag->getOption<bool>( "max_repeats" ) );
		set_max_repeats(max_repeats);
	}


	//void set_dock_ppk( bool dock_ppk_in );
	if ( tag->hasOption("dock_ppk" ) ) {
		bool const dock_ppk( tag->getOption<bool>( "dock_ppk" ) );
		set_dock_ppk(dock_ppk);
	}


	//void set_fullatom( bool const fullatom_in );
	if ( tag->hasOption("fullatom" ) ) {
		bool const fullatom( tag->getOption<bool>( "fullatom" ) );
		set_fullatom(fullatom);
	}


	//void set_local_refine( bool const local_refine_in );
	if ( tag->hasOption("local_refine" ) ) {
		bool const local_refine( tag->getOption<bool>( "local_refine" ) );
		set_local_refine(local_refine);
	}


	//void set_view( bool view_in );
	if ( tag->hasOption("view" ) ) {
		bool const view( tag->getOption<bool>( "view" ) );
		set_view(view);
	}


}//end parse_my_tag


// XRW TEMP std::string
// XRW TEMP SymDockProtocolCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SymDockProtocol::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SymDockProtocolCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SymDockProtocol() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SymDockProtocol::mover_name()
// XRW TEMP {
// XRW TEMP  return "SymDockProtocol";
// XRW TEMP }

std::string SymDockProtocol::get_name() const {
	return mover_name();
}

std::string SymDockProtocol::mover_name() {
	return "SymDockProtocol";
}

void SymDockProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "dock_rtmin", xsct_rosetta_bool, "Does rotamer trials with minimization." )
		+ XMLSchemaAttribute( "sc_min", xsct_rosetta_bool, "Does sidechain minimization of interface residues." )
		+ XMLSchemaAttribute( "max_repeats", xsct_rosetta_bool, "If a decoy does not pass the low- and high-resolution filters, how many attempts to make before failure." )
		+ XMLSchemaAttribute( "dock_ppk", xsct_rosetta_bool, "Docking prepack mode." )
		+ XMLSchemaAttribute( "fullatom", xsct_rosetta_bool, "Enable full-atom input of PDB or centroid structures." )
		+ XMLSchemaAttribute( "local_refine", xsct_rosetta_bool, "Do a local refinement of the docking position (high resolution)." )
		+ XMLSchemaAttribute( "view", xsct_rosetta_bool, "Decide whether to use the viewer (graphical) or not." ) ;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist , "docking_score_low" ) ;
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist , "docking_score_high" ) ;
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Symmetric oligomer docking.", attlist );
}

std::string SymDockProtocolCreator::keyname() const {
	return SymDockProtocol::mover_name();
}

protocols::moves::MoverOP
SymDockProtocolCreator::create_mover() const {
	return protocols::moves::MoverOP( new SymDockProtocol );
}

void SymDockProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymDockProtocol::provide_xml_schema( xsd );
}



} // symmetric_docking
} // protocols
