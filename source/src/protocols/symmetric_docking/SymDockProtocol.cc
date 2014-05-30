// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/symmetric_docking/SymDockProtocol.cc
///
/// @brief
/// @author Ingemar Andre


#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/symmetric_docking/SymSidechainMinMover.hh>
#include <protocols/symmetric_docking/SymDockProtocolCreator.hh>

////////////
#include <protocols/jd2/ScoreMap.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

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
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/util.hh>

#include <protocols/simple_filters/SAXSScoreFilter.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <core/pose/util.hh>
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
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/scoring/electron_density/util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/prof.hh>

#include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>


namespace protocols {
namespace symmetric_docking {

using namespace core;
using namespace ObjexxFCL;

static basic::Tracer TR("protocols.symmetric_docking.SymDockProtocol");

void SymDock_main() {
	using namespace protocols::simple_moves::symmetry;
  SetupForSymmetryMoverOP setup_mover = new SetupForSymmetryMover;
  SymDockProtocolOP dock_mover = new SymDockProtocol;
  protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
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

SymDockProtocol::~SymDockProtocol() {}

//clone
protocols::moves::MoverOP
SymDockProtocol::clone() const {
		return( new SymDockProtocol(  fullatom_, local_refine_, view_, docking_score_low_, docking_score_high_ ) );
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
	if (!fullatom_) passed_highres_filter_ = true;
}

void SymDockProtocol::set_local_refine( bool const local_refine_in )
{
	local_refine_ = local_refine_in;
	if (local_refine_){
		set_fullatom(true);
		passed_lowres_filter_ = true;
	}
}

void SymDockProtocol::set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_score_low_in )
{
	docking_score_low_ = docking_score_low_in;
	docking_low_ = 0;
}

void SymDockProtocol::set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_score_high_in )
{
	docking_score_high_ = docking_score_high_in;
	docking_score_high_min_ = docking_score_high_in;
	docking_high_ = 0;
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
	docking_score_pack_ = getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );

	if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ){
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
	/*	docking_score_low_ = ScoreFunctionFactory::create_score_function( "interchain_cen" ) ;
	docking_score_high_ = ScoreFunctionFactory::create_score_function( "docking" );
	docking_score_high_min_ = ScoreFunctionFactory::create_score_function( "docking", "docking_min" ) ;
	docking_score_pack_ = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) ;*/
	design(false); //??

	if ( dock_ppk_ ) set_local_refine(true);
	hurry( false );

	// score function setup
//	docking_score_low_ = ScoreFunctionFactory::create_score_function( "interchain_cen" ) ;
//	docking_score_high_ = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, DOCK_PATCH ) ;
/*	fullatom_ = option[ OptionKeys::out::file::fullatom ]();
	local_refine_ = option[ OptionKeys::docking::docking_local_refine ]();

	if (local_refine_) fullatom_=true;
	view_ = option[ OptionKeys::docking::view ]();

  // options
	protocol_ = "standard";
	docking_low_ = new protocols::symmetric_docking::SymDockingLowRes( scorefxn_lowres_ );
	docking_high_ = new protocols::symmetric_docking::SymDockingHiRes( scorefxn_hires_ );
*/
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
	if( basic::options::option[basic::options::OptionKeys::in::file::native].user() ){
		core::pose::PoseOP native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *native_pose, basic::options::option[basic::options::OptionKeys::in::file::native].value() );
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
	//using core::pose::datacache::CacheableDataType::SCORE_MAP;
	using utility::file::FileName;

	using namespace viewer;
	add_conformation_viewer( pose.conformation(), "start_pose", 450, 450 );

	//initialize docking protocol movers
	if (!docking_low_) docking_low_ = new SymDockingLowRes( docking_score_low_ );
	if (!docking_high_) docking_high_ = new SymDockingHiRes( docking_score_high_min_, docking_score_pack_ );

	// make sure the input pose has the right size
	core::pose::PoseOP input_pose = new core::pose::Pose();
	*input_pose = pose;
	docking_low_->set_input_pose( input_pose );
	docking_high_->set_input_pose( input_pose );
	set_input_pose( input_pose );


	if( init_task_factory_ ) {
		TR << "Setting non-default TaskFactory." << std::endl;
		docking_high_->task_factory( init_task_factory_ );
	}
	if( design_ ) {
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
	for (Size r = 1; r <= max_repeats; r++){
		pose = starting_pose;

	//MonteCarloOP mc;
		if ( !local_refine_ ) {
			docking_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_low_ ) ;

			// first convert to centroid mode
			to_centroid.apply( pose );
			if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ){
				protocols::simple_moves::ConstraintSetMover().apply(pose);
			}

			// make starting perturbations based on command-line flags over each movable jump (takes care of dock_pert and randomize)
			simple_moves::symmetry::SymDockingInitialPerturbation initial( true /*slide into contact*/ );
			initial.apply( pose );

			TR << "finished initial perturbation" << std::endl;

			if( !hurry_ && get_native_pose() ){
				Real st_rmsd = calc_rms( pose );
				score_map_["st_rmsd"] = st_rmsd;//jd1
				job->add_string_real_pair("st_rmsd", st_rmsd);
			}

			// low resolution docking

			TR.Debug << "DockingLowRes object created" << std::endl;
			docking_low_->apply( pose );

			//check low res filter to see if it should repeat low res or not
			if( !hurry_ ) passed_lowres_filter_ = docking_lowres_filter( pose );

			// add scores to map for output
			if( !hurry_ ) protocols::jd2::ScoreMap::nonzero_energies( score_map_, docking_scorefxn, pose );
			if( !hurry_ && get_native_pose() ){
				Real cen_rms = calc_rms( pose );
				score_map_["cen_rms"] = cen_rms; //jd1
				job->add_string_real_pair("cen_rms", cen_rms);
			}
		}

		// only do this is full atom is true
		if ( fullatom_ && passed_lowres_filter_ ) {

			docking_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *docking_score_high_ ) ;

			if (!local_refine_ || !pose.is_fullatom()){
				to_all_atom.apply( pose );
				if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ){
					protocols::simple_moves::ConstraintSetMover().apply(pose);
				}
				recover_sidechains.apply( pose );
			//	recover_sidechains( pose, *get_input_pose());
			}

			// run high resolution docking
			TR.Debug << "DockingHighRes object created" << std::endl;
			docking_high_->apply( pose );

			if( !hurry_ ) {
				Real interface_score = calc_interaction_energy( pose );
				// add scores to map for output
				protocols::jd2::ScoreMap::nonzero_energies( score_map_, docking_scorefxn, pose );//jd1
				score_map_["I_sc"] = interface_score; //jd1
				job->add_string_real_pair("I_sc", interface_score);
				// check highres docking filter
				if ( !option[ OptionKeys::docking::dock_ppk ]() ) {
					passed_highres_filter_ = docking_highres_filter( pose );
				}
			}

			if(passed_highres_filter_ && option[ OptionKeys::docking::kick_relax ].user()) {
				protocols::relax::relax_pose( pose, core::scoring::getScoreFunction(), "s" );
			}
		}

		if ( passed_lowres_filter_ && passed_highres_filter_ ) break; // defaults true. filter methods can set false
	  	else  TR<<"REPEAT STRUCTURE "<< r <<std::endl;
	}//end of repeat decoy creation

	if( !hurry_ ) {
		// calculate and store the rms no matter which mode was used
		if ( get_native_pose() ) {
			Real rms = calc_rms( pose );
			score_map_["rms"] = rms; //jd1
			job->add_string_real_pair("rms", rms);
		}

		if ( pose.is_fullatom() ){
			Real interface_score_ = calc_interaction_energy( pose );
			score_map_["I_sc"] = interface_score_; //jd1
			job->add_string_real_pair("I_sc", interface_score_);
		  if ( get_native_pose() ){
				//Real Irms = calc_rms(pose, *get_native_pose(), docking_score_high_, movable_jumps_ );
				//score_map_["Irms"] = Irms; //jd1
				//job->add_string_real_pair("Irms", Irms);
				//job->add_string_real_pair("Fnat", Fnat);
			}
		}else{
			score_map_["I_sc"] = calc_interaction_energy( pose );
		}
	}

	if ( option[ OptionKeys::run::debug ]() ) basic::prof_show();

	// cache the score map to the pose
	if( !hurry_ ) pose.data().set(core::pose::datacache::CacheableDataType::SCORE_MAP, new basic::datacache::DiagnosticData(score_map_));

}

bool
SymDockProtocol::docking_lowres_filter( core::pose::Pose & pose){

	using namespace core;
	using namespace core::scoring;
	using namespace basic::options;

	bool passed_filter = true;

	//criterion for failure
	Real interchain_contact_cutoff	= 10.0;
	Real interchain_vdw_cutoff = 1.0;
	Real distance_constraint_cutoff = 1.0; //distance constraint is soft for now

	if( option[ OptionKeys::docking::dock_lowres_filter ].user() ) {
		utility::vector1< Real > dock_filters = option[ OptionKeys::docking::dock_lowres_filter ]();
		interchain_contact_cutoff = dock_filters[1];
		interchain_vdw_cutoff = dock_filters[2];
		if (dock_filters.size() > 2) {
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

	if (icc_score >= interchain_contact_cutoff ) passed_filter = false;
	if (icvdw_score >= interchain_vdw_cutoff ) passed_filter = false;
	if (cst_score >= distance_constraint_cutoff ) passed_filter = false;

	if( ( option[basic::options::OptionKeys::filters::set_saxs_filter ].user() ) &&
	    ( option[basic::options::OptionKeys::score::saxs::ref_spectrum ].user() ) ) {
		protocols::simple_filters::SAXSScoreFilterOP saxs_filter = new protocols::simple_filters::SAXSScoreFilter();
    if( ! saxs_filter->apply(pose) )
			passed_filter = false;
		core::pose::setPoseExtraScores( pose, "saxs_score", saxs_filter->recent_score());
	}


	if (!passed_filter) {
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
	if (pose.energies().total_energy() >= score_cutoff) passed_filter = false;
	if (score_map_["I_sc"] >= 0.0) passed_filter = false;

	if (!passed_filter) TR << "STRUCTURE FAILED HIGH-RES FILTER " << pose.energies().total_energy() << " " << score_map_["I_sc"] << std::endl;

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

	// pack over each movable interface
		TaskFactoryOP tf = new TaskFactory;
		tf->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
		tf->push_back( new InitializeFromCommandline );
		tf->push_back( new IncludeCurrent );
		tf->push_back( new RestrictToRepacking );
		tf->push_back( new NoRepackDisulfides );
		if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) tf->push_back( new ReadResfile );
		//tf->push_back( new SymRestrictTaskForDocking( docking_score_pack_, true, 1000 ) );
		tf->push_back( new RestrictToInterface( 1 ) );

		protocols::simple_moves::PackRotamersMoverOP dock_pack = new protocols::simple_moves::symmetry::SymPackRotamersMover(docking_score_pack_);
		dock_pack->task_factory( tf );
		dock_pack->apply( pose );

		if (rtmin_){
	//		protocols::simple_moves::RotamerTrialsMinMoverOP rtmin_trial = new protocols::simple_moves::RotamerTrialsMinMover( docking_score_pack_, tf);
	//		rtmin_trial->apply( pose );
		}
		if (basic::options::option[ basic::options::OptionKeys::docking::sc_min ]()){
			SymInterfaceSidechainMinMoverOP scmin_mover = new SymInterfaceSidechainMinMover(docking_score_pack_, 1000);
			scmin_mover->apply(pose);
		}
}

std::string
SymDockProtocol::get_name() const {
	return "SymDockProtocol";
}


core::Real
SymDockProtocol::calc_interaction_energy( core::pose::Pose & pose ){
	using namespace scoring;
	using namespace core::conformation::symmetry;

	core::scoring::ScoreFunctionOP docking_scorefxn;
	core::pose::Pose complex_pose = pose;

	Real trans_magnitude = 1000;

	assert( core::pose::symmetry::is_symmetric( pose ) );
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	rigid::RigidBodyDofSeqTransMoverOP translate_away ( new rigid::RigidBodyDofSeqTransMover( dofs ) );
	translate_away->step_size( trans_magnitude );

	//Don't use patches for computer Isc, problematic with constraints for unbound
	if ( fullatom_ ){
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

	assert( core::pose::symmetry::is_symmetric( pose ) );
  SymmetricConformation & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

  FArray1D_bool superpos ( pose.total_residue(), false );
	for (Size res=1; res <= symm_info->num_total_residues_without_pseudo(); ++res )
	{
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

	if( tag->hasOption( "task_operations" ) ){
		task_factory(protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}

	if( tag->hasOption("docking_score_low" ) ){
		std::string const score_low( tag->getOption<std::string>( "docking_score_low" ) );
		set_lowres_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data, score_low ) );
	}
	if( tag->hasOption("docking_score_high" ) ){
		std::string const score_high( tag->getOption<std::string>( "docking_score_high" ) );
		set_highres_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data, score_high ) );
	}


	//initialize other flags to control behavior

	//void set_dock_rtmin( bool dock_rtmin_in );
	if( tag->hasOption("dock_rtmin" ) ){
		bool const dock_rtmin( tag->getOption<bool>( "dock_rtmin" ) );
		set_dock_rtmin(dock_rtmin);
	}


	//void set_sc_min( bool sc_min_in );
	if( tag->hasOption("sc_min" ) ){
		bool const sc_min( tag->getOption<bool>( "sc_min" ) );
		set_sc_min(sc_min);
	}


	//void set_max_repeats( Size const max_repeats_in );
	if( tag->hasOption("max_repeats" ) ){
		bool const max_repeats( tag->getOption<bool>( "max_repeats" ) );
		set_max_repeats(max_repeats);
	}


	//void set_dock_ppk( bool dock_ppk_in );
	if( tag->hasOption("dock_ppk" ) ){
		bool const dock_ppk( tag->getOption<bool>( "dock_ppk" ) );
		set_dock_ppk(dock_ppk);
	}


	//void set_fullatom( bool const fullatom_in );
	if( tag->hasOption("fullatom" ) ){
		bool const fullatom( tag->getOption<bool>( "fullatom" ) );
		set_fullatom(fullatom);
	}


	//void set_local_refine( bool const local_refine_in );
	if( tag->hasOption("local_refine" ) ){
		bool const local_refine( tag->getOption<bool>( "local_refine" ) );
		set_local_refine(local_refine);
	}


	//void set_view( bool view_in );
	if( tag->hasOption("view" ) ){
		bool const view( tag->getOption<bool>( "view" ) );
		set_view(view);
	}



}//end parse_my_tag



std::string
SymDockProtocolCreator::keyname() const
{
    return SymDockProtocolCreator::mover_name();
}

protocols::moves::MoverOP
SymDockProtocolCreator::create_mover() const {
    return new SymDockProtocol();
}

std::string
SymDockProtocolCreator::mover_name()
{
    return "SymDockProtocol";
}



} // symmetric_docking
} // protocols
