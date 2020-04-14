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
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>


#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/kinematics/FoldTree.hh>

#include <basic/datacache/DiagnosticData.hh>

#include <core/types.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <protocols/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/symmetric_docking/SymDockingHiRes.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
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
#include <protocols/simple_task_operations/RestrictToInterface.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/scoring/Interface.hh>

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
#include <protocols/symmetry/SetupForSymmetryMover.hh>
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
	using namespace protocols::symmetry;
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
	return( utility::pointer::make_shared< SymDockProtocol >(  fullatom_, local_refine_, view_, docking_score_low_, docking_score_high_ ) );
}

//set functions
void SymDockProtocol::set_dock_rtmin( bool dock_rtmin_in ) { rtmin_=dock_rtmin_in; }
void SymDockProtocol::set_sc_min( bool sc_min_in ) { sc_min_=sc_min_in; }
void SymDockProtocol::set_max_repeats( core::Size const max_repeats_in ) { max_repeats_=max_repeats_in; }
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

	// Ideally would like to have motif dock score as default low-res score function,
	// but cannot make it default as score tables need to be imported from Git-LFS
	// and it won't work without it, so keeping interchain_cen as placeholder.

	if ( option[ OptionKeys::docking::docking_low_res_score ].user() ) {
		docking_score_low_ = ScoreFunctionFactory::create_score_function( option[ OptionKeys::docking::docking_low_res_score ]() );
	} else if ( docking_score_low_ == nullptr ) {
		docking_score_low_ = ScoreFunctionFactory::create_score_function( "interchain_cen", option[ OptionKeys::score::patch ]() );
	}

	docking_score_high_  = core::scoring::get_score_function();
	docking_score_high_min_ = core::scoring::get_score_function();
	docking_score_pack_ = core::scoring::get_score_function();

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
	using namespace core::conformation::symmetry;
	using namespace core::scoring;

	using namespace viewer;
	add_conformation_viewer( pose.conformation(), "start_pose", 450, 450 );

	//initialize docking protocol movers
	if ( !docking_low_ ) docking_low_ = utility::pointer::make_shared< SymDockingLowRes >( docking_score_low_ );

	docking_score_high_min_->set_weight( core::scoring::fa_rep,  basic::options::option[ basic::options::OptionKeys::docking::SymDock_fa_rep_max ]() );
	docking_score_pack_->set_weight( core::scoring::fa_rep,  basic::options::option[ basic::options::OptionKeys::docking::SymDock_fa_rep_max ]() );
	docking_score_high_min_->set_weight( core::scoring::fa_sol,  basic::options::option[ basic::options::OptionKeys::docking::SymDock_fa_sol_max ]() );
	docking_score_pack_->set_weight( core::scoring::fa_sol,  basic::options::option[ basic::options::OptionKeys::docking::SymDock_fa_sol_max ]() );

	if ( !docking_high_ ) docking_high_ = utility::pointer::make_shared< SymDockingHiRes >( docking_score_high_min_, docking_score_pack_ );

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



	debug_assert( core::pose::symmetry::is_symmetric( pose ) );

	SymmetricConformation const & symm_conf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()) );

	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size const N ( symm_info->subunits() );

	if ( option[ OptionKeys::run::score_only ]() ) {
		score_only( pose );
		return;
	}

	core::Size const max_repeats( option[ OptionKeys::docking::max_repeats ]() );

	//start loop of decoy creation until filters are all passed
	for ( core::Size r = 1; r <= max_repeats; r++ ) {
		pose = starting_pose;

		//MonteCarloOP mc;
		if ( !local_refine_ ) {

			docking_scorefxn = docking_score_low_->clone();
			if ( option[ OptionKeys::docking::docking_low_res_score ]()=="motif_dock_score" && option[ OptionKeys::docking::SymDock_reduce_motif_dock_weights ]() ) {
				TR << "Setting weight of the score function to " << (Real) 1/N << " times the default value." << std::endl;
				docking_scorefxn->set_weight( score_type_from_name("motif_dock"), (Real) 1/N );
			}
			std::map < std::string, core::Real > lowres_scores;

			// first convert to centroid mode
			to_centroid.apply( pose );
			if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				protocols::constraint_movers::ConstraintSetMover().apply(pose);
			}

			// make starting perturbations based on command-line flags over each movable jump (takes care of dock_pert and randomize)
			symmetry::SymDockingInitialPerturbation initial( true /*slide into contact*/ );
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
			if ( !hurry_ ) {
				docking_scorefxn->score(pose);
				core::io::raw_data::ScoreMap::add_energies_data_from_scored_pose( pose,  score_map_ );
			}
			if ( !hurry_ && get_native_pose() ) {
				Real cen_rms = calc_rms( pose );
				score_map_["cen_rms"] = cen_rms; //jd1
				job->add_string_real_pair("cen_rms", cen_rms);
			}

			// store the low res scores for output
			if ( !hurry_ ) {
				core::io::raw_data::ScoreMap::add_energies_data_from_scored_pose( pose, lowres_scores );

				// output low-res terms in score file if low res was run
				if ( !lowres_scores.empty() ) {
					for ( std::map< std::string, core::Real >::const_iterator pair=lowres_scores.begin(); pair!=lowres_scores.end(); ++pair ) {
						if ( pair->first == "total_score" ) {
							// because total_score should be the high-res one
							job->add_string_real_pair( "total_lowres_score", pair->second );
						} else {
							job->add_string_real_pair( pair->first, pair->second );
						}
					}
				}
			}
		}

		// only do this is full atom is true
		if ( fullatom_ && passed_lowres_filter_ ) {

			docking_scorefxn = docking_score_high_->clone();

			if ( !local_refine_ || !pose.is_fullatom() ) {
				to_all_atom.apply( pose );
				if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_fa_file ].user() ) {
					protocols::constraint_movers::ConstraintSetMover().apply(pose);
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
				docking_scorefxn->score(pose);
				core::io::raw_data::ScoreMap::add_energies_data_from_scored_pose( pose, score_map_ );//jd1
				score_map_["I_sc"] = interface_score; //jd1
				job->add_string_real_pair("I_sc", interface_score);
				// check highres docking filter
				if ( !option[ OptionKeys::docking::dock_ppk ]() ) {
					passed_highres_filter_ = docking_highres_filter( pose );
				}
			}

			if ( passed_highres_filter_ && option[ OptionKeys::docking::kick_relax ].value() ) {
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
			Real Irms = calc_Irms( pose );
			score_map_["Irms"] = Irms;
			job->add_string_real_pair("Irms", Irms);
			Real Fnat = calc_fnat( pose, docking_score_pack_ );
			score_map_["Fnat"] = Fnat;
			job->add_string_real_pair("Fnat", Fnat);
			Real CAPRI_rank = calc_CAPRI_rank( Irms, rms, Fnat );
			score_map_["CAPRI_rank"] = CAPRI_rank;
			job->add_string_real_pair("CAPRI_rank", CAPRI_rank);
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
	if ( !hurry_ ) pose.data().set(core::pose::datacache::CacheableDataType::SCORE_MAP, utility::pointer::make_shared< basic::datacache::DiagnosticData >(score_map_) );

}

bool
SymDockProtocol::docking_lowres_filter( core::pose::Pose & pose){

	using namespace core;
	using namespace core::scoring;
	using namespace basic::options;

	bool passed_filter = true;

	using namespace core::conformation::symmetry;
	SymmetricConformation const & symm_conf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()) );

	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size const N ( symm_info->subunits() );

	//criterion for failure
	Real interchain_vdw_cutoff = 1.0;
	Real interchain_contact_cutoff = 10.0;
	Real distance_constraint_cutoff = 1.0; //distance constraint is soft for now

	//MDS requires only one filter and we set it to a looser bound

	if ( option[ OptionKeys::docking::docking_low_res_score ]()=="motif_dock_score" ) {
		Real base_interchain_vdw_cutoff = option[ OptionKeys::docking::SymDock_lowres_filter ]();
		TR << "Motif Dock Score: setting clash filter to "<< base_interchain_vdw_cutoff << " x " << N << " (number of interfaces per subunit)." << std::endl;
		interchain_vdw_cutoff = base_interchain_vdw_cutoff * N;
	}

	if ( option[ OptionKeys::docking::dock_lowres_filter ].user() ) {
		utility::vector1< Real > dock_filters = option[ OptionKeys::docking::dock_lowres_filter ]();
		interchain_contact_cutoff = dock_filters[1];
		interchain_vdw_cutoff = dock_filters[2];
		if ( dock_filters.size() > 2 ) {
			distance_constraint_cutoff = dock_filters[3];
		}
	}

	// fpd
	core::scoring::ScoreFunction icc_sf, icvdw_sf, cst_sf;
	cst_sf.set_weight( core::scoring::atom_pair_constraint , 1.0 );
	icc_sf.set_weight( core::scoring::interchain_contact , 1.0 );
	icvdw_sf.set_weight( core::scoring::interchain_vdw , 1.0 );

	core::Real icc_score, icvdw_score, cst_score;
	icc_score = icc_sf(pose);
	icvdw_score = icvdw_sf(pose);
	cst_score = cst_sf(pose);


	// For MDS we only care about interchain_vdw filter

	if ( option[ OptionKeys::docking::docking_low_res_score ]()=="motif_dock_score" ) {
		if ( interchain_vdw_cutoff < 0 ) {
			// a negative filter means no lowres filter
			passed_filter = true;
		} else {
			if ( icvdw_score >= interchain_vdw_cutoff ) passed_filter = false;
		}

	} else {

		if ( icc_score >= interchain_contact_cutoff ) passed_filter = false;
		if ( icvdw_score >= interchain_vdw_cutoff ) passed_filter = false;
		if ( cst_score >= distance_constraint_cutoff ) passed_filter = false;
	}

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
	using namespace protocols::simple_task_operations;
	using core::pack::task::operation::TaskOperationOP;

	// pack over each movable interface
	TaskFactoryOP tf( new TaskFactory );
	tf->push_back( utility::pointer::make_shared< OperateOnCertainResidues >( utility::pointer::make_shared< PreventRepackingRLT >(), utility::pointer::make_shared< ResidueLacksProperty >("PROTEIN") ) );
	tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >() );
	tf->push_back( utility::pointer::make_shared< IncludeCurrent >() );
	tf->push_back( utility::pointer::make_shared< RestrictToRepacking >() );
	tf->push_back( utility::pointer::make_shared< NoRepackDisulfides >() );
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) tf->push_back( utility::pointer::make_shared< ReadResfile >() );
	//tf->push_back( new SymRestrictTaskForDocking( docking_score_pack_, true, 1000 ) );
	tf->push_back( utility::pointer::make_shared< RestrictToInterface >( 1 ) );

	protocols::minimization_packing::PackRotamersMoverOP dock_pack( new protocols::minimization_packing::PackRotamersMover(docking_score_pack_) );
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

	std::map< core::Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	rigid::RigidBodyDofSeqTransMoverOP translate_away( new rigid::RigidBodyDofSeqTransMover( dofs ) );
	translate_away->step_size( trans_magnitude );

	//Don't use patches for computer Isc, problematic with constraints for unbound
	if ( fullatom_ ) {
		docking_scorefxn = docking_score_pack_->clone();
		docking_scorefxn->set_weight(core::scoring::atom_pair_constraint, 0.0);
		//docking_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "docking" );
	} else {
		docking_scorefxn = docking_score_low_->clone();
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
	for ( core::Size res=1; res <= symm_info->num_total_residues_without_pseudo(); ++res ) {
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



core::Real
SymDockProtocol::calc_Irms( core::pose::Pose & pose ){

	using namespace core::conformation::symmetry;
	using namespace basic::options;

	debug_assert( core::pose::symmetry::is_symmetric( pose ) );

	SymmetricConformation const & symm_conf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()) );

	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size const N ( symm_info->subunits() );
	Size const nres_per_monomer ( symm_info->num_total_residues_without_pseudo()/N );

	using namespace scoring;
	Real Irmsd( 1e3 ); // starting value

	core::pose::Pose native_docking_pose = *get_native_pose()->clone();

	debug_assert( native_docking_pose.size() == pose.size() - N );// one is symmterized pose with N virtual atoms

	if ( N > 2 ) {

		using namespace core::kinematics;
		auto const cutpoints = native_docking_pose.fold_tree().cutpoints();
		Size chainB_start = cutpoints[1] + 1;

		FoldTree ft;

		ft.add_edge(1, cutpoints[1], -1);
		ft.add_edge(1, chainB_start, 1);
		ft.add_edge(chainB_start, cutpoints[2], -1);

		for ( Size i = 2; i <= cutpoints.size(); ++i ) {
			ft.add_edge(chainB_start, cutpoints[i]+1, i);
			if ( i != cutpoints.size() ) {
				ft.add_edge( cutpoints[i] + 1, cutpoints[i+1], -1);
			} else {
				ft.add_edge( cutpoints[i] + 1, native_docking_pose.size(), -1);
			}
		}

		ft.check_fold_tree();

		TR << "Restructing native foldtree for interface detection of first chain. New foldtree:" << std::endl;
		TR << ft << std::endl;

		native_docking_pose.fold_tree( ft );
	} // no need to resturcture for C2

	// superimposing the first chains of pose and native_dockng_pose

	using namespace core::id;
	AtomID_Map< AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, native_docking_pose, AtomID::BOGUS_ATOM_ID() );

	for ( Size i=1; i<=nres_per_monomer; ++i ) {
		AtomID const id1( native_docking_pose.residue(i).atom_index("CA"), i );
		AtomID const id2( pose.residue(i).atom_index("CA"), i );
		atom_map[ id1 ] = id2;
	}


	using namespace core::scoring;


	utility::vector1< utility::vector1<Size> > individual_rmsds(N - 1);

	for ( Size i=2; i<=N; ++i ) {

		utility::vector1<Size> rmsd_chain_temp;
		for ( Size j=2; j<=N; ++j ) {

			utility::vector1< Size > native_residues;
			utility::vector1< Size > pose_residues;

			for ( Size k=1; k<=nres_per_monomer; ++k ) {
				native_residues.push_back(nres_per_monomer * (i-1) + k);
				pose_residues.push_back(nres_per_monomer * (j-1) + k);
			}

			rmsd_chain_temp.push_back( rmsd_no_super( native_docking_pose, pose, native_residues, pose_residues, is_protein_CA ) );
		}
		individual_rmsds[i-1] = rmsd_chain_temp;
	}

	// Verifying individual_rmsds
	for ( Size i = 1; i <= individual_rmsds.size(); ++i ) {

		TR.Debug << "RMSDs of all permutations:  ";

		for ( Size j = 1; j <= individual_rmsds[i].size(); ++j ) {

			TR.Debug << individual_rmsds[i][j] << "  ";
		}
		TR.Debug << std::endl;
	}

	// Obtaining optimal chain ordering based on lowest RMSD combinations

	std::map< Size, Size > chain_order;
	chain_order[ 1 ] = 1;

	utility::vector1<Size> native_excluded_chains;
	utility::vector1<Size> pose_excluded_chains;

	Real min_rmsd( 1e3 ); // starting value
	Size native_chain = 1;
	Size pose_chain = 1;

	for ( Size k=1; k<=N-1; ++k ) {

		for ( Size row=1; row<=N-1; ++row ) {

			if ( !pose_excluded_chains.has_value(row) ) {

				for ( Size col=1; col<=N-1; ++col ) {
					if ( !native_excluded_chains.has_value(col) ) {

						Real rmsd_temp = individual_rmsds[row][col];
						if ( rmsd_temp < min_rmsd ) {
							min_rmsd = rmsd_temp;
							pose_chain = row;
							native_chain = col;
						}
					}
				}
			}
		}

		// identified the lowest remaining value in the table, excluding these chains now
		chain_order[ native_chain + 1 ] = pose_chain + 1;
		native_excluded_chains.push_back(native_chain);
		pose_excluded_chains.push_back(pose_chain);
		min_rmsd = 1e3;
		native_chain = 1;
		pose_chain = 1;
	}


	// Arranging chains based on optimal ordering and calculating the rmsd
	AtomID_Map< AtomID > atom_map_final;
	core::pose::initialize_atomid_map( atom_map_final, native_docking_pose, AtomID::BOGUS_ATOM_ID() );

	for ( Size i=1; i<=N; ++i ) {
		Size j = chain_order[i];
		for ( Size k=1; k<=nres_per_monomer; ++k ) {

			TR.Debug << "Native: " << i << "\t" << k << "\t" << nres_per_monomer * (i-1) + k << std::endl;
			AtomID const id1( native_docking_pose.residue(nres_per_monomer * (i-1) + k).atom_index("CA"), nres_per_monomer * (i-1) + k );
			TR.Debug << "Pose: " << j << "\t" << k << "\t" << nres_per_monomer * (j-1) + k << std::endl;
			AtomID const id2( pose.residue(nres_per_monomer * (j-1) + k).atom_index("CA"), nres_per_monomer * (j-1) + k );
			atom_map_final[ id1 ] = id2;
		}
	}

	Irmsd = superimpose_pose( native_docking_pose, pose, atom_map_final );

	return Irmsd;


}


core::Real
SymDockProtocol::calc_fnat( core::pose::Pose & pose, core::scoring::ScoreFunctionOP dock_scorefxn ){

	using namespace core::conformation::symmetry;
	using namespace basic::options;
	using namespace core::conformation;

	debug_assert( core::pose::symmetry::is_symmetric( pose ) );

	SymmetricConformation const & symm_conf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()) );

	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size const N ( symm_info->subunits() );
	Size const nres_per_monomer ( symm_info->num_total_residues_without_pseudo()/N );

	using namespace scoring;
	Real fnat( 0 ); // starting value
	Real cutoff = 5.0;

	core::pose::Pose native_docking_pose = *get_native_pose()->clone();


	debug_assert( native_docking_pose.size() == pose.size() - N );// one is symmterized pose with N virtual atoms

	if ( N > 2 ) {

		using namespace core::kinematics;
		auto const cutpoints = native_docking_pose.fold_tree().cutpoints();
		Size chainB_start = cutpoints[1] + 1;

		FoldTree ft;

		ft.add_edge(1, cutpoints[1], -1);
		ft.add_edge(1, chainB_start, 1);
		ft.add_edge(chainB_start, cutpoints[2], -1);

		for ( Size i = 2; i <= cutpoints.size(); ++i ) {
			ft.add_edge(chainB_start, cutpoints[i]+1, i);
			if ( i != cutpoints.size() ) {
				ft.add_edge( cutpoints[i] + 1, cutpoints[i+1], -1);
			} else {
				ft.add_edge( cutpoints[i] + 1, native_docking_pose.size(), -1);
			}
		}

		ft.check_fold_tree();

		TR << "Restructing native foldtree for interface detection of first chain. New foldtree:" << std::endl;
		TR << ft << std::endl;

		native_docking_pose.fold_tree( ft );
	}

	// score to set up interface object
	// scoring only happened here to update the residue neighbors
	core::scoring::ScoreFunctionOP scorefxn = dock_scorefxn->clone();

	( *scorefxn )( native_docking_pose );

	protocols::scoring::Interface interface( 1 ); //only considering interfaces A makes
	interface.distance( 8.0 );
	interface.calculate( native_docking_pose );

	utility::vector1< utility::vector1<Size> > native_interface_residues(N); // empty 2-D vector array, each row is for one subunit

	Size chain_number = 0;
	for ( Size i = 1; i <= native_docking_pose.size(); ++i ) {
		if ( interface.is_interface( i ) ) {
			chain_number = native_docking_pose.chain(i);
			native_interface_residues[chain_number].push_back(i);
		}
	}

	// creating contact pair list (not using ObjexxFCL::FArray3D_bool as we need variable sized arrays)
	// populating contact list with falses
	utility::vector1< utility::vector1<utility::vector1<bool>> > contact_list( N-1 );

	for ( Size i = 2; i <= N; ++i ) {
		utility::vector1<bool> temp_interface_list( native_interface_residues[i].size() );
		utility::vector1< utility::vector1<bool> > temp_2D_bool_list;
		for ( Size j = 1; j <= native_interface_residues[1].size(); ++j ) {
			utility::vector1<bool> temp_bool_list = temp_interface_list; // deep copying
			temp_2D_bool_list.push_back(temp_bool_list);
		}
		contact_list[i-1] = temp_2D_bool_list;
		temp_2D_bool_list.clear();
		temp_interface_list.clear();
	}

	//identify native contacts across the interface
	utility::vector1< utility::vector1<utility::vector1<bool>> > native_contact_list = contact_list;

	for ( Size i = 2; i<= N; ++i ) {
		for ( Size j = 1; j <= native_interface_residues[1].size(); ++j ) {
			ResidueOP rsd1( new Residue( native_docking_pose.residue( native_interface_residues[1][j] ) ) );
			for ( Size k = 1; k <= native_interface_residues[i].size(); ++k ) {
				ResidueOP rsd2( new Residue( native_docking_pose.residue( native_interface_residues[i][k] ) ) );
				native_contact_list[i-1][j][k] = calc_res_contact( rsd1, rsd2, cutoff ) ;
			}
		}
	}

	// manually calculating permutations for all chains except the very first chain
	// then calculating rmsd with superposition and taking the minimum to get correct chain orientation

	std::string secondary_chains;
	for ( Size i = 1; i < N; ++i ) {
		secondary_chains += std::to_string(i);
	}

	utility::vector1< std::string > chain_permut;
	do {
		chain_permut.push_back(secondary_chains);
	} while(std::next_permutation(secondary_chains.begin(), secondary_chains.end()));

	// finding best permutation of chains to get highest Fnat
	for ( auto & chain_order : chain_permut ) {

		utility::vector1< utility::vector1<utility::vector1<bool>> > pose_contact_list = contact_list;
		Real native_ncontact = 0;
		Real pose_ncontact = 0;

		for ( Size i = 2; i<= N; ++i ) {
			for ( Size j = 1; j <= native_interface_residues[1].size(); ++j ) {
				ResidueOP rsd1( new Residue( pose.residue( native_interface_residues[1][j] ) ) );
				for ( Size k = 1; k <= native_interface_residues[i].size(); ++k ) {
					Size orig_res_number = native_interface_residues[i][k];
					Size chain_number = orig_res_number / nres_per_monomer;
					if ( orig_res_number % nres_per_monomer == 0 ) chain_number = chain_number - 1; // this happens for the C-terminus residue
					Size chain_number_on_first = orig_res_number - ( chain_number * nres_per_monomer );
					Size new_residue_number = ( (int) chain_order[ chain_number - 1 ] - '0' ) * nres_per_monomer+ chain_number_on_first ;
					ResidueOP rsd2( new Residue( pose.residue( new_residue_number ) ) );
					pose_contact_list[i-1][j][k] = calc_res_contact( rsd1, rsd2, cutoff ) ;
				}
			}
		}

		//identify which native contacts are recovered in the decoy
		for ( Size i = 1; i<= native_contact_list.size(); ++i ) {
			for ( Size j = 1; j<= native_contact_list[i].size(); ++j ) {
				for ( Size k = 1; k<= native_contact_list[i][j].size(); ++k ) {
					if ( native_contact_list[i][j][k] ) {
						++native_ncontact;
						if ( pose_contact_list[i][j][k] ) {
							++pose_ncontact;
						}
					}

				}
			}
		}
		if ( fnat < pose_ncontact/native_ncontact ) fnat = pose_ncontact/native_ncontact;
	}

	return fnat;

}

bool SymDockProtocol::calc_res_contact(
	conformation::ResidueOP const rsd1,
	conformation::ResidueOP const rsd2,
	Real const dist_cutoff ) {
	Real dist_cutoff2 = dist_cutoff * dist_cutoff;
	double dist2 = 9999.0;

	for ( Size m = 1; m <= rsd1->nheavyatoms(); m++ ) {
		for ( Size n = 1; n <= rsd2->nheavyatoms(); n++ ) {
			dist2 = rsd1->xyz( m ).distance_squared( rsd2->xyz( n ) );
			if ( dist2 <= dist_cutoff2 ) {
				return true;
			}
		}
	}
	return false;
}

core::Real
SymDockProtocol::calc_CAPRI_rank( Real const Irmsd, Real const Lrmsd, Real const Fnat ) {

	if ( (Fnat >= 0.5 ) && ( Lrmsd <= 1 || Irmsd <= 1 ) ) {
		return 3; // high quality
	} else if ( (Fnat >= 0.5 && Lrmsd > 1 && Irmsd > 1) || (Fnat >= 0.3 && (Lrmsd <= 5 || Irmsd <= 2)) ) {
		return 2; // medium quality
	} else if ( (Fnat >= 0.3 && Lrmsd > 5 && Irmsd > 2) || (Fnat >= 0.1 && (Lrmsd <= 10 || Irmsd <= 4)) ) {
		return 1; // acceptable quality
	} else if ( (Fnat < 0.1) || (Lrmsd > 10 && Irmsd > 4) ) {
		return 0; // incorrect structure
	} else {
		TR.Debug << "Unreachable code" << std::endl;
		return -1; // should never be reached
	}

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
		core::scoring::ScoreFunctionOP high_score = docking_score_high_->clone();
		simple_moves::ScoreMover score_and_exit( high_score ) ;
		score_and_exit.insert_rms( calc_rms( pose) );
		score_and_exit.apply( pose );
	} else {
		simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( pose );
		core::scoring::ScoreFunctionOP low_score = docking_score_low_->clone();
		simple_moves::ScoreMover score_and_exit( low_score );
		score_and_exit.insert_rms( calc_rms( pose ) );
		score_and_exit.apply( pose );
	}
}

void
SymDockProtocol::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
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


	//void set_max_repeats( core::Size const max_repeats_in );
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
	return utility::pointer::make_shared< SymDockProtocol >();
}

void SymDockProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymDockProtocol::provide_xml_schema( xsd );
}



} // symmetric_docking
} // protocols
