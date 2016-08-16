// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/StepWiseMinimizer.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/StepWiseMinimizer.hh>
#include <protocols/stepwise/modeler/align/StepWiseClusterer.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/scoring_util.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/polar_hydrogens/util.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/legacy/modeler/protein/util.hh> // for output_pose_list, maybe should deprecate soon
#include <protocols/stepwise/monte_carlo/util.hh> // for output_to_silent_file, for RNA
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/magnesium/util.hh>
#include <protocols/simple_moves/ConstrainToIdealMover.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <ObjexxFCL/format.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.functions.hh>
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh>
#include <fstream>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.StepWiseMinimizer" );
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using utility::tools::make_vector1;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Single, simplified minimizer for stepwise assembly and monte carlo of proteins and RNA.
//
// TO DO: restore native_edensity_score_checker, if it is needed for ERRASER. Also skip_minimize? -- need to check with Fang.
// TO DO: put in CCD closing of protein loops, in the same manner as RNA loops.
// TO DO: figure out which side-chains are neighbors based on moving_partition, not moving_res!
//
// not carried over from StepWiseRNA_Minimizer/StepWiseProteinMinimizer:
//    coordinate_constraints, rescore_only_, linear_chainbreak ramping.
//    exensive pass-screen framework (could reconstitute with StepWiseScreeners, though)
//////////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {

//Constructor
StepWiseMinimizer::StepWiseMinimizer( utility::vector1< pose::PoseOP > const & pose_list,
	working_parameters::StepWiseWorkingParametersCOP working_parameters,
	options::StepWiseModelerOptionsCOP options,
	core::scoring::ScoreFunctionCOP scorefxn):
	pose_list_( pose_list ),
	options_( options ),
	scorefxn_( scorefxn ),
	num_pose_minimize_( options_->num_pose_minimize() ),
	working_moving_res_( get_all_working_moving_res( working_parameters ) ),
	working_fixed_res_( working_parameters->working_fixed_res() ),
	working_calc_rms_res_( working_parameters->working_calc_rms_res() ), // only for output -- may deprecate.
	allow_virtual_o2prime_hydrogens_( options->allow_virtual_o2prime_hydrogens() && !options->o2prime_legacy_mode() ),
	protein_ccd_closer_( protein::loop_close::StepWiseProteinCCD_CloserOP( new protein::loop_close::StepWiseProteinCCD_Closer( working_parameters ) ) ),
	working_parameters_( working_parameters ) // needed only for legacy SWA RNA main output.
{
	set_native_pose( working_parameters->working_native_pose() );
	runtime_assert( !options_->skip_coord_constraints() );
	runtime_assert( !options_->skip_minimize() );
	runtime_assert( !options_->move_jumps_between_chains() );
	if ( pose_list_.size() > 0 ) {
		pose::Pose const & pose = *pose_list[1];
		working_minimize_res_ = figure_out_working_minimize_res( pose );
		// note this choice -- found that 5 really is necessary for protein stepwise monte carlo. can't minimize all due to slowdown.
		if ( options_->choose_random() && num_pose_minimize_ == 0 /*asking this class for default*/ ) {
			num_pose_minimize_ =  ( protein::contains_protein( pose, working_minimize_res_ ) ? 5 : 1 );
		}
	} else {
		runtime_assert( options_->rna_legacy_output_mode() );
	}

}

//Destructor
StepWiseMinimizer::~StepWiseMinimizer()
{}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::apply( core::pose::Pose & pose ) {

	if ( !check_pose_list( pose ) ) return;

	setup_minimizers();
	setup_scorefxns( pose );

	do_full_minimizing( pose );
	do_clustering( pose );

	output_minimized_pose_list();

}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::do_clustering( pose::Pose & pose ){
	using namespace protocols::stepwise::modeler::align;
	StepWiseClusterer clusterer( options_ );
	clusterer.set_pose_list( pose_list_ );
	clusterer.set_max_decoys( pose_list_.size() );
	runtime_assert( working_minimize_res_ == figure_out_working_minimize_res( pose ) );
	clusterer.set_calc_rms_res( working_minimize_res_ );

	if ( pose_list_.size() > 1 && clusterer.cluster_rmsd() > 0.0 ) TR << "Will cluster "  << pose_list_.size() << " poses with cluster radius " << clusterer.cluster_rmsd() << std::endl;
	clusterer.cluster();

	pose_list_ = clusterer.pose_list();
	pose = *pose_list_[ 1 ];
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::do_full_minimizing( pose::Pose & pose ){

	kinematics::MoveMap mm;
	utility::vector1< pose::PoseOP > output_pose_list;
	if ( pose_list_.size() > 0 && num_pose_minimize_ > 1 ) TR << "Will minimize "  << std::min<Size>( pose_list_.size(), num_pose_minimize_ )  << " poses." << std::endl;

	for ( Size n = 1; n <= pose_list_.size(); n++ ) {

		if ( num_pose_minimize_ > 0 && n > num_pose_minimize_ ) break;

		pose = *pose_list_[ n ];
		Real const score_original = pose.energies().total_energy();

		if ( options_->rm_virt_phosphate() ) rna::remove_all_virtual_phosphates( pose ); // ERRASER.

		// The movemap has all dofs for "non-fixed residues" free to move.
		get_move_map_and_atom_level_domain_map( mm, pose );

		// We can also let sidechains minimize in fixed-residues -- for
		// speed only look at neighbors of moving residues.
		let_neighboring_side_chains_minimize( mm, pose );

		Real const score_before_min = (*minimize_scorefxn_)( pose );
		do_minimize( pose, mm );

		close_chainbreaks( pose, mm );

		TR << "Score minimized from " << F(8,3, score_before_min) << " to " << F(8,3,(*minimize_scorefxn_)( pose )) << "   [original: " << F(8,3,score_original);
		if ( hasPoseExtraScore( pose, "cluster_size" ) ) TR << " with cluster_size " << I( 4, getPoseExtraScore( pose, "cluster_size" ) );
		TR <<  "]" <<  std::endl;
		( *final_scorefxn_ )( pose );
		output_pose_list.push_back( pose.clone() );
	}

	pose_list_ = output_pose_list;
	working_pack_res_.clear();
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::setup_minimizers() {
	using namespace core::optimization;
	atom_tree_minimizer_ = core::optimization::AtomTreeMinimizerOP( new AtomTreeMinimizer );
	cartesian_minimizer_ = core::optimization::CartesianMinimizerOP( new CartesianMinimizer );

	bool const use_nblist( true );
	minimizer_options_ = core::optimization::MinimizerOptionsOP( new MinimizerOptions( options_->min_type() /*default */, options_->min_tolerance() /* default */, use_nblist, false, false ) );
	minimizer_options_->nblist_auto_update( true );
}

void
StepWiseMinimizer::setup_scorefxns( pose::Pose const & pose ) {
	using namespace core::scoring;
	minimize_scorefxn_ = get_minimize_scorefxn( pose, scorefxn_, options_ );

	// kind of ridiculous -- only here because rna_chem_map is so slow right now.
	final_scorefxn_ = minimize_scorefxn_->clone();
	final_scorefxn_->set_weight( rna_chem_map, scorefxn_->get_weight( rna_chem_map ) );
}

//////////////////////////////////////////////////////////////////////////
// quick hack -- testing water movement around Mg(2+)
void
freeze_waters( core::pose::Pose const & pose, core::kinematics::MoveMap & mm )
{
	utility::vector1< Size> const water_res = magnesium::get_water_res( pose );
	for ( Size n = 1; n <= water_res.size(); n++ ) {
		mm.set_jump( pose.fold_tree().get_jump_that_builds_residue( water_res[ n ] ), false );
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::do_minimize( pose::Pose & pose, kinematics::MoveMap & mm ){
	rna::o2prime_trials( pose, minimize_scorefxn_, working_pack_res_, allow_virtual_o2prime_hydrogens_ );
	if ( options_->vary_polar_hydrogen_geometry() ) polar_hydrogens::pack_polar_hydrogens( pose, allow_virtual_o2prime_hydrogens_ );
	if ( options_->hydrate_magnesiums() ) magnesium::hydrate_magnesiums( pose, true, options_->test_all_mg_hydration_frames() ); // this might be better as part of MasterPacker, in connection sampler.

	core::scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	kinematics::MoveMap mm_save = mm;
	setup_vary_bond_geometry( pose, mm ); // careful -- must only do once, or constraints will keep getting added...
	rna::add_syn_anti_chi_constraints( pose );
	if ( options_->cart_min() ) {
		cartesian_minimizer_->run( pose, mm, *minimize_scorefxn_, *minimizer_options_ );
	} else {
		atom_tree_minimizer_->run( pose, mm, *minimize_scorefxn_, *minimizer_options_ );
	}
	pose.constraint_set( save_pose_constraints );
	mm = mm_save;


	rna::o2prime_trials( pose, minimize_scorefxn_, working_pack_res_, allow_virtual_o2prime_hydrogens_ );

}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::close_chainbreaks( pose::Pose & pose, kinematics::MoveMap & mm ){

	utility::vector1< Size > const moving_chainbreaks = figure_out_moving_cutpoints_closed_from_moving_res( pose, working_moving_res_ ); // would be better to use moving_res.
	if ( moving_chainbreaks.size() == 0 ) return;

	farna::movers::RNA_LoopCloser rna_loop_closer;
	rna_loop_closer.apply( pose, rna::just_rna( moving_chainbreaks, pose ) ); // only handles RNA chainbreaks for now...

	utility::vector1< Size > const protein_cutpoints_closed = protein::just_protein( moving_chainbreaks, pose );
	for ( Size n = 1; n <= protein_cutpoints_closed.size(); n++ ) {
		protein_ccd_closer_->set_ccd_close_res( protein_cutpoints_closed[ n ] );
		protein_ccd_closer_->init( pose );
		protein_ccd_closer_->apply( pose );
	}

	do_minimize( pose, mm );
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::get_move_map_and_atom_level_domain_map( core::kinematics::MoveMap & mm, pose::Pose const & pose ){
	runtime_assert( working_minimize_res_ == figure_out_working_minimize_res( pose ) );
	bool const move_takeoff_torsions = !options_->disable_sampling_of_loop_takeoff();
	atom_level_domain_map_ = toolbox::AtomLevelDomainMapOP( new toolbox::AtomLevelDomainMap( pose ) ); // can come in handy later...
	movemap::figure_out_stepwise_movemap( mm, atom_level_domain_map_, pose, working_minimize_res_, move_takeoff_torsions );
	if ( !options_->minimize_waters() ) freeze_waters( pose, mm );
	output_movemap( mm, pose, TR.Debug );
}

//////////////////////////////////////////////////////////////////////////
// "side-chains" refers to proteins and 2'-OH for RNA
//
// A problem: what side-chains really should move? ones that form new contacts,
//  but also ones that *previously* were in contact before modeler -- they should relax too.
//
void
StepWiseMinimizer::let_neighboring_side_chains_minimize( core::kinematics::MoveMap & mm,
	core::pose::Pose const & pose ) {
	using namespace core::scoring;
	if ( options_->global_optimize() ) working_pack_res_ = get_all_residues( pose );
	if ( working_pack_res_.size() == 0 ) {
		working_pack_res_ = packer::figure_out_working_interface_res( pose, working_moving_res_ );
	}
	for ( Size n = 1; n <= working_pack_res_.size(); n++ ) {
		move_side_chain( mm, pose, working_pack_res_[n] );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::move_side_chain( core::kinematics::MoveMap & mm,
	pose::Pose const & pose,
	Size const j ){
	if ( mm.get_chi( j ) ) return;
	if ( pose.residue(j).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) return;
	if ( pose.residue(j).is_protein() ) {
		mm.set_chi( j, true );
	} else if ( pose.residue(j).is_RNA() ) {
		mm.set( id::TorsionID( j, id::CHI, 4), true ); // 2'-OH.
	}
	// how about terminal phosphates?
	if ( options_->vary_polar_hydrogen_geometry() ) {
		utility::vector1< Size > const & Hpos_polar = pose.residue( j ).Hpos_polar();
		for ( Size q = 1; q <= Hpos_polar.size(); q++ ) atom_level_domain_map_->set( id::AtomID( Hpos_polar[q], j ), true );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::setup_vary_bond_geometry( core::pose::Pose & pose, core::kinematics::MoveMap & mm ) {
	if ( options_->vary_rna_bond_geometry() ) {
		TR << TR.Magenta << "Performing variable geometry minimization for RNA backbone..." << TR.Reset << std::endl;
		if ( !minimize_scorefxn_->has_nonzero_weight( core::scoring::bond_geometry ) ) minimize_scorefxn_->set_weight( core::scoring::bond_geometry, 1.0 );
		simple_moves::setup_vary_rna_bond_geometry( mm, pose, atom_level_domain_map_, core::scoring::bond_geometry );
	}
	if ( options_->vary_polar_hydrogen_geometry() ) {
		TR << TR.Magenta << "Performing variable geometry minimization for polar hydrogens..." << TR.Reset << std::endl;
		if ( !minimize_scorefxn_->has_nonzero_weight( core::scoring::bond_geometry ) ) minimize_scorefxn_->set_weight( core::scoring::bond_geometry, 1.0 );
		simple_moves::setup_vary_polar_hydrogen_geometry( mm, pose, atom_level_domain_map_ );
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::output_empty_minimizer_silent_file() const {

	if ( utility::file::file_exists( options_->silent_file() ) ) utility_exit_with_message( options_->silent_file() + " already exist!" );

	std::ofstream outfile;
	outfile.open( options_->silent_file().c_str() ); //Opening the file with this command removes all prior content..

	outfile << "StepWiseMinimizer:: num_pose_outputted == 0, empty silent_file!\n"; //specific key signal to SWA_modeler_post_process.py

	outfile.flush();
	outfile.close();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////                                                                         void
bool
StepWiseMinimizer::check_pose_list( core::pose::Pose const & pose ){
	if ( pose_list_.size() == 0 ) {
		// Required for legacy mode
		if ( options_->rna_legacy_output_mode() ) {
			TR.Debug << "pose_list_.size() == 0, early exit from StepWiseMinimizer::apply" << std::endl;
			output_empty_minimizer_silent_file();
			return false;
		} else {
			pose_list_ = make_vector1( pose.clone() );
			return true;
		}
	}
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////                                                                         void
// the minimizer really should not be outputting anything -- this is historical from Parin's
// stepwise assembly code from 2009-2010.
// could probably deprecate soon along with swa_rna_main and swa_protein_main.
void
StepWiseMinimizer::output_minimized_pose_list() const{
	if ( options_->output_minimized_pose_list() ) {
		if ( options_->rna_legacy_output_mode() ) {
			// copied from StepWiseRNA_VirtualSugarSamplerFromStringList -- trying to be consistent with that crazy old thing.
			//    output_data( silent_file_out_, pose_tag, false, pose, working_parameters_->working_native_pose(), working_parameters_ );
			for ( Size n = 1; n <= pose_list_.size(); n++ ) {
				Pose & pose = ( *pose_list_[n] ); //set viewer_pose;
				std::string const tag = "S_"+ ObjexxFCL::string_of( n-1 /* start with zero */);
				stepwise::modeler::rna::output_data( options_->silent_file(),
					tag,
					false /* write score only*/,
					pose,
					get_native_pose(),
					working_parameters_,
					true /*NAT_rmsd*/  );
			}
		} else if ( !protein::contains_protein( *pose_list_[1] ) ) {
			for ( Size n = 1; n <= pose_list_.size(); n++ ) {
				Pose & pose = ( *pose_list_[n] ); //set viewer_pose;
				std::string const tag = "S_"+ ObjexxFCL::string_of( n-1 /* start with zero */);
				stepwise::monte_carlo::output_to_silent_file( tag,
					options_->silent_file(),
					pose,
					get_native_pose(),
					false /*superimpose_over_all*/,
					true /*rms_fill*/ );
			}
		} else { // well its protein legacy output mode then...
			stepwise::legacy::modeler::protein::output_pose_list( pose_list_,
				get_native_pose(),
				options_->silent_file(),
				working_calc_rms_res_  );

		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::Size >
StepWiseMinimizer::figure_out_working_minimize_res( core::pose::Pose const & pose ) {
	using namespace core::pose::full_model_info;
	utility::vector1< core::Size > working_minimize_res;
	FullModelInfo const & full_model_info( const_full_model_info( pose ) );
	utility::vector1< core::Size > const working_extra_minimize_res( full_model_info.full_to_sub( full_model_info.extra_minimize_res() ) );

	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( pose.residue_type( n ).name3() == "HOH" ) continue; // allowing waters to vary in different poses (Mg hydration)
		if ( !working_fixed_res_.has_value( n ) ||
				working_extra_minimize_res.has_value( n ) ) {
			working_minimize_res.push_back( n );
		}
	}
	return working_minimize_res;
}



} //modeler
} //stepwise
} //protocols
