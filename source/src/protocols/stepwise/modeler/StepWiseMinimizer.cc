// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/stepwise/legacy/modeler/protein/util.hh> // for output_pose_list, maybe should deprecate soon
#include <protocols/farna/RNA_LoopCloser.hh>
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
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <ObjexxFCL/format.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.StepWiseMinimizer" );
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
	protein_ccd_closer_( protein::loop_close::StepWiseProteinCCD_CloserOP( new protein::loop_close::StepWiseProteinCCD_Closer( working_parameters ) ) )
{
	set_native_pose( working_parameters->working_native_pose() );
	runtime_assert( !options_->skip_coord_constraints() );
	runtime_assert( !options_->skip_minimize() );
	runtime_assert( !options_->move_jumps_between_chains() );
	// note this choice -- found that 5 really is necessary for protein stepwise monte carlo. can't minimize all due to slowdown.
	if ( pose_list_.size() > 0 && options_->choose_random() && num_pose_minimize_ == 0 /*asking this class for default*/ ) {
		num_pose_minimize_ =  ( protein::contains_protein( *pose_list[1] ) ? 5 : 1 );
	}
}

//Destructor
StepWiseMinimizer::~StepWiseMinimizer()
{}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::apply( core::pose::Pose & pose ) {
	if ( pose_list_.size() == 0 ) pose_list_ = make_vector1( pose.clone() );

	setup_minimizers();
	setup_scorefxns( pose );

	do_full_minimizing( pose );
	do_clustering( pose );

	// could probably deprecate soon, or at least use a function unified with RNA.
	if ( options_->output_minimized_pose_list() ) {
		stepwise::legacy::modeler::protein::output_pose_list( pose_list_, get_native_pose(),
			options_->silent_file(), working_calc_rms_res_ );
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::do_clustering( pose::Pose & pose ){
	using namespace protocols::stepwise::modeler::align;
	StepWiseClusterer clusterer( options_ );
	clusterer.set_pose_list( pose_list_ );
	clusterer.set_max_decoys( pose_list_.size() );
	clusterer.set_calc_rms_res( working_moving_res_ );

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
	if ( pose_list_.size() > 0 && num_pose_minimize_ != 1 ) TR << "Will minimize "  << pose_list_.size() << " poses." << std::endl;

	for ( Size n = 1; n <= pose_list_.size(); n++ ) {

		if ( num_pose_minimize_ > 0 &&  n > num_pose_minimize_ ) break;

		pose = *pose_list_[ n ];
		Real const score_original = pose.energies().total_energy();

		if ( options_->rm_virt_phosphate() ) rna::remove_all_virtual_phosphates( pose ); // ERRASER.

		// The movemap has all dofs for "non-fixed residues" free to move.
		get_move_map_and_allow_insert( mm, pose );

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
	//std::string const min_type = "dfpmin_armijo_nonmonotone";
	//Real const min_tolerance = 0.000025;
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
void
StepWiseMinimizer::do_minimize( pose::Pose & pose, kinematics::MoveMap & mm ){
	rna::o2prime_trials( pose, minimize_scorefxn_, working_pack_res_, allow_virtual_o2prime_hydrogens_ );
	if ( options_->vary_polar_hydrogen_geometry() ) polar_hydrogens::pack_polar_hydrogens( pose, allow_virtual_o2prime_hydrogens_ );

	core::scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	kinematics::MoveMap mm_save = mm;
	setup_vary_bond_geometry( pose, mm ); // careful -- must only do once, or constraints will keep getting added...

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

	farna::RNA_LoopCloser rna_loop_closer;
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
StepWiseMinimizer::get_move_map_and_allow_insert( core::kinematics::MoveMap & mm, pose::Pose const & pose ){
	utility::vector1< Size > working_minimize_res;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) { if ( !working_fixed_res_.has_value( n ) ) working_minimize_res.push_back( n );}
	bool const move_takeoff_torsions = !options_->disable_sampling_of_loop_takeoff();
	allow_insert_ = toolbox::AllowInsertOP( new toolbox::AllowInsert( pose ) ); // can come in handy later...
	movemap::figure_out_stepwise_movemap( mm, allow_insert_, pose, working_minimize_res, move_takeoff_torsions );
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
		for ( Size q = 1; q <= Hpos_polar.size(); q++ ) allow_insert_->set( id::AtomID( Hpos_polar[q], j ), true );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMinimizer::setup_vary_bond_geometry( core::pose::Pose & pose, core::kinematics::MoveMap & mm ) {
	if ( options_->vary_rna_bond_geometry() ) {
		TR << TR.Magenta << "Performing variable geometry minimization for RNA backbone..." << TR.Reset << std::endl;
		if ( !minimize_scorefxn_->has_nonzero_weight( core::scoring::bond_geometry ) ) minimize_scorefxn_->set_weight( core::scoring::bond_geometry, 1.0 );
		simple_moves::setup_vary_rna_bond_geometry( mm, pose, allow_insert_, core::scoring::bond_geometry );
	}
	if ( options_->vary_polar_hydrogen_geometry() ) {
		TR << TR.Magenta << "Performing variable geometry minimization for polar hydrogens..." << TR.Reset << std::endl;
		if ( !minimize_scorefxn_->has_nonzero_weight( core::scoring::bond_geometry ) ) minimize_scorefxn_->set_weight( core::scoring::bond_geometry, 1.0 );
		simple_moves::setup_vary_polar_hydrogen_geometry( mm, pose, allow_insert_ );
	}
}


} //modeler
} //stepwise
} //protocols
