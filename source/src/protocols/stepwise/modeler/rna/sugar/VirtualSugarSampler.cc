// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file VirtualSugarSampler
/// @brief
/// @details
/// @author Parin Sripakdeevong, Rhiju Das (rhiju@stanford.edu)


//////////////////////////////////
#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarSampler.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/sampler/rna/RNA_SuiteStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_NucleosideStepWiseSampler.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/full_model_info/util.hh> // needed for chain/gap-size.
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

using namespace core;
using namespace core::pose;
using utility::tools::make_vector1;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.sugar.VirtualSugarSampler" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Parin's SWA code would leave residues with the following virtualization status after
//   'floating base' ( a.k.a., dinucleotide, or "skip bulge") moves:
//
//
//       Anchor residue    Bulge residue         Moving residue             partition that needs to be virtualized.
//                      (virtual) (virtual)     (virtual)   (virtual)  ------------------------->
//     5'   --Sugar --  Phosphate -- Sugar  -x- Phosphate -- Sugar --  Phosphate -- Sugar --
//              |                     |      |                |                      |
//          Reference               Bulge    |             Floating                Distal
//             Base                (virtual) |                  Base!               Base
//              |                          cutpoint            |
//              |______________________________________________|
//                     Jump between anchor & moving
//
//  This class instantiates and samples the virtual ribose, including some minimization of the ribose &
//   moving base.
//
//  There is also a virtualization of the 'other partition' - that use case involves
//   merging two poses -- the moving_suite connects the poses, and this VirtualSugarSampler is called
//   to deal with virtual sugars right before or after that suite.
//
// A lots of this stuff is deprecated...
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

//Constructor
VirtualSugarSampler::VirtualSugarSampler( working_parameters::StepWiseWorkingParametersCOP & working_parameters,
	SugarModeling & sugar_modeling ) :
	working_parameters_( working_parameters ),
	sugar_modeling_( sugar_modeling ), // holds both job parameters and pose list.
	tag_( "" ),
	use_phenix_geo_( false ),
	keep_base_fixed_( false ),
	choose_random_( false ),
	do_minimize_( true ),
	do_screens_( true ),
	integration_test_mode_( false ),
	virtual_sugar_is_from_prior_step_( true ),
	legacy_mode_( false ), //for comparison to original code -- may deprecate in 2014.
	do_chain_closure_( true || legacy_mode_ ),
	first_minimize_with_fixed_base_( false ),
	max_tries_for_random_sugar_setup_( 1 ),
	sugar_setup_success_( false ),
	moving_phosphate_virtualized_( false )
{}

//Destructor
VirtualSugarSampler::~VirtualSugarSampler()
{}

/////////////////////
std::string
VirtualSugarSampler::get_name() const {
	return "VirtualSugarSampler";
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::apply( pose::Pose & viewer_pose ){

	utility::vector1< PoseOP > pose_list;

	Size const max_tries = choose_random_ ? max_tries_for_random_sugar_setup_ : 1;
	for ( Size ntries = 1; ntries <= max_tries; ntries++ ) {
		TR.Debug << TR.Green << "Will sample residue " << sugar_modeling_.moving_res << " ( reference res: " << sugar_modeling_.reference_res <<  " ) " << TR.Reset << std::endl;
		apply( pose_list, viewer_pose );
		if ( pose_list.size() > 0 || !sugar_setup_success_ ) break;
	}
	sugar_modeling_.pose_list = pose_list; // stores the solutions.

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::apply( utility::vector1< PoseOP > & pose_list,
	pose::Pose & viewer_pose ){
	using namespace ObjexxFCL;
	using namespace core::io::silent;

	output_title_text( "Enter sample_virtual_sugar_and_bulge_and_close_chain()", TR.Debug );
	clock_t const time_start( clock() );
	pose::Pose const viewer_pose_save = viewer_pose;

	tag_into_pose( viewer_pose, tag_ );
	sugar_modeling_.check_compatibility( viewer_pose.total_residue() );
	original_constraint_set_ = viewer_pose.constraint_set()->clone();

	//Virtualize any residues that move 'with' the moving residue -- their locations have not been sampled yet.
	if ( do_screens_ ) {
		if ( virtual_sugar_is_from_prior_step_ ) virtualize_distal_partition( viewer_pose );
		else runtime_assert( working_parameters_->gap_size() );
	}

	// if ( !choose_random_ ) setup_VDW_bin_checker( viewer_pose ); // this takes too long.

	setup_sugar_conformations( pose_list, viewer_pose );

	if ( ( sugar_modeling_.bulge_res > 0 ) && do_chain_closure_ ) do_chain_closure_modeler( pose_list, viewer_pose );

	reinstate_original_constraints( pose_list );
	if ( do_screens_ && virtual_sugar_is_from_prior_step_ ) reinstantiate_distal_partition( pose_list );

	viewer_pose = viewer_pose_save;

	TR.Debug << "Time in sample_virtual_sugar_and_bulge_and_close_chain(): " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Instantiates the virtual ribose and goes through a few pucker and chi conformations
//
// Then do a minimize of residue without instantiated sugar, including base (!).
//
// Note also further pose variant changes:
//  Instantiate the bulged ribose that connect the virtual ribose to the anchor residue.
//  Virtualize phosphate 5' to any chainbreak. [I thought this would already be done...]
//
void
VirtualSugarSampler::setup_sugar_conformations( utility::vector1< PoseOP > & pose_list, pose::Pose & pose ) {

	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	using namespace core::id;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring::rna;
	using namespace sampler::rna;
	using namespace protocols::farna::movers;
	using namespace core::pose;

	clock_t const time_start( clock() );
	//July 28th, 2011 Could set extra_chi here, BUT necessary?
	RNA_NucleosideStepWiseSampler sampler( sugar_modeling_.moving_res,
		sugar_modeling_.moving_res_pucker_state,
		sugar_modeling_.moving_res_base_state );
	sampler.set_idealize_coord(   use_phenix_geo_ );
	sampler.set_skip_same_pucker( use_phenix_geo_ );
	sampler.set_random( choose_random_ );
	sampler.init();

	checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker;
	checker::RNA_AtrRepCheckerOP atr_rep_checker;
	Size const five_prime_gap_res  = sugar_modeling_.is_prepend ? sugar_modeling_.moving_res    : sugar_modeling_.reference_res;
	Size const three_prime_gap_res = sugar_modeling_.is_prepend ? sugar_modeling_.reference_res : sugar_modeling_.moving_res;

	if ( do_screens_ ) {
		runtime_assert( sugar_modeling_.reference_res > 0 );

		// Deprecate soon:
		//  Size const moving_suite = sugar_modeling_.is_prepend ? sugar_modeling_.moving_res : ( sugar_modeling_.moving_res - 1);
		//TR << TR.Red << "MOVING " << sugar_modeling_.moving_res << "  BULGE " << sugar_modeling_.bulge_res << "  REFERENCE " << sugar_modeling_.reference_res << " IS_PREPEND " << sugar_modeling_.is_prepend << " MOVING_SUITE " << moving_suite << " sugar_modeling:FIVE_PRIME_CHAINBREAK " << sugar_modeling_.five_prime_chain_break << " FIVE_PRIME_GAP_RES " << five_prime_gap_res << " DO_MINIMIZE " << do_minimize_ << TR.Reset << std::endl;
		using namespace core::pose::full_model_info;
		utility::vector1< Size > const chains = figure_out_chain_numbers_from_full_model_info( pose );
		if ( chains[ five_prime_gap_res ] == chains[ three_prime_gap_res ] ) {
			utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
			Size const gap_size = res_list[ three_prime_gap_res ] -  res_list[ five_prime_gap_res ];
			checker::RNA_ChainClosableGeometryChecker chain_closable_geometry_checker( five_prime_gap_res, three_prime_gap_res, gap_size );
		}

		utility::vector1< Size > const & moving_partition_res = working_parameters_->working_moving_partition_res();
		if ( moving_partition_res.has_value( sugar_modeling_.moving_res ) == moving_partition_res.has_value( sugar_modeling_.reference_res ) ) {
			atr_rep_checker = checker::RNA_AtrRepCheckerOP( new checker::RNA_AtrRepChecker( pose, sugar_modeling_.moving_res, sugar_modeling_.reference_res, 0 /*gap_size = 0 forces loose rep_cutoff (10.0 )*/, scorefxn_->energy_method_options().clone() ) );
		}
	}

	pose::Pose pose_init = pose;//Hard copy
	pose::remove_variant_type_from_pose_residue( pose_init,
		core::chemical::VIRTUAL_RIBOSE, sugar_modeling_.moving_res );

	Size count( 0 );
	sugar_setup_success_ = false;
	pose_list.clear();
	for ( sampler.reset(); sampler.not_end(); ++sampler ) {
		count++;
		if ( choose_random_ && count > max_tries_for_random_sugar_setup_ ) {
			pose_list.push_back( pose.clone() ); // this routine must return something with an instantiated sugar, even if its a fail.
			break;
		}

		pose = pose_init;
		// following *must* be mistake... should occur below inside actual chain_closure..
		if ( legacy_mode_ ) {
			add_harmonic_chain_break_constraint( pose, sugar_modeling_.five_prime_chain_break );
			if  ( do_chain_closure_ ) setup_chain_break_variants( pose, sugar_modeling_.five_prime_chain_break );
		} else if ( do_screens_ ) {
			add_fade_chain_break_constraint_across_gap( pose, five_prime_gap_res, three_prime_gap_res, 1 /*gap_size -- one intervening bulge residue -- perhaps should not be hard-wired! */ );
		}

		sampler.apply( pose );

		if ( do_minimize_ ) minimize_sugar( pose ); // remove clashes; show in graphics too.
		if ( ( chain_closable_geometry_checker != 0 ) && !chain_closable_geometry_checker->check_screen( pose ) ) continue;
		if ( ( atr_rep_checker != 0 ) &&  !atr_rep_checker->check_screen( pose ) ) continue;

		tag_into_pose( pose, tag_from_pose( pose_init ) + '_' + ObjexxFCL::string_of( count ) );
		pose_list.push_back( pose.clone() );
		sugar_setup_success_ = true;
		if ( choose_random_ ) break; // take the first random sample.
	}

	output_title_text( "", TR.Debug );
	TR.Debug << " pose_list.size() = " << pose_list.size() << std::endl;
	TR.Debug << "Total time in bulge_chain_closure SETUP: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	pose = pose_init;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::minimize_sugar( pose::Pose & pose_with_sugar ){

	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::id;
	using namespace core::kinematics;

	ScoreFunctionOP sugar_scorefxn, sugar_scorefxn_without_ch_bond, rescaled_sugar_score_fxn_without_ch_bond;
	get_sugar_setup_scorefxns( sugar_scorefxn, sugar_scorefxn_without_ch_bond, rescaled_sugar_score_fxn_without_ch_bond );

	AtomTreeMinimizer minimizer;
	bool const use_nblist( true );
	float const tolerance = 0.000000000000025; //Sept 21, 2011, same result for  dfp_min and dfpmin_atol| converge to identical energy score with 0.00000025 of dfpmin_atol!
	//Fang: This number is insanely small....
	MinimizerOptions options_standard( "lbfgs_armijo_atol", tolerance, use_nblist, false, false );      //Switch to absolute tolerance on Sept 21, 2011
	MinimizerOptions options_armijo( "lbfgs_armijo_atol", tolerance, use_nblist, false, false ); //Add this on Sept 21, 2011

	TR.Debug << "options_standard: min_type = " << options_standard.min_type() << " minimize_tolerance = "  << options_standard.minimize_tolerance() << std::endl;
	TR.Debug << "options_armijo  : min_type = " << options_armijo.min_type()   << " minimize_tolerance = "  << options_armijo.minimize_tolerance()   << std::endl;

	options_standard.nblist_auto_update( true );
	options_armijo.nblist_auto_update( true );

	MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	if ( first_minimize_with_fixed_base_ || keep_base_fixed_ ) {
		mm.set( TorsionID( sugar_modeling_.moving_res, id::CHI, 1 ), true ); // nucleosidic chi
		mm.set( TorsionID( sugar_modeling_.moving_res, id::CHI, 2 ), true ); // sugar torsion
		mm.set( TorsionID( sugar_modeling_.moving_res, id::CHI, 3 ), true ); // sugar torsion
		minimizer.run( pose_with_sugar, mm, *( sugar_scorefxn_without_ch_bond ), options_standard );
	}

	if ( keep_base_fixed_ ) return;

	bool found_desired_jump_ID = false;
	for ( Size jump_ID = 1; jump_ID <= pose_with_sugar.fold_tree().num_jump(); jump_ID++ ) {
		Size const jump_pos1( pose_with_sugar.fold_tree().upstream_jump_residue( jump_ID ) );
		Size const jump_pos2( pose_with_sugar.fold_tree().downstream_jump_residue( jump_ID ) );
		if ( ( jump_pos1 == sugar_modeling_.moving_res ) || ( jump_pos2 == sugar_modeling_.moving_res ) ) {
			found_desired_jump_ID = true;
			mm.set_jump( jump_ID, true );
		}
	}
	runtime_assert( found_desired_jump_ID );

	/////////////Switch to armijo on Sept 21, 2011///////////////////////////////////////////////////////////////////////////////////////
	/////////////My understanding is that dfpmin_armijo is a "inexact" line search whereas the standard dfpmin is a exact line search///////
	/////////////It seem to indicate the dfpmin should be slower (require more function evaluation) but at the same time more accurate//////
	/////////////See http://www.rosettacommons.org/manuals/archive/rosetta3.3_user_guide/minimization_overview.html for details/////////////
	/////////////However standard dfpmin seem to lead cases where the floating base just "explode" and more far away from starting point////
	/////////////This sometimes lead to the Hbond tripped error/////////////////////////////////////////////////////////////////////////////
	/////////////Side note: switching to dfpmin_atol (atol-> absolute tolerance didn't help!)///////////////////////////////////////////////
	/////////////So switching to dfpmin_armijo which doesn't seem to exhibit this behavior//////////////////////////////////////////////////
	/////////////Note that there is a currently a bug in in dfpmin_armijo:
	/////////////core.optimization.LineMinimizer: Inaccurate G! step= 9.53674e-07 Deriv= -0.0226443 Finite Diff= 0.00628252/////////////////
	/////////////Rhiju mention that this bug is fixed in the latest Rosetta version in trunk////////////////////////////////////////////////
	/////////////So will keep using standard dfp_min except at the first minimiziation step/////////////////////////////////////////////////
	/////////////Also tried two round minimizations with the first using options_armijo. This fixes the "explode" bug but led to worst score!/
	minimize_with_constraints( pose_with_sugar, mm, rescaled_sugar_score_fxn_without_ch_bond, options_armijo ); //Add this round on Sept 20, 2011, Switch to armijo on Sept 21, 2011
	minimizer.run( pose_with_sugar, mm, *( rescaled_sugar_score_fxn_without_ch_bond ), options_armijo );        //Add this round on Sept 20, 2011, Switch to armijo on Sept 21, 2011
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	minimize_with_constraints( pose_with_sugar, mm, sugar_scorefxn_without_ch_bond, options_standard );
	minimizer.run( pose_with_sugar, mm, *( sugar_scorefxn_without_ch_bond ), options_standard );

	minimize_with_constraints( pose_with_sugar, mm, sugar_scorefxn, options_standard );
	minimizer.run( pose_with_sugar, mm, *( sugar_scorefxn ), options_standard );

}


/////////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::get_sugar_setup_scorefxns( scoring::ScoreFunctionOP & sugar_scorefxn, scoring::ScoreFunctionOP & sugar_scorefxn_without_ch_bond, scoring::ScoreFunctionOP & rescaled_sugar_score_fxn_without_ch_bond ) const {

	using namespace core::scoring;

	sugar_scorefxn = scorefxn_->clone();
	sugar_scorefxn->set_weight( linear_chainbreak, 0.0 );
	sugar_scorefxn->set_weight( angle_constraint, 0.0 );
	sugar_scorefxn->set_weight( atom_pair_constraint, 0.0 );
	sugar_scorefxn->set_weight( coordinate_constraint, 0.1 );

	sugar_scorefxn_without_ch_bond = sugar_scorefxn->clone();
	sugar_scorefxn_without_ch_bond->set_weight( ch_bond, 0.0 );
	//This makes sure that there are no chain_break score involved.

	//Sept 20, 2011 To solve the problem with the floating base res exploding when minimizing
	//Problem case occur when floating base res is the 1st working res
	//Note that this error doesn't seem to occur if virtual_sugar is sampled in same step as SAMPLER/ (no Hbond_tripped!)
	//The non-rescale scorefxn does however causes the floating base from moving far away from the starting
	//point even in the same step as SAMPLER case
	//Also generally the minimizer same and separate step virtual sampler doesn't give the same results!
	TR.Debug << "--------------START Creating rescaled one_tenth_sugar_score_fxn_without_ch_bond--------------" << std::endl;
	rescaled_sugar_score_fxn_without_ch_bond = rescale_scorefxn( sugar_scorefxn_without_ch_bond, 0.1 );
	TR.Debug << "--------------FINISH Creating rescaled one_tenth_sugar_score_fxn_without_ch_bond--------------" << std::endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////
}


/////////////////////////////////////////////////////////////////////////////////////
// chain closure -- sample the bulge suite, i.e. the suite connecting the anchor to the bulge.
// Then close at the other suite, connecting the bulge to the 'moving residue' [which just had its virtual sugar instantiated]
//
//                   SAMPLE                  CCD-CLOSE
//                     |                        |
//       Anchor residue    Bulge residue        Moving residue
//     5'   --Sugar --  Phosphate -- Sugar  -- Phosphate -- Sugar --
//              |                     |                       |
//          Reference               Bulge                  Floating
//             Base                (virtual)                 Base!
//
// Could probably also use analytical loop closure here.
//
void
VirtualSugarSampler::do_chain_closure_modeler( utility::vector1< PoseOP > & pose_list,
	pose::Pose & viewer_pose ){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring::rna;
	using namespace ObjexxFCL;

	if ( pose_list.size() == 0 ) return; //early return

	initialize_pose_variants_for_chain_closure( pose_list );
	bulge_chain_closure( pose_list, viewer_pose );
	restore_pose_variants_after_chain_closure( pose_list );

}

/////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::bulge_chain_closure( utility::vector1< PoseOP > & pose_list,
	pose::Pose & viewer_pose ){
	if ( legacy_mode_ ) { // deprecate legacy code soon
		bulge_chain_closure_legacy( pose_list, viewer_pose );
		bulge_chain_minimize_legacy( pose_list, viewer_pose );
	} else {
		bulge_chain_closure_complete( pose_list, viewer_pose );
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
// through the use of KIC-based modeler, this bulge modeling is both more comprehensive and
// faster than the legacy code [which was CCD-based].
void
VirtualSugarSampler::bulge_chain_closure_complete( utility::vector1< PoseOP > & pose_list,
	pose::Pose & viewer_pose ){
	utility::vector1< PoseOP > output_pose_list;
	Size const rebuild_res = sugar_modeling_.bulge_res;
	StepWiseModelerOP stepwise_modeler( new StepWiseModeler( rebuild_res, scorefxn_ ) );
	options::StepWiseModelerOptionsOP options( new options::StepWiseModelerOptions );
	options->set_use_phenix_geo( use_phenix_geo_ );
	options->set_kic_modeler_if_relevant( true /*erraser modeler for speed & completeness*/ );
	options->set_num_pose_minimize( 1 );
	options->set_integration_test_mode( true /*fast modeler*/ );
	options->set_disallow_realign( true );
	stepwise_modeler->set_options( options );
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		pose::Pose & pose( *pose_list[n] );
		viewer_pose = pose;
		stepwise_modeler->set_moving_res_and_reset( rebuild_res );
		stepwise_modeler->set_working_minimize_res( make_vector1( rebuild_res ) );
		stepwise_modeler->apply( viewer_pose );
		if ( stepwise_modeler->get_num_sampled() > 0 ) output_pose_list.push_back( viewer_pose.clone() );
		if ( integration_test_mode_ && output_pose_list.size() > 0 ) break;
	}
	pose_list = output_pose_list;
}

/////////////////////////////////////////////////////////////////////////////////////////
// deprecate this after march 2014 if new 'bulge_chain_closure_complete' is
// satisfactory for stepwise assembly & stepwise monte carlo runs.
void
VirtualSugarSampler::bulge_chain_closure_legacy( utility::vector1< PoseOP > & pose_list,
	pose::Pose & viewer_pose ){

	using namespace sampler::rna;
	using namespace protocols::stepwise::modeler::rna::checker;
	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	using namespace protocols::farna::movers;

	pose::Pose screening_pose = *(pose_list[1]);

	output_title_text( "bulge_chain_closure SAMPLING", TR.Debug );
	clock_t const time_start_modeler( clock() );

	// in legacy mode chainbreak constraint addition occurs above in 'setup'. but it really should be here.
	if ( !legacy_mode_ ) {
		for ( Size n = 1; n <= pose_list.size(); n++ ) {
			Pose & pose = *(pose_list[n]);
			add_harmonic_chain_break_constraint( pose, sugar_modeling_.five_prime_chain_break );
		}
	}

	Size const bulge_suite = sugar_modeling_.bulge_suite;
	Size const bulge_rsd   = sugar_modeling_.bulge_res;
	utility::vector1<bool> sample_sugar( 2, false );
	utility::vector1<ChiState> base_state  ( 2, ANY_CHI );
	utility::vector1<PuckerState> pucker_state( 2, ANY_PUCKER );

	for ( Size i = 1; i <= 2; ++i ) {
		Size const curr_rsd( bulge_suite + i - 1 );
		if ( bulge_rsd == curr_rsd ) {
			sample_sugar[i] = true;
			pucker_state[i] = ANY_PUCKER;
			base_state[i] = ANY_CHI;
		} else {
			sample_sugar[i] = false;
			//FANG: Assuming all pose in the list have same pucker for the unmoved
			//residue. Not sure if this is correct but it only affects
			//epsilon modeler.
			pucker_state[i] = assign_pucker( screening_pose, curr_rsd ); // NOT TRUE! FIX!
			base_state[i]   = NO_CHI;
		}
	}

	RNA_SuiteStepWiseSampler sampler( bulge_suite, pucker_state[1], pucker_state[2], base_state[1], base_state[2] );
	sampler.set_skip_same_pucker( use_phenix_geo_ );
	sampler.set_idealize_coord( use_phenix_geo_ );
	sampler.set_fast( integration_test_mode_ );
	sampler.set_extra_epsilon( true ); //Parin's original option. Necessary?
	sampler.set_sample_nucleoside_lower( sample_sugar[1] );
	sampler.set_sample_nucleoside_upper( sample_sugar[2] );
	sampler.init();

	StepWiseRNA_CountStruct count_data;
	RNA_LoopCloser rna_loop_closer;
	Size num_closed_chain_pose = 0;

	utility::vector1< RNA_AtrRepCheckerOP > atr_rep_checkers_;
	for ( Size n = 1; n <= pose_list.size(); n++ ) atr_rep_checkers_.push_back( rna::checker::RNA_AtrRepCheckerOP( new RNA_AtrRepChecker( *(pose_list[n]), bulge_suite, bulge_rsd, 0 /*gap_size*/ ) ) );

	RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker_( new RNA_ChainClosableGeometryChecker( sugar_modeling_.five_prime_chain_break, 0 /*gap_size*/ ) );
	RNA_ChainClosureCheckerOP chain_closure_checker_( new RNA_ChainClosureChecker( screening_pose, sugar_modeling_.five_prime_chain_break ) );
	chain_closure_checker_->set_reinitialize_CCD_torsions( true );

	utility::vector1< bool > is_chain_close( pose_list.size(), false );
	utility::vector1< pose::PoseOP > output_pose_list;

	// following is not a rigorous enumeration and kind of slow.
	// first, the sampler is not subdivided by pucker of the actual pose.
	// And... each backbone torsion is checked once, looking for any pose, rather
	// than trying all backbone against all poses.
	for ( sampler.reset(); sampler.not_end(); ++sampler ) {
		if ( integration_test_mode_ && num_closed_chain_pose > 1 ) break;

		if ( num_closed_chain_pose == pose_list.size() ) break;

		count_data.tot_rotamer_count++;
		sampler.apply( screening_pose );

		if ( VDW_bin_checker_ && !VDW_bin_checker_->VDW_rep_screen( screening_pose, sugar_modeling_.bulge_res ) ) continue;
		count_data.good_bin_rep_count++;

		for ( Size n = 1; n <= pose_list.size(); n++ ) {

			if ( is_chain_close[ n ] ) continue;
			Pose & current_pose = *(pose_list[n]);

			// don't copy torsions into current pose just quite yet -- check if chain is closable.
			if  ( !chain_closable_geometry_checker_->check_screen( screening_pose   /* moving_pose */,
					current_pose     /* reference_pose */,
					!sugar_modeling_.is_prepend,
					true /*strict*/ ) ) continue;
			count_data.chain_closable_geometry_count++;

			sampler.apply( current_pose );

			if ( !fast_full_atom_VDW_repulsion_screen( current_pose,
					sugar_modeling_.bulge_res, sugar_modeling_.moving_res,
					sugar_modeling_.is_prepend ) ) continue;
			count_data.fast_full_atom_VDW_repulsion_screen++;

			if ( !atr_rep_checkers_[n]->check_screen( current_pose ) ) continue;
			if ( chain_closure_checker_->check_screen( current_pose ) ) continue;

			is_chain_close[ n ] = true;
			num_closed_chain_pose++;
			viewer_pose = current_pose; // show in graphics.
			output_pose_list.push_back( pose_list[n] );
			if ( integration_test_mode_ ) break;
		}
	}

	TR.Debug << " bin_rep_count = " << count_data.good_bin_rep_count;
	TR.Debug << " fast_rep_count = " << count_data.fast_full_atom_VDW_repulsion_screen;
	TR.Debug << " chain_closable_geometry = " << count_data.chain_closable_geometry_count;
	TR.Debug << " angle_n = " << count_data.good_angle_count << " dist_n = " << count_data.good_distance_count;
	TR.Debug << " rep = " << count_data.good_rep_rotamer_count;
	TR.Debug << " rmsd = " << count_data.rmsd_count << " tot = " << count_data.tot_rotamer_count << std::endl;
	TR.Debug << " " << num_closed_chain_pose << " out of " << pose_list.size() << " pose were closable" << std::endl;
	TR.Debug << "Total time in bulge_chain_closure SAMPLING: " << static_cast< Real > ( clock() - time_start_modeler ) / CLOCKS_PER_SEC << std::endl;

	pose_list = output_pose_list;
}

/////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::bulge_chain_minimize_legacy( utility::vector1< PoseOP > & pose_list,
	pose::Pose & viewer_pose ){

	using namespace core::optimization;
	using namespace core::id;
	using namespace core::scoring;

	////////////////////////////////////////////////////////////////
	// post-process! -- why not just use existing Minimizer class?
	///////////////////////////////////////////////////////////////
	output_title_text( "bulge_chain_closure POST_PROCESSING", TR.Debug );

	//Quick minimize to reduce error in CCD///////////////////////////
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025 );
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	core::kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );

	//The whole richardson's suite ... avoid delta, nu2, nu1, WARNING VIRTUAL_PHOSPHATE is ON at this point...
	mm.set( TorsionID( sugar_modeling_.bulge_res - 1, id::BB,  5 ), true ); //epsilon
	mm.set( TorsionID( sugar_modeling_.bulge_res - 1, id::BB,  6 ), true ); //zeta

	mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  1 ), true ); //alpha
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  2 ), true ); //beta
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  3 ), true ); //gamma
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  5 ), true ); //epsilon
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  6 ), true ); //zeta
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::CHI, 1 ), true ); //chi (torsion between base and sugar sugar)

	mm.set( TorsionID( sugar_modeling_.bulge_res + 1, id::BB,  1 ), true ); //alpha
	mm.set( TorsionID( sugar_modeling_.bulge_res + 1, id::BB,  2 ), true ); //beta
	mm.set( TorsionID( sugar_modeling_.bulge_res + 1, id::BB,  3 ), true ); //gamma

	ScoreFunctionOP bulge_chain_closure_scorefxn( new ScoreFunction );
	bulge_chain_closure_scorefxn->set_weight( fa_rep , 0.12 );
	bulge_chain_closure_scorefxn->set_weight( angle_constraint, 1.0 );
	bulge_chain_closure_scorefxn->set_weight( atom_pair_constraint, 1.0 );
	bulge_chain_closure_scorefxn->set_weight( linear_chainbreak, 5.0 );
	bulge_chain_closure_scorefxn->set_weight( rna_sugar_close, 70 ); //100X normal weight

	mm.set( TorsionID( sugar_modeling_.bulge_res, id::CHI, 2 ), true ); //nu2
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::CHI, 3 ), true ); //nu1

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< PoseOP > output_pose_list;
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		viewer_pose = *(pose_list[n]);
		TR.Debug << "POST_PROCESSING pose # " << n << " out of " << pose_list.size() << " " << std::endl;;
		minimizer.run( viewer_pose, mm, *( bulge_chain_closure_scorefxn ), options );
		output_pose_list.push_back( viewer_pose.clone() );
	}
	pose_list = output_pose_list;
}


/////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::initialize_pose_variants_for_chain_closure( utility::vector1< pose::PoseOP > & pose_list ){
	// instantiate bulge.
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		Pose & pose = *(pose_list[n]);
		// in legacy mode this occurs above in 'setup'. but it really should be here.
		if ( !legacy_mode_ ) setup_chain_break_variants( pose, sugar_modeling_.five_prime_chain_break );
		pose::add_variant_type_to_pose_residue(  pose,
			core::chemical::VIRTUAL_PHOSPHATE, sugar_modeling_.five_prime_chain_break + 1 );
		core::pose::rna::remove_virtual_rna_residue_variant_type( pose, sugar_modeling_.bulge_res );

		// this used to happen inside StepWiseModeler (and not get reverted) -- special for bulge cases, so moved out here.
		reinstantiate_backbone_at_moving_res( pose, sugar_modeling_.bulge_res, sugar_modeling_.five_prime_chain_break );
	}
}

/////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::restore_pose_variants_after_chain_closure( utility::vector1< pose::PoseOP > & pose_list ){
	using namespace core::scoring;
	// reverse instantiation of bulge.
	ScoreFunctionOP modeler_scorefxn = get_modeler_scorefxn( scorefxn_ );
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		Pose & pose = *(pose_list[n]);
		remove_chain_break_variants( pose, sugar_modeling_.five_prime_chain_break );
		pose::remove_variant_type_from_pose_residue( pose,
			core::chemical::VIRTUAL_PHOSPHATE,  sugar_modeling_.five_prime_chain_break + 1 );
		core::pose::rna::apply_virtual_rna_residue_variant_type( pose, sugar_modeling_.bulge_res );
		( *modeler_scorefxn )( pose ); //for output purposes...
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::reinstantiate_backbone_at_moving_res( core::pose::Pose & pose, core::Size const rebuild_res,
	Size const five_prime_chain_break_res ) {
	pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, rebuild_res ); //May 31, 2010
	pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_O2PRIME_HYDROGEN, rebuild_res );
	if ( rebuild_res == five_prime_chain_break_res + 1 ) {
		pose::remove_variant_type_from_pose_residue( pose,
			core::chemical::VIRTUAL_PHOSPHATE, rebuild_res ); //this virtual_phosphate was added to pose at the beginning of this function.
	}
}

/////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::setup_VDW_bin_checker( pose::Pose const & input_pose ){

	core::kinematics::Stub reference_stub;
	reference_stub.v = core::chemical::rna::get_rna_base_centroid( input_pose.residue( sugar_modeling_.reference_res ), false );
	reference_stub.M = core::chemical::rna::get_rna_base_coordinate_system( input_pose.residue( sugar_modeling_.reference_res ), reference_stub.v );

	utility::vector1 < core::Size > ignore_res_list = make_vector1( sugar_modeling_.moving_res, sugar_modeling_.bulge_res, sugar_modeling_.reference_res );

	VDW_bin_checker_ = checker::RNA_VDW_BinCheckerOP( new checker::RNA_VDW_BinChecker() );
	VDW_bin_checker_->create_VDW_screen_bin( input_pose, ignore_res_list, sugar_modeling_.is_prepend, reference_stub.v, false /*verbose*/ );
}

/////////////////////////////////////////////////////////////////
// want to instantiate sugar of moving residue, but may not know
//  yet what is happening to a chunk of stuff that is attached to
//  it but distal from the reference (and will eventually be sampled).
//  This function recognizes that other stuff and virtualizes it.
void
VirtualSugarSampler::virtualize_distal_partition( pose::Pose & viewer_pose ){

	Pose pose = viewer_pose; // using a 'scratch' pose, because applying variants can confuse graphics viewers.
	distal_partition_res_.clear();
	already_virtualized_res_list_.clear();


	Size const nres = working_parameters_->working_sequence().size();
	ObjexxFCL::FArray1D < bool > partition( nres, false );
	Size const jump = pose.fold_tree().jump_nr( sugar_modeling_.moving_res, sugar_modeling_.reference_res );
	runtime_assert( jump > 0 );
	pose.fold_tree().partition_by_jump( jump, partition );

	bool const moving_res_partition    = partition( sugar_modeling_.moving_res );
	bool const reference_res_partition = partition( sugar_modeling_.reference_res );
	runtime_assert( moving_res_partition != reference_res_partition );

	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
		if ( seq_num == sugar_modeling_.moving_res ) continue;
		if ( partition( seq_num ) == moving_res_partition ) distal_partition_res_.push_back( seq_num );
	}

	// moving residues's phosphate could be part of distal partition...
	moving_phosphate_virtualized_ = false;
	if ( !pose.residue( sugar_modeling_.moving_res ).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) ) {
		pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, sugar_modeling_.moving_res );
		// If phosphate instantiated, I thought there should always be a connection to a distal 5' residue..
		// but that's not the case at a cutpoint_closed. -- rhiju.
		//  runtime_assert( distal_partition_res_.has_value( sugar_modeling_.moving_res - 1 ) )
		moving_phosphate_virtualized_ = true;
	}

	pose_with_original_terminal_phosphates_ = pose.clone();
	phosphate::remove_terminal_phosphates( pose, distal_partition_res_ );

	for ( Size ii = 1; ii <= distal_partition_res_.size(); ii++ ) {
		Size const seq_num = distal_partition_res_[ii];
		if ( pose.residue( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
			already_virtualized_res_list_.push_back( seq_num );
		} else {
			pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RNA_RESIDUE, seq_num );
		}
	}

	viewer_pose = pose;
}


/////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::reinstantiate_distal_partition( utility::vector1< PoseOP > & final_pose_list ){
	for ( Size n = 1; n <= final_pose_list.size(); n++ ) {
		pose::Pose & current_pose = ( *final_pose_list[n] );
		if ( virtual_sugar_is_from_prior_step_ ) reinstantiate_distal_partition( current_pose );
	}
}

/////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::reinstantiate_distal_partition( pose::Pose & current_pose ){

	if ( moving_phosphate_virtualized_ ) {
		pose::remove_variant_type_from_pose_residue( current_pose,
			core::chemical::VIRTUAL_PHOSPHATE, sugar_modeling_.moving_res );
	}

	for ( Size ii = 1; ii <= distal_partition_res_.size(); ii++ ) {
		Size const seq_num = distal_partition_res_[ii];
		if ( already_virtualized_res_list_.has_value( seq_num ) ) continue;
		pose::remove_variant_type_from_pose_residue( current_pose, core::chemical::VIRTUAL_RNA_RESIDUE, seq_num );
	}

	phosphate::copy_over_phosphate_variants( current_pose, *pose_with_original_terminal_phosphates_, distal_partition_res_ );
}

//////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
void
VirtualSugarSampler::reinstate_original_constraints( utility::vector1< pose::PoseOP >  & pose_list ){
	for ( Size n = 1; n <= pose_list.size(); n++ ) pose_list[ n ]->constraint_set( original_constraint_set_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
VirtualSugarSampler::fast_full_atom_VDW_repulsion_screen( core::pose::Pose const & pose, core::Size const res_1, core::Size const res_2, bool const is_prepend ){

	conformation::Residue const & rsd_1 = pose.residue( res_1 );
	conformation::Residue const & rsd_2 = pose.residue( res_2 );


	for ( Size n_1 = 1; n_1 <= rsd_1.natoms(); n_1++ ) {

		//atom 1-4 are " P  ", " OP2", " OP1" and " O5'"
		Size const act_res_1 = ( is_prepend && n_1 <= 4 ) ? res_1 + 1: res_1;

		if ( pose.residue( act_res_1 ).is_virtual( n_1 )  ) continue;

		for ( Size n_2 = 1; n_2 <= rsd_2.natoms(); n_2++ ) {

			Size const act_res_2 = ( is_prepend && n_2 <= 4 ) ? res_2 + 1: res_2;

			if ( pose.residue( act_res_2 ).is_virtual( n_2 )  ) continue;

			Real const VDW_radius_1 = pose.residue( act_res_1 ).atom_type( n_1 ).lj_radius();
			Real const VDW_radius_2 = pose.residue( act_res_2 ).atom_type( n_2 ).lj_radius();

			Real const clash_dist_cutoff = 0.8; //Fail van der Waals repulsion screen if two atoms radius within 0.5 Angstrom of each other

			Real const clash_radius = VDW_radius_1 + VDW_radius_2 - clash_dist_cutoff;

			if ( ( pose.residue( act_res_1 ).xyz( n_1 ) - pose.residue( act_res_2 ).xyz( n_2 ) ).length_squared() < clash_radius*clash_radius ) {
				return false; //OK consider fail screening if find even one clash...
			}

		}
	}
	return true;
}

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols
