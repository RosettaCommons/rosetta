// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_VirtualSugarSampler
/// @brief
/// @detailed
/// @author Parin Sripakdeevong, Rhiju Das (rhiju@stanford.edu)


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh> // for nice, encapsulated chain closer.
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/rna/screener/AtrRepScreener.hh>
#include <protocols/swa/rna/screener/ChainClosableScreener.hh>
#include <protocols/swa/rna/screener/ChainBreakScreener.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/rotamer_sampler/rna/RNA_SuiteRotamer.hh>
#include <protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

using namespace core;
using utility::tools::make_vector1;

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_VirtualSugarSampler" );

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

namespace protocols {
namespace swa {
namespace rna {

//Constructor
StepWiseRNA_VirtualSugarSampler::StepWiseRNA_VirtualSugarSampler( StepWiseRNA_JobParametersCOP & job_parameters,
																																	SugarModeling & sugar_modeling	):
	job_parameters_( job_parameters ),
	sugar_modeling_( sugar_modeling ), // holds both job parameters and pose list.
	tag_( "" ),
	use_phenix_geo_( false ),
	integration_test_mode_( false ),
	virtual_sugar_is_from_prior_step_( true ),
	legacy_mode_( false ), //for comparison to original code -- may deprecate in 2014.
	keep_base_fixed_( false ),
	do_chain_closure_( true || legacy_mode_ ),
	first_minimize_with_fixed_base_( false ),
	max_tries_for_random_overall_( 12 ),
	max_tries_for_random_sugar_setup_( 1 ),
	sugar_setup_success_( false ),
	moving_phosphate_virtualized_( false )
{
}

//Destructor
StepWiseRNA_VirtualSugarSampler::~StepWiseRNA_VirtualSugarSampler()
{}

/////////////////////
std::string
StepWiseRNA_VirtualSugarSampler::get_name() const {
	return "StepWiseRNA_VirtualSugarSampler";
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::apply( pose::Pose & viewer_pose ){

	utility::vector1< PoseOP > pose_list;

	Size const max_tries = choose_random_ ? max_tries_for_random_sugar_setup_ : 1;
	for ( Size ntries = 1; ntries <= max_tries; ntries++ ){
		TR << TR.Green << "On sampling try: " << ntries << TR.Reset << std::endl;
		TR << TR.Green << "Will sample residue " << sugar_modeling_.moving_res << " -- reference res: " << sugar_modeling_.reference_res <<  TR.Reset << std::endl;
		TR << TR.Green << viewer_pose.fold_tree() << TR.Reset << std::endl;
		apply( viewer_pose, pose_list );
		if ( pose_list.size() > 0 || !sugar_setup_success_ ) break;
	}
	sugar_modeling_.pose_list = pose_list; // stores the solutions.

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////July 21, 2011 Moved from StepWiseRNA_ResidueSampler.cc ////////////////////////////////////////////////
////////Code used to be part of the function StepWiseRNA_ResidueSampler::previous_floating_base_chain_closure();

// RD 2013 -- following is Parin's warning. I think it is deprecated now, all cases are handled... no need to back-recurse through all virtual sugars.
//Dec 9, 2010
//Warning, currently the code doesn't handle the case of i-1, i-3 (and potentially i-5 and so on) bulges very well. Basically will sample the sugar conformation of base i. But even though base i-2, i-4 have virtual sugars, these sugars are not sampled (use the one that currently exist in the pose) which is NOT a complete fail, since the sugar did pass distance check of chain_closable_screener (correctly implement this on Dec 9,2010) and VDW screening.
//Possible ways to fix this (HOWEVER lets not fix this until we actual find a REAL existing case where i-1 and i-3 are bulge res):

//Will have to at least rewrite both floating_base_chain_closure_setup() and floating_base_chain_closure_sampling() to sample both i and i-2 sugar.
//Hardest 1. To change this would have to recursively sample the conformation of the bulge and corresponding virtual sugar starting from the bulge that is futhest away. If assume that each bulge conformation is independent of the other bulges then this is not very computationally expensive. Although this would require sampling two sugar for each bulge except for the bulge that is futhest again.   A trace back is then needed to find combinations of all the bulge conformation which is closable.
//Medium 2. Another possibility is to simple assume that i-2 base can assume all bulge conformation subjected to chain_closable_screener.
//Easiest 3. Just set chain_closure_sampling to false in this function. This will just keep all possible bulge conformation in pose_list that pass O3i_C5iplus2_distance>O3I_C5I_PLUS_TWO_MAX_DIST (in setup)
//However in Easiest 3. case, there is a problem is problem in that O3i_C5iplus2_distance has determined with fixed i-2 sugar.
void
StepWiseRNA_VirtualSugarSampler::apply( pose::Pose & viewer_pose,
																				utility::vector1< PoseOP > & pose_list ){
	using namespace ObjexxFCL;
	using namespace core::io::silent;

	output_title_text( "Enter sample_virtual_sugar_and_bulge_and_close_chain()", TR.Debug );
	clock_t const time_start( clock() );
	pose::Pose const viewer_pose_save = viewer_pose;

	tag_into_pose( viewer_pose, tag_ );
	sugar_modeling_.check_compatibility( viewer_pose.total_residue() );
	original_constraint_set_ = viewer_pose.constraint_set()->clone();

	//Virtualize any residues that move 'with' the moving residue -- their locations have not been sampled yet.
	if ( virtual_sugar_is_from_prior_step_ ) virtualize_distal_partition( viewer_pose );
	else runtime_assert( job_parameters_->gap_size() > 0 );

	setup_VDW_bin_screener( viewer_pose );

	setup_sugar_conformations( pose_list, viewer_pose );

	if ( ( sugar_modeling_.bulge_res > 0 ) && do_chain_closure_ ) do_chain_closure_sampling( pose_list, viewer_pose );

	reinstate_original_constraints( pose_list );
	if ( virtual_sugar_is_from_prior_step_ )	reinstantiate_distal_partition( pose_list );

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
StepWiseRNA_VirtualSugarSampler::setup_sugar_conformations( utility::vector1< PoseOP > & pose_list, pose::Pose & pose ){

	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	using namespace core::id;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring::rna;
	using namespace rotamer_sampler::rna;
	using namespace protocols::rna;
	using namespace core::pose;

	clock_t const time_start( clock() );
	//July 28th, 2011 Could set extra_chi here, BUT necessary?
	RNA_NucleosideRotamer sampler( sugar_modeling_.moving_res,
																 sugar_modeling_.moving_res_pucker_state,
																 sugar_modeling_.moving_res_base_state );
	sampler.set_idealize_coord(   use_phenix_geo_ );
	sampler.set_skip_same_pucker( use_phenix_geo_ );
	sampler.set_random( choose_random_ );
	sampler.init();

	Size const moving_suite = sugar_modeling_.is_prepend ? sugar_modeling_.moving_res : ( sugar_modeling_.moving_res - 1);
	bool const do_minimize = !integration_test_mode_; // for speed.
	Size const five_prime_gap_res  = sugar_modeling_.is_prepend ? sugar_modeling_.moving_res    : sugar_modeling_.reference_res;
	Size const three_prime_gap_res = sugar_modeling_.is_prepend ? sugar_modeling_.reference_res : sugar_modeling_.moving_res;

	// Deprecate soon:
	//	TR << TR.Red << "MOVING " << sugar_modeling_.moving_res << "  BULGE " << sugar_modeling_.bulge_res << "  REFERENCE " << sugar_modeling_.reference_res << " IS_PREPEND " << sugar_modeling_.is_prepend << " MOVING_SUITE " << moving_suite << " sugar_modeling:FIVE_PRIME_CHAINBREAK " << sugar_modeling_.five_prime_chain_break << " FIVE_PRIME_GAP_RES " << five_prime_gap_res << " DO_MINIMIZE " << do_minimize << TR.Reset << std::endl;

	screener::ChainClosableScreener chain_closable_screener( five_prime_gap_res, three_prime_gap_res, 1 /*gap_size -- problem!?*/);
	screener::AtrRepScreener atr_rep_screener( pose, sugar_modeling_.moving_res, sugar_modeling_.reference_res, 0 /*gap_size = 0 forces loose rep_cutoff (10.0 )*/ );

	pose::Pose pose_init = pose;//Hard copy
	pose::remove_variant_type_from_pose_residue( pose_init, "VIRTUAL_RIBOSE", sugar_modeling_.moving_res );

	Size count( 0 );
	sugar_setup_success_ = false;
	pose_list.clear();
	for ( sampler.reset(); sampler.not_end(); ++sampler ) {
		count++;
		if ( choose_random_ && count > max_tries_for_random_sugar_setup_ ){
			pose_list.push_back( pose.clone() ); // this routine must return something with an instantiated sugar, even if its a fail.
			break;
		}

		pose = pose_init;
		// following *must* be mistake... should occur below inside actual chain_closure..
		if ( legacy_mode_) {
			add_harmonic_chain_break_constraint( pose, sugar_modeling_.five_prime_chain_break );
			if  ( do_chain_closure_ ) setup_chain_break_variants( pose, sugar_modeling_.five_prime_chain_break );
		} else {
			add_fade_chain_break_constraint_across_gap( pose, five_prime_gap_res, three_prime_gap_res,	1 /*gap_size -- one intervening bulge residue -- perhaps should not be hard-wired! */ );
		}

		sampler.apply( pose );

		if ( do_minimize ) minimize_sugar( pose ); // remove clashes; show in graphics too.
		if ( ! chain_closable_screener.check_screen( pose ) ) continue;
		if ( ! atr_rep_screener.check_screen( pose ) ) continue;

		tag_into_pose( pose, tag_from_pose( pose_init ) + '_' + ObjexxFCL::string_of( count ) );
		pose_list.push_back( pose.clone() );
		sugar_setup_success_ = true;
		if ( choose_random_ ) break; // take the first random sample.
	}

	output_title_text( "", TR.Debug );
	TR.Debug << " pose_list.size() = " << pose_list.size() << std::endl;
	TR.Debug << "Total time in Floating_base_chain_closure SETUP: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	pose = pose_init;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::minimize_sugar( pose::Pose & pose_with_sugar ){

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
	MinimizerOptions options_standard( "dfpmin_atol", tolerance, use_nblist, false, false );      //Switch to absolute tolerance on Sept 21, 2011
	MinimizerOptions options_armijo( "dfpmin_armijo_atol", tolerance, use_nblist, false, false ); //Add this on Sept 21, 2011

	TR.Debug << "options_standard: min_type = " << options_standard.min_type() << " minimize_tolerance = "  << options_standard.minimize_tolerance() << std::endl;
	TR.Debug << "options_armijo  : min_type = " << options_armijo.min_type()   << " minimize_tolerance = "  << options_armijo.minimize_tolerance()   << std::endl;

	options_standard.nblist_auto_update( true );
	options_armijo.nblist_auto_update( true );

	MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	if ( first_minimize_with_fixed_base_ || keep_base_fixed_ ){
		mm.set( TorsionID( sugar_modeling_.moving_res, id::CHI, 1 ), true ); // nucleosidic chi
		mm.set( TorsionID( sugar_modeling_.moving_res, id::CHI, 2 ), true ); // sugar torsion
		mm.set( TorsionID( sugar_modeling_.moving_res, id::CHI, 3 ), true ); // sugar torsion
		minimizer.run( pose_with_sugar, mm, *( sugar_scorefxn_without_ch_bond ), options_standard );
	}

	if ( keep_base_fixed_ ) return;

	bool found_desired_jump_ID = false;
	for ( Size jump_ID = 1; jump_ID <= pose_with_sugar.fold_tree().num_jump(); jump_ID++ ){
		Size const jump_pos1( pose_with_sugar.fold_tree().upstream_jump_residue( jump_ID ) );
		Size const jump_pos2( pose_with_sugar.fold_tree().downstream_jump_residue( jump_ID ) );
		Size const cutpoint = pose_with_sugar.fold_tree().cutpoint( jump_ID );
		if ( ( jump_pos1 == sugar_modeling_.moving_res ) || ( jump_pos2 == sugar_modeling_.moving_res ) ){
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

} // do_minimize.


/////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::get_sugar_setup_scorefxns( scoring::ScoreFunctionOP & sugar_scorefxn, scoring::ScoreFunctionOP & sugar_scorefxn_without_ch_bond, scoring::ScoreFunctionOP & rescaled_sugar_score_fxn_without_ch_bond ) const {

	using namespace core::scoring;

	sugar_scorefxn =	scorefxn_->clone();
	sugar_scorefxn->set_weight( linear_chainbreak, 0.0 );
	sugar_scorefxn->set_weight( angle_constraint, 0.0 );
	sugar_scorefxn->set_weight( atom_pair_constraint, 0.0 );
	sugar_scorefxn->set_weight( coordinate_constraint, 0.1 );

	sugar_scorefxn_without_ch_bond =	sugar_scorefxn->clone();
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
StepWiseRNA_VirtualSugarSampler::do_chain_closure_sampling( utility::vector1< PoseOP > & pose_list,
																														pose::Pose & viewer_pose ){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring::rna;
	using namespace ObjexxFCL;

	if ( pose_list.size() == 0 ) return; //early return

	initialize_pose_variants_for_chain_closure( pose_list );
	floating_base_chain_closure( pose_list, viewer_pose );
	restore_pose_variants_after_chain_closure( pose_list );

}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::initialize_pose_variants_for_chain_closure( utility::vector1< pose::PoseOP > & pose_list ){
	// instantiate bulge.
	for ( Size n = 1; n <= pose_list.size(); n++ ){
		Pose & pose = *(pose_list[n]);
		// in legacy mode this occurs above in 'setup'. but it really should be here.
		if ( !legacy_mode_)	setup_chain_break_variants( pose, sugar_modeling_.five_prime_chain_break );
		remove_virtual_rna_residue_variant_type( pose, sugar_modeling_.bulge_res );
		pose::add_variant_type_to_pose_residue(  pose, "VIRTUAL_PHOSPHATE", sugar_modeling_.five_prime_chain_break + 1 );
	}
}

/////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::restore_pose_variants_after_chain_closure( utility::vector1< pose::PoseOP > & pose_list ){
	using namespace core::scoring;
	// reverse instantiation of bulge.
	ScoreFunctionOP sampling_scorefxn = get_sampling_scorefxn( scorefxn_ );
	for ( Size n = 1; n <= pose_list.size(); n++ ){
		Pose & pose = *(pose_list[n]);
		remove_chain_break_variants( pose, sugar_modeling_.five_prime_chain_break );
		apply_virtual_rna_residue_variant_type( pose, sugar_modeling_.bulge_res );
		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE",  sugar_modeling_.five_prime_chain_break + 1 );
		( *sampling_scorefxn )( pose ); //for output purposes...
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::floating_base_chain_closure( utility::vector1< PoseOP > & pose_list,
																															pose::Pose & viewer_pose ){
	if ( legacy_mode_ ){ // deprecate legacy code soon
		floating_base_chain_closure_legacy( pose_list, viewer_pose );
		floating_base_chain_minimize_legacy( pose_list, viewer_pose );
	} else {
		floating_base_chain_closure_complete( pose_list, viewer_pose );
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
// through the use of KIC-based sampling, this bulge modeling is both more comprehensive and
// faster than the legacy code [which was CCD-based].
void
StepWiseRNA_VirtualSugarSampler::floating_base_chain_closure_complete( utility::vector1< PoseOP > & pose_list,
																																			 pose::Pose & viewer_pose ){
	utility::vector1< PoseOP > output_pose_list;
	StepWiseRNA_ModelerOP stepwise_rna_modeler = new StepWiseRNA_Modeler( sugar_modeling_.bulge_res,
																																				scorefxn_ );
	stepwise_rna_modeler->set_use_phenix_geo( use_phenix_geo_ );
	stepwise_rna_modeler->set_kic_sampling_if_relevant( true /*erraser sampling for speed & completeness*/ );
	stepwise_rna_modeler->set_num_pose_minimize( 1 );
	stepwise_rna_modeler->set_integration_test_mode( true /*fast sampling*/ );
	for ( Size n = 1; n <= pose_list.size(); n++ ){
		Pose & pose = *(pose_list[n]);
		stepwise_rna_modeler->apply( pose);
		if ( stepwise_rna_modeler->get_num_sampled() > 0 ) output_pose_list.push_back( pose_list[n] );
	}
	pose_list = output_pose_list;
}

/////////////////////////////////////////////////////////////////////////////////////////
// deprecate this after march 2014 if new 'floating_base_chain_closure_complete' is
// satisfactory for stepwise assembly & stepwise monte carlo runs.
void
StepWiseRNA_VirtualSugarSampler::floating_base_chain_closure_legacy( utility::vector1< PoseOP > & pose_list,
																																		 pose::Pose & viewer_pose ){

	using namespace rotamer_sampler::rna;
	using namespace protocols::swa::rna::screener;
	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	using namespace protocols::rna;

	pose::Pose screening_pose = *(pose_list[1]);

	output_title_text( "Floating_base_chain_closure SAMPLING", TR.Debug );
	clock_t const time_start_sampling( clock() );

	// in legacy mode chainbreak constraint addition occurs above in 'setup'. but it really should be here.
	if ( !legacy_mode_ ){
		for ( Size n = 1; n <= pose_list.size(); n++ ){
			Pose & pose = *(pose_list[n]);
			add_harmonic_chain_break_constraint( pose, sugar_modeling_.five_prime_chain_break );
		}
	}

	Size const bulge_suite = sugar_modeling_.bulge_suite;
	Size const bulge_rsd   = sugar_modeling_.bulge_res;
	utility::vector1<bool> sample_sugar( 2, false );
	utility::vector1<Size> base_state  ( 2, WHATEVER );
	utility::vector1<Size> pucker_state( 2, WHATEVER );

	for ( Size i = 1; i <= 2; ++i ) {
		Size const curr_rsd( bulge_suite + i - 1 );
		if ( bulge_rsd == curr_rsd ) {
			sample_sugar[i] = true;
			pucker_state[i] = WHATEVER;
			base_state[i] = WHATEVER;
		} else {
			sample_sugar[i] = false;
			//FANG: Assuming all pose in the list have same pucker for the unmoved
			//residue. Not sure if this is correct but it only affects
			//epsilon sampling.
			pucker_state[i] = assign_pucker( screening_pose, curr_rsd ); // NOT TRUE! FIX!
			base_state[i]   = NONE;
		}
	}

	RNA_SuiteRotamer sampler( bulge_suite, pucker_state[1], pucker_state[2], base_state[1], base_state[2] );
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

	utility::vector1< AtrRepScreenerOP > atr_rep_screeners_;
	for ( Size n = 1; n <= pose_list.size(); n++ )	atr_rep_screeners_.push_back( new AtrRepScreener( *(pose_list[n]), bulge_suite, bulge_rsd, 0 /*gap_size*/ ) );

	ChainClosableScreenerOP chain_closable_screener_ = new ChainClosableScreener( sugar_modeling_.five_prime_chain_break, 0 /*gap_size*/ );
	ChainBreakScreenerOP chain_break_screener_       = new ChainBreakScreener( screening_pose, sugar_modeling_.five_prime_chain_break );
	chain_break_screener_->set_reinitialize_CCD_torsions( true );

	utility::vector1< bool > is_chain_close( pose_list.size(), false );
	utility::vector1< pose::PoseOP > output_pose_list;

	// following is not a rigorous enumeration and kind of slow.
	// first, the sampler is not subdivided by pucker of the actual pose.
	// And... each backbone torsion is checked once, looking for any pose, rather
	// than trying all backbone against all poses.
	for ( sampler.reset(); sampler.not_end(); ++sampler ) {
		if ( integration_test_mode_ && num_closed_chain_pose > 1 ) break;

		if ( num_closed_chain_pose == pose_list.size() )	break;

		count_data.tot_rotamer_count++;
		sampler.apply( screening_pose );

		if ( !VDW_bin_screener_->VDW_rep_screen( screening_pose, sugar_modeling_.bulge_res ) ) continue;
		count_data.good_bin_rep_count++;

		for ( Size n = 1; n <= pose_list.size(); n++ ){

			if ( is_chain_close[ n ] ) continue;
			Pose & current_pose = *(pose_list[n]);

			// don't copy torsions into current pose just quite yet -- check if chain is closable.
			if  ( !chain_closable_screener_->check_screen( screening_pose   /* moving_pose */,
																										 current_pose     /* reference_pose */,
																										 !sugar_modeling_.is_prepend,
																										 true /*strict*/ ) ) continue;
			count_data.chain_closable_count++;

			sampler.apply( current_pose );

			if ( !fast_full_atom_VDW_repulsion_screen( current_pose,
																								 sugar_modeling_.bulge_res, sugar_modeling_.moving_res,
																								 sugar_modeling_.is_prepend ) ) continue;
			count_data.fast_full_atom_VDW_repulsion_screen++;

			if ( !atr_rep_screeners_[n]->check_screen( current_pose ) ) continue;
			if ( chain_break_screener_->check_screen( current_pose ) ) continue;

			is_chain_close[ n ] = true;
			num_closed_chain_pose++;
			viewer_pose = current_pose; // show in graphics.
			output_pose_list.push_back( pose_list[n] );
			if ( integration_test_mode_ ) break;
		}
	}

	TR.Debug << " bin_rep_count = " << count_data.good_bin_rep_count;
	TR.Debug << " fast_rep_count = " << count_data.fast_full_atom_VDW_repulsion_screen;
	TR.Debug << " chain_closable = " << count_data.chain_closable_count;
	TR.Debug << " angle_n = " << count_data.good_angle_count << " dist_n = " << count_data.good_distance_count;
	TR.Debug << " rep = " << count_data.good_rep_rotamer_count;
	TR.Debug << " rmsd = " << count_data.rmsd_count << " tot = " << count_data.tot_rotamer_count << std::endl;
	TR.Debug << " " << num_closed_chain_pose << " out of " << pose_list.size() << " pose were closable" << std::endl;
	TR.Debug << "Total time in Floating_base_chain_closure SAMPLING: " << static_cast< Real > ( clock() - time_start_sampling ) / CLOCKS_PER_SEC << std::endl;

	pose_list = output_pose_list;
}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::floating_base_chain_minimize_legacy( utility::vector1< PoseOP > & pose_list,
																																			pose::Pose & viewer_pose ){

	using namespace core::optimization;
	using namespace core::id;
	using namespace core::scoring;

	////////////////////////////////////////////////////////////////
	// post-process! -- why not just use existing Minimizer class?
	///////////////////////////////////////////////////////////////
	output_title_text( "Floating_base_chain_closure POST_PROCESSING", TR.Debug );

	//Quick minimize to reduce error in CCD///////////////////////////
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025 );
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
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

	ScoreFunctionOP bulge_chain_closure_scorefxn = new ScoreFunction;
	bulge_chain_closure_scorefxn->set_weight( fa_rep , 0.12 );
	bulge_chain_closure_scorefxn->set_weight( angle_constraint, 1.0 );
	bulge_chain_closure_scorefxn->set_weight( atom_pair_constraint, 1.0 );
	bulge_chain_closure_scorefxn->set_weight( linear_chainbreak, 5.0 );
	bulge_chain_closure_scorefxn->set_weight( rna_sugar_close, 70 ); //100X normal weight

	mm.set( TorsionID( sugar_modeling_.bulge_res, id::CHI, 2 ), true ); //nu2
	mm.set( TorsionID( sugar_modeling_.bulge_res, id::CHI, 3 ), true ); //nu1

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< PoseOP > output_pose_list;
	for ( Size n = 1; n <= pose_list.size(); n++ ){
		viewer_pose = *(pose_list[n]);
		TR.Debug << "POST_PROCESSING pose # " << n << " out of " << pose_list.size() << " " << std::endl;;
		minimizer.run( viewer_pose, mm, *( bulge_chain_closure_scorefxn ), options );
		output_pose_list.push_back( viewer_pose.clone() );
	}
	pose_list = output_pose_list;
}


/////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::setup_VDW_bin_screener( pose::Pose const & input_pose ){

	core::kinematics::Stub reference_stub;
	reference_stub.v = core::chemical::rna::get_rna_base_centroid( input_pose.residue( sugar_modeling_.reference_res ), false );
	reference_stub.M = core::chemical::rna::get_rna_base_coordinate_system( input_pose.residue( sugar_modeling_.reference_res ), reference_stub.v );

	utility::vector1 < core::Size > ignore_res_list = make_vector1( sugar_modeling_.moving_res, sugar_modeling_.bulge_res, sugar_modeling_.reference_res );

	VDW_bin_screener_ = new screener::StepWiseRNA_VDW_BinScreener();
	VDW_bin_screener_->create_VDW_screen_bin( input_pose, ignore_res_list, sugar_modeling_.is_prepend, reference_stub.v, false /*verbose*/ );
}

/////////////////////////////////////////////////////////////////
// want to instantiate sugar of moving residue, but may not know
//  yet what is happening to a chunk of stuff that is attached to
//  it but distal from the reference (and will eventually be sampled).
//  This function recognizes that other stuff and virtualizes it.
void
StepWiseRNA_VirtualSugarSampler::virtualize_distal_partition( pose::Pose & viewer_pose ){

	Pose pose = viewer_pose; // using a 'scratch' pose, because applying variants can confuse graphics viewers.
	distal_partition_pos_.clear();
	already_virtualized_res_list_.clear();
	Size const nres = job_parameters_->working_sequence().size();

	ObjexxFCL::FArray1D < bool > partition( nres, false );
	Size const jump = pose.fold_tree().jump_nr( sugar_modeling_.moving_res, sugar_modeling_.reference_res );
	runtime_assert( jump > 0 );
	pose.fold_tree().partition_by_jump( jump, partition );

	bool const moving_res_partition    = partition( sugar_modeling_.moving_res );
	bool const reference_res_partition = partition( sugar_modeling_.reference_res );
	runtime_assert( moving_res_partition != reference_res_partition );

	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
		if ( seq_num == sugar_modeling_.moving_res ) continue;
		if ( partition( seq_num ) == moving_res_partition ) distal_partition_pos_.push_back( seq_num );
	}

	// moving residues's phosphate could be part of distal partition...
	moving_phosphate_virtualized_ = false;
	if ( !pose.residue( sugar_modeling_.moving_res ).has_variant_type( "VIRTUAL_PHOSPHATE" ) &&
			 !pose.residue( sugar_modeling_.moving_res ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" /*same as virtual phosphate*/ ) ){
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", sugar_modeling_.moving_res );
		// If phosphate instantiated, I thought there should always be a connection to a distal 5' residue..
		// but that's not the case at a cutpoint_closed. -- rhiju.
		//		runtime_assert( distal_partition_pos_.has_value( sugar_modeling_.moving_res - 1 ) )
		moving_phosphate_virtualized_ = true;
	}

	for ( Size ii = 1; ii <= distal_partition_pos_.size(); ii++ ){
		Size const seq_num = distal_partition_pos_[ii];
		if ( pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
			already_virtualized_res_list_.push_back( seq_num );
		} else {
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", seq_num );
		}
	}

	//	if ( job_parameters_->gap_size() == 0 ){
	//		runtime_assert ( !pose.residue( job_parameters_->five_prime_chain_break_res() + 1 ).has_variant_type( "VIRTUAL_PHOSPHATE" ) );
	//	}
	viewer_pose = pose;
}


/////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::reinstantiate_distal_partition( utility::vector1< PoseOP > & final_pose_list ){
	for ( Size n = 1; n <= final_pose_list.size(); n++ ){
		pose::Pose & current_pose = ( *final_pose_list[n] );
		if ( virtual_sugar_is_from_prior_step_ )	reinstantiate_distal_partition( current_pose );
	}
}

/////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::reinstantiate_distal_partition( pose::Pose & current_pose ){

	if ( moving_phosphate_virtualized_ ) pose::remove_variant_type_from_pose_residue( current_pose, "VIRTUAL_PHOSPHATE", sugar_modeling_.moving_res );

	for ( Size ii = 1; ii <= distal_partition_pos_.size(); ii++ ){
		Size const seq_num = distal_partition_pos_[ii];
		if ( already_virtualized_res_list_.has_value( seq_num ) ) continue;
		pose::remove_variant_type_from_pose_residue( current_pose, "VIRTUAL_RNA_RESIDUE", seq_num );
	}
	//	if ( job_parameters_->gap_size() == 0 ) {
	//	//		pose::remove_variant_type_from_pose_residue( current_pose, "VIRTUAL_PHOSPHATE", job_parameters_->five_prime_chain_break_res() + 1 );
	//	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarSampler::reinstate_original_constraints( utility::vector1< pose::PoseOP >  & pose_list ){
	for ( Size n = 1; n <= pose_list.size(); n++ ) pose_list[ n ]->constraint_set( original_constraint_set_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarSampler::fast_full_atom_VDW_repulsion_screen( core::pose::Pose const & pose, core::Size const res_1, core::Size const res_2, bool const is_prepend ){

	conformation::Residue const & rsd_1 = pose.residue( res_1 );
	conformation::Residue const & rsd_2 = pose.residue( res_2 );


	for ( Size n_1 = 1; n_1 <= rsd_1.natoms(); n_1++ ){

		//atom 1-4 are " P  ", " OP2", " OP1" and " O5'"
		Size const act_res_1 = ( is_prepend && n_1 <= 4 ) ? res_1 + 1: res_1;

		if ( pose.residue( act_res_1 ).is_virtual( n_1 )  ) continue;

		for ( Size n_2 = 1; n_2 <= rsd_2.natoms(); n_2++ ){

			Size const act_res_2 = ( is_prepend && n_2 <= 4 ) ? res_2 + 1: res_2;

			if ( pose.residue( act_res_2 ).is_virtual( n_2 )  ) continue;

			Real const VDW_radius_1 = pose.residue( act_res_1 ).atom_type( n_1 ).lj_radius();
			Real const VDW_radius_2 = pose.residue( act_res_2 ).atom_type( n_2 ).lj_radius();

			Real const clash_dist_cutoff = 0.8; //Fail van der Waals repulsion screen if two atoms radius within 0.5 Angstrom of each other

			Real const clash_radius = VDW_radius_1 + VDW_radius_2 - clash_dist_cutoff;

			if ( ( pose.residue( act_res_1 ).xyz( n_1 ) - pose.residue( act_res_2 ).xyz( n_2 ) ).length_squared() < clash_radius*clash_radius ){
				return false; //OK consider fail screening if find even one clash...
			}

		}
	}
	return true;
}

} //rna
} //swa
} // protocols

