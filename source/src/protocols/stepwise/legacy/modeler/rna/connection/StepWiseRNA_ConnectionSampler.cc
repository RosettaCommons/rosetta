// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/rigid_body/StepWiseRNA_ConnectionSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/legacy/modeler/rna/connection/StepWiseRNA_ConnectionSampler.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.hh>
#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/rigid_body/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/StepWiseRNA_VirtualSugarJustInTimeInstantiator.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/legacy/screener/RNA_AtrRepScreener.hh>
#include <protocols/stepwise/screener/BaseBinMapUpdater.hh>
#include <protocols/stepwise/screener/BaseCentroidScreener.hh>
#include <protocols/stepwise/screener/BulgeApplier.hh>
#include <protocols/stepwise/screener/CentroidDistanceScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosableGeometryScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosableGeometryResidueBasedScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosureScreener.hh>
#include <protocols/stepwise/screener/FastForwardToNextRigidBody.hh>
#include <protocols/stepwise/screener/FastForwardToNextResidueAlternative.hh>
#include <protocols/stepwise/screener/IntegrationTestBreaker.hh>
#include <protocols/stepwise/screener/NativeRMSD_Screener.hh>
#include <protocols/stepwise/screener/O2PrimeScreener.hh>
#include <protocols/stepwise/screener/PhosphateScreener.hh>
#include <protocols/stepwise/screener/PoseSelectionScreener.hh>
#include <protocols/stepwise/screener/ResidueContactScreener.hh>
#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/StubApplier.hh>
#include <protocols/stepwise/screener/StubDistanceScreener.hh>
#include <protocols/stepwise/screener/SugarInstantiator.hh>
#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/stepwise/screener/VDW_BinScreener.hh>
#include <protocols/stepwise/StepWiseSampleAndScreen.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/rna/util.hh>
#include <protocols/stepwise/sampler/rigid_body/util.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/types.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

using namespace core;
using ObjexxFCL::string_of;
using ObjexxFCL::lead_zero_string_of;

static basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.connection.StepWiseRNA_ConnectionSampler" ) ;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Enumerates over suite backbone, sugar pucker, and base chi angle; and packs 2' hydroxyls.
//    or...
// enumerate rigid body degrees of freedom for a chunk of the pose.
//
// Decides which one to do based on what residue is upstream of the supplied moving_residue.
//
//
//        Reference residue        Moving residue                 Distal residue
//                                (virtual*)   (virtual*)        (virtual*)
//     5'   --Sugar -- ...       -- Phosphate -- Sugar -- ... -- Phosphate -- Sugar -- Phosphate -- 3'
//              |                       |                           |
//          Reference                Moving                      Distal
//             Base                   Base!                       Base
//              |_______________________________________________|
//                         jump or suite connection
//
//
//  * may be virtual -- the exceptions are if connected through chain closure to reference and/or distal residue.
//
// Note that sugar of the floating base will be virtualized unless it needs to be instantiated to close chain to
//  anchor or distal residues. Nevertheless, code below carries out a geometric 'sanity check' that involves temporarily
//  instantiating the sugar (within the screening_pose).
//
// The Sample-and-screen uses several poses and checkers to reduce the number of variant changes &
// energy recomputations -- faster code but somewhat complicated. Note that some information on, say, multiple
//  options for a sugar can be passed into  the sampler in 'residue_alternative_set'; some of the first screens
//  are based on base centroid distances and are independent of sugar conformation, and so some computation
//  that would be redundant across anchor sugar conformations can be avoided.
//
//   -- Rhiju, 2014.
//    [after massive refactoring of code from Parin Sripakdeevong & Rhiju Das, 2009-2013]
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {
namespace connection {

//Constructor
StepWiseRNA_ConnectionSampler::StepWiseRNA_ConnectionSampler( stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP & working_parameters ):
	moving_res_( working_parameters->working_moving_res() ),
	reference_res_( 0 ), // updated below.
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "farna/rna_hires.wts" ) ), // can be replaced from the outside
	silent_file_( "silent_file.txt" ),
	kic_modeler_( false ), // updated below
	rebuild_bulge_mode_( working_parameters->rebuild_bulge_mode() ),
	max_distance_squared_( 0.0 ), // updated below
	rigid_body_modeler_( false ), // updated below.
	try_sugar_instantiation_( false ),
	o2prime_instantiation_distance_cutoff_( 6.0 ),
	extra_tag_( "" ),
	virt_sugar_atr_rep_screen_( false ),
	build_pose_from_scratch_( working_parameters->working_sequence().length() == ( working_parameters->working_moving_res_list().size() + 1 ) ), // somewhat hacky, used for rna puzzle
	working_parameters_( working_parameters ) // may deprecate
{
	set_native_pose( working_parameters->working_native_pose() );
}

//Destructor
StepWiseRNA_ConnectionSampler::~StepWiseRNA_ConnectionSampler()
{}

/////////////////////
std::string
StepWiseRNA_ConnectionSampler::get_name() const {
	return "StepWiseRNA_ConnectionSampler";
}

///////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::apply( core::pose::Pose & pose ){

	stepwise::modeler::rna::output_title_text( "Enter StepWiseRNA_ConnectionSampler::floating_base_modeler", TR.Debug );

	clock_t const time_start( clock() );
	figure_out_reference_res( pose );
	if ( !initialize_poses_and_checkers( pose ) ) return;
	initialize_sampler();
	initialize_screeners( pose );
	check_working_parameters( pose );

	bool const verbose = ( !options_->choose_random() || options_->integration_test_mode() );

	StepWiseSampleAndScreen sample_and_screen( sampler_, screeners_ );
	sample_and_screen.set_max_ntries( get_max_ntries() );
	sample_and_screen.set_num_random_samples( options_->num_random_samples() );

	if ( verbose ) TR << "Running SampleAndScreen... " << std::endl;
	sample_and_screen.run();

	if ( verbose ) sample_and_screen.output_counts();
	sample_and_screen.output_info_on_random_trials();

	pose_selection_->finalize( !build_pose_from_scratch_ /*do_clustering*/ );
	pose_list_ = pose_selection_->pose_list();
	if ( rigid_body_modeler_ && options_->verbose() ) analyze_base_bin_map( base_bin_map_, "test/" );

	TR.Debug << "Sampling time: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::figure_out_reference_res( pose::Pose const & pose ){

	Size moving_res = moving_res_;
	rigid_body_modeler_ = ( sampler::rigid_body::figure_out_reference_res_for_jump( pose, moving_res ) > 0 );
	if ( rigid_body_modeler_ ){
		figure_out_reference_res_with_rigid_body_rotamer( pose );
	} else {
		reference_res_ = figure_out_reference_res_for_suite( pose, moving_res );
		moving_partition_res_ = figure_out_moving_partition_res_for_suite( pose, moving_res, reference_res_ );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::figure_out_reference_res_with_rigid_body_rotamer( pose::Pose const & pose ){
	rigid_body_rotamer_ = new sampler::rigid_body::RigidBodyStepWiseSampler( pose, moving_res_ );
	reference_res_ = rigid_body_rotamer_->reference_res();
	moving_partition_res_ = rigid_body_rotamer_->moving_partition_res();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_screeners( pose::Pose & pose )
{
	screeners_.clear();
	if ( rigid_body_modeler_ ) initialize_residue_level_screeners( pose ); // do not keep this commented out!
	initialize_pose_level_screeners( pose );
}

////////////////////////////////////////////////////////////////////////////////
// geometry checks that require base conformation (but not sugar) at moving_res
////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_residue_level_screeners( pose::Pose & pose ) {

	using namespace screener;
	using namespace core::conformation;
	using utility::tools::make_vector1;

	runtime_assert( rigid_body_rotamer_ != 0 );

	screeners_.push_back( new StubApplier( moving_res_base_stub_ ) ); // will pull stub out of the sampler

	screeners_.push_back( new StubDistanceScreener( moving_res_base_stub_, rigid_body_rotamer_->reference_stub(),	max_distance_squared_ ) );

	if ( base_centroid_checker_ ) screeners_.push_back( new BaseCentroidScreener( base_centroid_checker_,
																																								moving_res_base_stub_ ) );

	tag_definition_ = new TagDefinition( pose, screeners_[ screeners_.size() ] );
	screeners_.push_back( tag_definition_ );

	for ( Size n = 1; n <= five_prime_chain_breaks_.size(); n++ ){
		screeners_.push_back( new RNA_ChainClosableGeometryResidueBasedScreener( chain_closable_geometry_checkers_[ n ] ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// clash checks
	ResidueCOP screening_moving_rsd_at_origin = rigid_body_rotamer_->get_residue_at_origin( moving_res_ ).clone();
	if ( VDW_bin_checker_ != 0 ){
		screeners_.push_back( new VDW_BinScreener( VDW_bin_checker_, *virt_sugar_screening_pose_, moving_res_,
																							 screening_moving_rsd_at_origin, moving_res_base_stub_ ) );
	}
	// User-input VDW: Does not work for chain_closure move and is_internal_ move yet, since the checker does not know that
	// moving residue atoms can bond to previous or next residues.
	if ( user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ) {
		screeners_.push_back( new VDW_BinScreener( user_input_VDW_bin_checker_, *virt_sugar_screening_pose_, moving_res_,
																							 screening_moving_rsd_at_origin, moving_res_base_stub_ ) );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// 'Standard screens' that look at entire pose.
//////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_pose_level_screeners( pose::Pose & pose ) {

	using namespace screener;
	using namespace core::conformation;
	using utility::tools::make_vector1;

	//////////////////////////////////////////////////////////////////////////////////////////////
	// finally copy in the base conformation into the pose for RMSD and VDW screens
	// at first, don't copy in dofs for sugar/backbone (i.e., residue alternative), because those copy_dofs take extra computation;
	//  just apply rigid_body transformation.
	screeners_.push_back( new stepwise::screener::SampleApplier( *screening_pose_, false /*apply_residue_alternative_sampler*/ ) );

	NativeRMSD_ScreenerOP native_rmsd_screener;
	if ( get_native_pose() ){
		bool do_screen = ( ( options_->rmsd_screen() > 0.0 ) && !options_->integration_test_mode() ); // gets toggled to true in integration tests.
		native_rmsd_screener = new NativeRMSD_Screener( *get_native_pose(), *screening_pose_,
																										working_parameters_, options_->rmsd_screen(),
																										do_screen );
		screeners_.push_back( native_rmsd_screener );
	}

	// For KIC closure, immediate check that a closed loop solution was actually found.
	if ( kic_modeler_ ) {
		screeners_.push_back( new RNA_ChainClosureScreener( chain_closure_checkers_[ 1 ], *screening_pose_, true /*just do closure check*/ ) );
	}

	// Following may still work, but has not been tested.
	//	if ( options_->combine_long_loop_mode()  && ( cutpoints_closed_.size() == 0 ) ) {
	//		screeners_.push_back( new ResidueContactScreener( *screening_pose_, last_append_res_,  last_prepend_res_, atom_atom_overlap_dist_cutoff_ ) );
	//	}

	if ( !rigid_body_modeler_ && base_centroid_checker_ ){
		bool const force_centroid_interaction = ( rigid_body_modeler_ || options_->force_centroid_interaction()
																							|| ( cutpoints_closed_.size() == 0 ) );
		screeners_.push_back( new BaseCentroidScreener( base_centroid_checker_, screening_pose_, force_centroid_interaction ) );
	}

	if ( VDW_bin_checker_ ){
		screeners_.push_back( new VDW_BinScreener( VDW_bin_checker_, *screening_pose_, moving_res_ ) );
	}
	if ( user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ){
		screeners_.push_back( new VDW_BinScreener( user_input_VDW_bin_checker_, *screening_pose_, moving_res_ ) );
	}

	// This one can actually filters out models compared to the other atr-rep screener below -- since it
	// involves a virtualized sugar, can fail atr screen.
	if ( rigid_body_modeler_ && virt_sugar_atr_rep_screen_ && options_->atr_rep_screen() ) {
		screeners_.push_back( new stepwise::screener::SampleApplier( *virt_sugar_screening_pose_, false /*apply_residue_alternative_sampler*/ ) );
		screeners_.push_back( new RNA_AtrRepScreener( virt_sugar_atr_rep_checker_, *virt_sugar_screening_pose_ ) );
	}

	screeners_.push_back( new stepwise::screener::SampleApplier( *screening_pose_, true /*apply_residue_alternative_sampler*/ ) );
	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) screeners_.push_back( new RNA_ChainClosableGeometryScreener( chain_closable_geometry_checkers_[ n ], screening_pose_ ) );

	for ( Size n = 1; n <= five_prime_chain_breaks_.size(); n++ ) {
		bool strict = rigid_body_modeler_ && cutpoints_closed_.has_value( five_prime_chain_breaks_[n] ) && (cutpoints_closed_.size() < 3);
		screeners_.push_back( new RNA_ChainClosableGeometryScreener( chain_closable_geometry_checkers_[ n ], screening_pose_, strict  /*strict*/ ) );
	}

	RNA_AtrRepScreenerOP atr_rep_screener;
	if ( options_->atr_rep_screen() ) atr_rep_screener = new RNA_AtrRepScreener( atr_rep_checker_, *screening_pose_ );
	screeners_.push_back( atr_rep_screener );

	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) {
		if ( kic_modeler_ && n == 1 ) continue; // if KIC, first one is screened above, actually.
		screeners_.push_back( new RNA_ChainClosureScreener( chain_closure_checkers_[ n ] ) );
	}

	if ( options_->sampler_perform_phosphate_pack() ) screeners_.push_back( new PhosphateScreener( phosphate_sampler_ ) );

	if ( rigid_body_modeler_ && try_sugar_instantiation_ )	screeners_.push_back( new SugarInstantiator( *screening_pose_, moving_res_,
																																																			 o2prime_instantiation_distance_cutoff_ ) );

	if ( options_->sampler_perform_o2prime_pack() )	screeners_.push_back( new O2PrimeScreener( o2prime_packer_ ) );

	/////////////////////////////////////////////////////
	screeners_.push_back( new stepwise::screener::SampleApplier( pose ) );

	if ( !tag_definition_ ){ // may have been defined above in residue level modeler.
		tag_definition_ = new TagDefinition( pose, screeners_[1], options_->sampler_include_torsion_value_in_tag(),
																				 moving_res_, reference_res_, extra_tag_ );
		screeners_.push_back( tag_definition_ );
	}

	if ( !rigid_body_modeler_ && !rebuild_bulge_mode_ &&
			 options_->allow_bulge_at_chainbreak() && moving_partition_res_.size() == 1 &&
			 ( cutpoints_closed_.size() > 0 )  ) {
		screeners_.push_back( new BulgeApplier( atr_rep_checker_, base_centroid_checker_, moving_res_ ) ); // apply bulge at the last minute.
	}

	screeners_.push_back( new PoseSelectionScreener( pose_selection_, pose /*const reference*/, tag_definition_,
																									 options_->verbose(), silent_file_, get_native_pose(), working_parameters_ ) );

	if ( rigid_body_modeler_ ) {
		screeners_.push_back( new BaseBinMapUpdater( base_bin_map_ ) );
		// As long as we're not instantiating any sugars in the final pose, break once we find any valid sugar rotamer [really?]
		if ( chain_closure_checkers_.size() == 0 && !try_sugar_instantiation_ ) {
			screeners_.push_back( new FastForwardToNextRigidBody ); // generalize to non-rigid-body rotamer case.
		} else {
			screeners_.push_back( new FastForwardToNextResidueAlternative( moving_res_ ) ); // generalize -- should fast forward past all virtual residue alternatives.
		}
	}

	if ( options_->integration_test_mode() ) {
		screeners_.insert( screeners_.begin() /*right at beginning!*/,
											 new IntegrationTestBreaker( atr_rep_screener, screeners_[ screeners_.size() ], native_rmsd_screener ) );
	}

}


/////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_ConnectionSampler::initialize_poses_and_checkers( pose::Pose & pose  ){

	bool ready_for_more = presample_virtual_sugars( pose ); //defines alternatives for residues at chainbreaks with a virtual sugar.
	if ( !ready_for_more ) return false;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// initialization of variants for actual pose.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( options_->sampler_perform_phosphate_pack() ){
		phosphate_sampler_ = new phosphate::MultiPhosphateSampler( pose, moving_res_ );
		runtime_assert(  moving_partition_res_ == phosphate_sampler_->moving_partition_res()/*determined inside*/ );
	}

	if ( options_->sampler_perform_o2prime_pack() ) {
		remove_virtual_O2Prime_hydrogen( pose );
		o2prime_packer_ = new o2prime::O2PrimePacker( pose, scorefxn_, moving_partition_res_ );
	} else {
		add_virtual_O2Prime_hydrogen( pose );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// RNA Base Centroid stuff -- could soon generalize to find protein partners too. That would be cool.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	base_centroid_checker_ = new RNA_BaseCentroidChecker ( pose, working_parameters_,
																												 options_->tether_jump());
	base_centroid_checker_->set_floating_base( working_parameters_->floating_base() &&
																						 working_parameters_->working_moving_partition_res().size() == 1  );
	base_centroid_checker_->set_allow_base_pair_only_screen( options_->allow_base_pair_only_centroid_screen() );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// info on chain-breaks -- note update to constraints [?] on pose.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	figure_out_moving_chain_breaks( pose, moving_partition_res_,
																	cutpoints_closed_,
																	five_prime_chain_breaks_, three_prime_chain_breaks_, chain_break_gap_sizes_ );
	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) reinstantiate_backbone_and_add_constraint_at_moving_res( pose, cutpoints_closed_[ n ] );
	kic_modeler_ = ( !rigid_body_modeler_ && options_->kic_modeler_if_relevant() && ( cutpoints_closed_.size() > 0 ) );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// sets up residue alternatives for moving_res. Later move this *out* of ConnectionSampler along with other sugar stuff?
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( rigid_body_modeler_ ) initialize_moving_residue_pose_list( pose );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// set up screening pose -- do not change pose itself.
	// get rid of stuff that will be CCD-closed or packed (2'-OH, terminal phosphates) at last stages.
	// The idea is that if the screening pose fails basic stereochemistry checks, then we don't
	// have to carry out expensive CCD closure or packing.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// sets up cutpoints_closed_, five_prime_chain_break_res_, three_prime_chain_break_res_ and chainbreak_gaps_
	screening_pose_ = pose.clone();
	add_virtual_O2Prime_hydrogen( *screening_pose_ );
	phosphate::remove_terminal_phosphates( *screening_pose_ );
	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) {
		add_variant_type_to_pose_residue( *screening_pose_, chemical::VIRTUAL_PHOSPHATE, cutpoints_closed_[n] + 1 ); // PS May 31, 2010 -- updated to all cutpoints by rhiju, feb. 2014
	}
	// VDW bin checker can take a while to set up... becomes rate-limiting in random.
	if ( !options_->choose_random() ){
		VDW_bin_checker_ = new checker::RNA_VDW_BinChecker();
		VDW_bin_checker_->setup_using_working_pose( *screening_pose_, working_parameters_ );
	}
	if ( !user_input_VDW_bin_checker_ /* could be externally defined for speed */ ) {
		user_input_VDW_bin_checker_ = new RNA_VDW_BinChecker();
		if ( options_->VDW_rep_screen_info().size() > 0 ){
			options_->setup_options_for_VDW_bin_checker( user_input_VDW_bin_checker_ );
			user_input_VDW_bin_checker_->setup_using_user_input_VDW_pose( options_->VDW_rep_screen_info(),
																																		pose, working_parameters_ );
		}
	}


	virt_sugar_screening_pose_ = screening_pose_->clone(); //Hard copy. Used for trying out sugar at moving residue.
	// virtual sugars even at residues that have instantiated sugars -- we can quickly screen this pose,
	// and it provides the appropriate baseline atr/rep for checking contacts and clashes.
	for ( Size n = 1; n <= residue_alternative_sets_.size(); n++ ){
		pose::add_variant_type_to_pose_residue( *virt_sugar_screening_pose_, chemical::VIRTUAL_RIBOSE,
																							residue_alternative_sets_[ n ].representative_seqpos() );
	}
	// following is to check atr/rep even on sugars that will remain virtualized.
	bool const use_loose_rep_cutoff = ( kic_modeler_ || moving_partition_res_.size() > 1 );
	atr_rep_checker_ = new checker::RNA_AtrRepChecker( *screening_pose_, working_parameters_, use_loose_rep_cutoff );
	virt_sugar_atr_rep_checker_ = new checker::RNA_AtrRepChecker( *virt_sugar_screening_pose_, working_parameters_, use_loose_rep_cutoff );

	// we will be checking clashes of even virtual sugars compared to no-sugar baseline.
	for ( Size n = 1; n <= residue_alternative_sets_.size(); n++ ){
		pose::remove_variant_type_from_pose_residue( *screening_pose_, chemical::VIRTUAL_RIBOSE,
																								 residue_alternative_sets_[ n ].representative_seqpos() );
	}

	for ( Size n = 1; n <= five_prime_chain_breaks_.size(); n++ ) {
		chain_closable_geometry_checkers_.push_back( new checker::RNA_ChainClosableGeometryChecker( five_prime_chain_breaks_[n], three_prime_chain_breaks_[n], chain_break_gap_sizes_[n] ) );
	}

	screening_pose_->remove_constraints(); // chain closure actually assumes no constraints in pose.
	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) {
		chain_closure_checkers_.push_back( new checker::RNA_ChainClosureChecker( *screening_pose_, cutpoints_closed_[n] ) );
		chain_closure_checkers_[n]->set_reinitialize_CCD_torsions( options_->reinitialize_CCD_torsions() );
	}

	Size num_pose_kept_to_use =  get_num_pose_kept();
	pose_selection_ = new StepWiseRNA_PoseSelection( working_parameters_, scorefxn_ );
	pose_selection_->set_num_pose_kept( num_pose_kept_to_use );
	if ( options_->cluster_rmsd() > 0.0 ) pose_selection_->set_cluster_rmsd( options_->cluster_rmsd() );
	pose_selection_->set_PBP_clustering_at_chain_closure( options_->PBP_clustering_at_chain_closure() );
	pose_selection_->set_distinguish_pucker( options_->distinguish_pucker() );
	pose_selection_->set_pose_list( pose_list_ );

	return true;
}

/////////////////////////////////////////////////////////////////////////////////
Size
StepWiseRNA_ConnectionSampler::get_max_ntries() {
	Size max_ntries( 0 );
	if ( rigid_body_modeler_ ) {
		max_ntries = std::max( 100000, 1000 * int( options_->num_random_samples() ) );
		if ( options_->rmsd_screen() && !options_->integration_test_mode() ) max_ntries *= 10;
	} else {
		max_ntries = std::max( 10000, 100 * int( options_->num_random_samples() ) );
		if ( chain_closure_checkers_.size() > 0 ) max_ntries *= 10;
		if ( kic_modeler_ ) max_ntries = 5 * options_->num_random_samples(); // some chains just aren't closable.
	}
	return max_ntries;
}


/////////////////////////////////////////////////////////////////////////////////
Size
StepWiseRNA_ConnectionSampler::get_num_pose_kept(){
	Size num_pose_kept( 108 );
	if ( options_->sampler_num_pose_kept() > 0 ) num_pose_kept = options_->sampler_num_pose_kept();
	//	if ( build_pose_from_scratch_ )	 num_pose_kept = 36* num_pose_kept;
	if ( !rigid_body_modeler_ && working_parameters_ && working_parameters_->sample_both_sugar_base_rotamer() ) num_pose_kept *= 12;
	if ( rigid_body_modeler_ && base_centroid_checker_->allow_base_pair_only_screen() ) num_pose_kept *= 4;
	return num_pose_kept;
}

/////////////////////////////////////////////////////////////////////////////////////
Size
StepWiseRNA_ConnectionSampler::which_residue_alternative_set_is_moving_residue() const {
	Size which_set_is_moving_res( 0 );
	for ( Size n = 1; n <= residue_alternative_sets_.size(); n++ ){
		if ( residue_alternative_sets_[n].representative_seqpos() == moving_res_ ) which_set_is_moving_res = n;
	}
	return which_set_is_moving_res;
}

/////////////////////////////////////////////////////////////////////////////////////
// sets up a 'bare-bones' sugar modeling object with moving_residue information.
// in the (near) future change SugarModeling object to something more general with
// a pose_list, res_map, and representative residue -- ResidueAlternativeSet [?]
void
StepWiseRNA_ConnectionSampler::initialize_moving_residue_pose_list( pose::Pose const & pose ){
	if ( which_residue_alternative_set_is_moving_residue() > 0 ) return; // already initialized.
	utility::vector1< pose::PoseOP > pose_list;
	if ( rigid_body_modeler_ && moving_partition_res_.size() == 1 ){ // single floating base [classic]
		pose_list = setup_pose_with_moving_residue_alternative_list( pose, moving_res_, options_->extra_chi(), options_->use_phenix_geo() );
	} else {
		pose_list = utility::tools::make_vector1( pose.clone() ); // no alternatives.
	}
	sampler::copy_dofs::ResidueAlternativeSet residue_alternative_set( pose_list, moving_res_ );
	residue_alternative_sets_.push_back( residue_alternative_set );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::reinstantiate_backbone_and_add_constraint_at_moving_res(	pose::Pose & pose, Size const & five_prime_chain_break_res )
{ //harmonic angle and distance constraints are used ONLY by chainbreak_screening
	pose::remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_RIBOSE, moving_res_ ); //May 31, 2010
	pose::remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_O2PRIME_HYDROGEN, moving_res_ );
	if ( moving_res_ == ( five_prime_chain_break_res + 1 ) ){
		pose::remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_PHOSPHATE, moving_res_ ); //this virtual_phosphate was added to pose at the beginning of this function.
	}
	add_harmonic_chain_break_constraint( pose, five_prime_chain_break_res );
	runtime_assert( !pose.residue( five_prime_chain_break_res   ).has_variant_type( chemical::VIRTUAL_RIBOSE ) );
	runtime_assert( !pose.residue( five_prime_chain_break_res+1 ).has_variant_type( chemical::VIRTUAL_RIBOSE ) );
}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_euler_angle_grid_parameters(){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// definition of euler angle grid search parameters
	// Following are set by #define in StepWiseRNA_FloatingBase_Samper_util_hh
	sampler::rigid_body::RigidBodyStepWiseSamplerValueRange & value_range = rigid_body_rotamer_->value_range();
	value_range.set_euler_angle_bin_size( STANDARD_EULER_ANGLE_BIN_SIZE );
	value_range.set_euler_z_bin_size( STANDARD_EULER_Z_BIN_SIZE );
	value_range.set_centroid_bin_size( STANDARD_CENTROID_BIN_SIZE );
	if ( options_->integration_test_mode() ){ // use coarser search for speed
		value_range.set_euler_angle_bin_size( STANDARD_EULER_ANGLE_BIN_SIZE * 4 );
		value_range.set_centroid_bin_size( STANDARD_CENTROID_BIN_SIZE * 4);
	}
}

//////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_xyz_grid_parameters(){
	Distance max_distance = options_->sampler_max_centroid_distance(); // if unspecified (0.0), will be replaced
	if ( options_->tether_jump() && max_distance == 0.0 ) max_distance = 8.0;
	int centroid_bin_min, centroid_bin_max;
	initialize_xyz_parameters( max_distance, max_distance_squared_,
														 centroid_bin_min, centroid_bin_max,
														 get_moving_rsd_list(), truly_floating_base() );
	sampler::rigid_body::RigidBodyStepWiseSamplerValueRange & value_range = rigid_body_rotamer_->value_range();
	value_range.set_max_distance( max_distance );
	value_range.set_centroid_bin_min( centroid_bin_min );
	value_range.set_centroid_bin_max( centroid_bin_max );
	base_bin_map_.clear();
}

//////////////////////////////////////////////////////////////////////
// used in setting max_distance -- can use a tighter tether if moving res is covalently connected to the reference (anchor) residue.
Size  // Size instead of bool for historical reasons.
StepWiseRNA_ConnectionSampler::truly_floating_base(){
	if ( cutpoints_closed_.has_value( moving_res_    ) && reference_res_ == moving_res_+1 ) return 0;
	if ( cutpoints_closed_.has_value( reference_res_ ) && reference_res_ == moving_res_-1 ) return 0;
	return 1;
}

//////////////////////////////////////////////////////////////////////
// kind of here for historical reasons -- set max distance based on
// list of all the options for moving_residue. This will probably
// be unnecessary soon when we shift to a mode where max_distance
// is a fixed number.
utility::vector1< core::conformation::ResidueOP >
StepWiseRNA_ConnectionSampler::get_moving_rsd_list() const {
	runtime_assert( rigid_body_modeler_ );
	Size const n = which_residue_alternative_set_is_moving_residue();
	runtime_assert( n > 0 );
	utility::vector1< core::conformation::ResidueOP > moving_rsd_list;
	sampler::copy_dofs::ResidueAlternativeSet const & moving_res_modeling = residue_alternative_sets_[ n ];
	for ( Size n = 1; n <= moving_res_modeling.pose_list().size(); n++ ){
		moving_rsd_list.push_back( moving_res_modeling.pose( n )->residue( moving_res_ ).clone() );
	}
	return moving_rsd_list;
}

//////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_sampler(){
	if ( rigid_body_modeler_ ){
		initialize_euler_angle_grid_parameters();
		initialize_xyz_grid_parameters();
		initialize_full_rigid_body_sampler();
	} else {
		sampler_ = get_full_bond_sampler();
	}
}

//////////////////////////////////////////////////////////////////////
sampler::StepWiseSamplerBaseOP
StepWiseRNA_ConnectionSampler::get_full_bond_sampler(){
	using namespace sampler;
	StepWiseSamplerBaseOP sampler_ = sampler::rna::setup_sampler( *screening_pose_, options_,
																																								working_parameters_, build_pose_from_scratch_,
																																								kic_modeler_, (cutpoints_closed_.size() > 0) );
	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer = get_rsd_alternatives_rotamer();
	if ( rsd_alternatives_rotamer == 0 ) return sampler_;

	StepWiseSamplerCombOP sampler = new StepWiseSamplerComb;
	sampler->add_external_loop_rotamer( rsd_alternatives_rotamer );
	sampler->add_external_loop_rotamer( sampler_ );
	sampler->set_random( options_->choose_random() );
	sampler->init();
	return sampler;
}

//////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::initialize_full_rigid_body_sampler(){

	using namespace sampler::rigid_body;
	using namespace sampler::copy_dofs;

	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer = get_rsd_alternatives_rotamer();
	sampler_ = new RigidBodyStepWiseSamplerWithResidueAlternatives( rsd_alternatives_rotamer, rigid_body_rotamer_ );
	sampler_->set_random( options_->choose_random() );
	sampler_->init();
}

//////////////////////////////////////////////////////////////////////
sampler::copy_dofs::ResidueAlternativeStepWiseSamplerCombOP
StepWiseRNA_ConnectionSampler::get_rsd_alternatives_rotamer(){

	if ( residue_alternative_sets_.size() == 0 ) return 0;
	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer =	new ResidueAlternativeStepWiseSamplerComb();
	// note that following will include moving_res_ for sampler, as well as any other chunks that might move...
	for ( Size n = 1; n <= residue_alternative_sets_.size(); n++ ){
		ResidueAlternativeSet const & residue_alternative_set = residue_alternative_sets_[n];
		ResidueAlternativeStepWiseSamplerOP rsd_alt_rotamer;
		if ( rigid_body_rotamer_ != 0 ){
			rsd_alt_rotamer = new ResidueAlternativeStepWiseSampler( residue_alternative_set,
																																									 *rigid_body_rotamer_->pose_at_origin() /*take representative residues from this pose after applying copy_dofs*/);
		} else {
			rsd_alt_rotamer = new ResidueAlternativeStepWiseSampler( residue_alternative_set );
		}
		rsd_alternatives_rotamer->add_residue_alternative_rotamer( rsd_alt_rotamer );
	}
	return rsd_alternatives_rotamer;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::set_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

void
StepWiseRNA_ConnectionSampler::add_residue_alternative_set( sampler::copy_dofs::ResidueAlternativeSet const & residue_alternative_set ){
	residue_alternative_sets_.push_back( residue_alternative_set );
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseRNA_ConnectionSampler::get_pose_list(){
	return pose_list_;
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::set_pose_list( utility::vector1< pose::PoseOP > &	pose_list ){
	pose_list_ = pose_list;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ConnectionSampler::set_options( options::StepWiseModelerOptionsCOP options ){
	options_ = options;
}


/////////////////////////////////////////////////////////////////////////////////////////
// This function does nothing of consequence. Used as a consistency check with old-style
// working_parameters-inputted run.
void
StepWiseRNA_ConnectionSampler::check_working_parameters( pose::Pose const & pose ){

	using namespace core::pose::full_model_info;

	if ( !working_parameters_ ) return;

	runtime_assert( moving_res_ ==  working_parameters_->working_moving_res() );
	runtime_assert( moving_partition_res_ == working_parameters_->working_moving_partition_res() );

	bool const is_prepend_ = working_parameters_->is_prepend();
	Size const five_prime_chain_break_res_( working_parameters_->five_prime_chain_break_res() );
	if ( five_prime_chain_break_res_ ) runtime_assert( five_prime_chain_breaks_.has_value( five_prime_chain_break_res_ ) );

	if ( rigid_body_modeler_ ){
		if ( working_parameters_->floating_base_anchor_res() > 0 ) runtime_assert( reference_res_ == working_parameters_->full_to_sub( working_parameters_->floating_base_anchor_res() ) );
		if ( moving_partition_res_.size() == 1 ) {
			runtime_assert( reference_res_ == working_parameters_->working_reference_res() );
			bool const same_chain = ( get_chain_for_resnum( moving_res_, pose ) ==
																get_chain_for_resnum( reference_res_, pose ) );
			if ( same_chain ) {
				Size const floating_base_five_prime_chain_break_ ( ( is_prepend_ ) ? moving_res_   : reference_res_ );
				runtime_assert( five_prime_chain_breaks_.has_value( floating_base_five_prime_chain_break_ ) );
			}
		}
	}

	//	Size const gap_size_ = working_parameters_->gap_size(); /* If this is zero or one, need to screen or closable chain break */
	Size const gap_size_to_anchor_ = working_parameters_->gap_size_to_anchor();
	if ( rigid_body_modeler_ && gap_size_to_anchor_ == 0 )	runtime_assert( pose.fold_tree().is_cutpoint( moving_res_ - 1 ) );

	bool const 	is_dinucleotide_( gap_size_to_anchor_ == 1 );
	runtime_assert ( !is_dinucleotide_ || !working_parameters_->is_internal() );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_ConnectionSampler::presample_virtual_sugars( pose::Pose & pose ){
	using namespace modeler::rna::sugar;
	StepWiseRNA_VirtualSugarJustInTimeInstantiatorOP virtual_sugar_just_in_time_instantiator =
		instantiate_any_virtual_sugars( pose, working_parameters_, scorefxn_, options_ );
	if ( !virtual_sugar_just_in_time_instantiator->success() )	return false;
	virtual_sugar_just_in_time_instantiator->instantiate_sugars_at_cutpoint_closed( pose );

	// in backwards order to match some old runs.
	for ( Size n = virtual_sugar_just_in_time_instantiator->num_sets(); n >= 1; n-- ) {
	 	add_residue_alternative_set( virtual_sugar_just_in_time_instantiator->residue_alternative_set( n ) );
	}
	return true;
}


} //connection
} //rna
} //modeler
} //legacy
} //stepwise
} //protocols
