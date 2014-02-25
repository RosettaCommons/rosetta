// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/legacy/rigid_body/StepWiseRNA_RigidBodyConnectionSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/sampling/rna/legacy/rigid_body/StepWiseRNA_RigidBodyConnectionSampler.hh>
#include <protocols/stepwise/sampling/rna/rigid_body/StepWiseRNA_FloatingBaseSamplerUtil.hh>
#include <protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosureChecker.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/sampling/rna/phosphate/PhosphateUtil.hh>
#include <protocols/stepwise/screener/AtrRepScreener.hh>
#include <protocols/stepwise/screener/BaseBinMapUpdater.hh>
#include <protocols/stepwise/screener/BaseCentroidScreener.hh>
#include <protocols/stepwise/screener/BulgeApplier.hh>
#include <protocols/stepwise/screener/CentroidDistanceScreener.hh>
#include <protocols/stepwise/screener/ChainClosableGeometryScreener.hh>
#include <protocols/stepwise/screener/ChainClosableGeometryResidueBasedScreener.hh>
#include <protocols/stepwise/screener/ChainClosureScreener.hh>
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
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamer.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerComb.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

using namespace core;
using ObjexxFCL::string_of;
using ObjexxFCL::lead_zero_string_of;

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_RigidBodyConnectionSampler" ) ;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// "Floating base", a.k.a. dinucleotide or "skip bulge" sampling, developed by Parin Sripakdeevong.
//
// Intead of enumerating torsions for a residue directly attached to the fixed pose (as in StepWiseRNA_ClassicResidueSampler),
//  this code samples rigid body degrees of freedom for a residue's base, allowing for one or more (or even zero!)
//  nucleotides intervening between the residue and the pose.
//
// The main use case in SWA is 'dinucleotide' (single intervening residue that
//  is virtualized). Example
//
//        Anchor residue                                 Moving residue                 Distal residue
//                       (virtual) (virtual)     (virtual*)   (virtual*)        (virtual*)
//     5'   --Sugar -- [ Phosphate -- Sugar ] -- Phosphate -- Sugar -- ... -- Phosphate -- Sugar -- Phosphate -- 3'
//              |      [               |    ]                   |                           |
//          Reference  [             Bulge  ]                Floating                      Distal
//             Base                (virtual)                   Base!                        Base
//                                            |                          |
//                     If no bulge, close_chain_to_anchor             If no residues to distal, close_chain_to_distal
//
//              | <--------  gap_size_to_anchor + 1  ----------->|<------- gap_size + 1 ----->|
//
//
//  * usually virtual -- the exceptions are if connected through chain closure to anchor and/or distal residue.
//
// As with StepWiseRNA_ClassicResidueSampler, this Sampler uses several poses and checkers to reduce the number of variant changes &
// energy recomputations -- faster code but somewhat complicated.
//
// Note that sugar of the floating base will be virtualized unless it needs to be instantiated to close chain to
//  anchor or distal residues. Nevertheless, code below carried out a geometric 'sanity check' that involves temporarily
//  instantiating the sugar.
//
// On the way to dinucleotide sampling, we often start with poses in which a virtual sugar in the attachment point in the fixed
//  pose was instantiated. One option would be to run this dinucleotide sampler separately for each of those instantiation cases.
//  That was the choice made in ClassicResidueSampler -- and the code is simpler there. But here, that information is passed into
//  the sampler in 'residue_alternative_set'; some of the first screens are based on base centroid distances and are independent of
//  sugar conformation, and so some computation that would be redundant across anchor sugar conformations can be avoided.
//
//   -- Rhiju, 2013, 2014.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {
namespace rigid_body {

//Constructor
StepWiseRNA_RigidBodyConnectionSampler::StepWiseRNA_RigidBodyConnectionSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
	job_parameters_( job_parameters ),
	moving_res_(  job_parameters_->working_moving_res() ), // Might not correspond to user input.
	moving_suite_(  job_parameters_->working_moving_suite() ), // dofs betweeen this value and value+1 actually move.
	is_prepend_(  job_parameters_->is_prepend() ),
	is_internal_(  job_parameters_->is_internal() ), // no cutpoints before or after moving_res.
	gap_size_( job_parameters_->gap_size() ), /* If this is zero or one, need to screen or closable chain break */
	gap_size_to_anchor_( job_parameters_->gap_size_to_anchor() ),
	five_prime_chain_break_res_( job_parameters_->five_prime_chain_break_res() ),
	chain_break_reference_res_( ( is_prepend_ ) ? five_prime_chain_break_res_ : five_prime_chain_break_res_ + 1 ),
	reference_res_( job_parameters_->working_reference_res() ), //the last static_residues that this attach to the moving residues
	floating_base_five_prime_chain_break_ ( ( is_prepend_ ) ? moving_res_   : reference_res_ ),
	floating_base_three_prime_chain_break_( ( is_prepend_ ) ? reference_res_: moving_res_ ),
	is_dinucleotide_( gap_size_to_anchor_ == 1 ),
	close_chain_to_distal_( gap_size_ == 0 ),
	close_chain_to_anchor_( gap_size_to_anchor_ == 0 ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ), // can be replaced from the outside
	silent_file_( "silent_file.txt" ),
	native_rmsd_screen_( false ), // updated below
	num_pose_kept_to_use_( 0 ),   // updated below
	euler_angle_bin_min_( 0 ),    // updated below
	euler_angle_bin_max_( 0 ),    // updated below
	euler_z_bin_min_( 0 ),        // updated below
	euler_z_bin_max_( 0 ),        // updated below
	centroid_bin_min_( 0 ),       // updated below
	centroid_bin_max_( 0 ),       // updated below
	max_distance_( 0.0 ), // updated below
	max_distance_squared_( 0.0 ), // updated below
	try_sugar_instantiation_( false ),
	o2prime_instantiation_distance_cutoff_( 6.0 ),
	extra_tag_( "" ),
	rigid_body_sampling_( true ), // will not be true if we unify with suite sampling
	residue_level_screening_( true ),
	full_pose_level_screening_( true ) // will make unification with suite sampling simpler.
{
	set_native_pose( job_parameters_->working_native_pose() );
	runtime_assert ( !is_dinucleotide_ || !is_internal_ );

	// remove following soon.
	TR << "GAP_SIZE_TO_ANCHOR " << gap_size_to_anchor_ << "  REFERENCE RES        " << reference_res_ << "  MOVING_RES " << moving_res_ << std::endl;
	TR << "GAP_SIZE_TO_DISTAL " << gap_size_           << "  FIVE' CHAINBREAK RES " << five_prime_chain_break_res_ << "  MOVING_RES " << moving_res_ << std::endl;
}

//Destructor
StepWiseRNA_RigidBodyConnectionSampler::~StepWiseRNA_RigidBodyConnectionSampler()
{}

/////////////////////
std::string
StepWiseRNA_RigidBodyConnectionSampler::get_name() const {
	return "StepWiseRNA_RigidBodyConnectionSampler";
}

///////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::apply( core::pose::Pose & pose ){
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace protocols::rotamer_sampler::rigid_body;

	output_title_text( "Enter StepWiseRNA_RigidBodyConnectionSampler::floating_base_sampling", TR.Debug );

	clock_t const time_start( clock() );

	initialize_poses_and_stubs_and_checkers( pose );
	initialize_euler_angle_grid_parameters();
	initialize_xyz_grid_parameters();
	initialize_rigid_body_sampler( pose );
	initialize_screeners( pose );

	StepWiseSampleAndScreen sample_and_screen( sampler_, screeners_ );
	sample_and_screen.set_max_ntries( max_ntries_ );

	// Do it!
	TR << "KICKING OFF SAMPLE AND SCREEN!!!!!!!!!! " << std::endl;
	sample_and_screen.run();
	sample_and_screen.output_counts();
	sample_and_screen.output_info_on_random_trials();

	pose_selection_->finalize();
	pose_list_ = pose_selection_->pose_list();
	if ( options_->verbose() ) analyze_base_bin_map( base_bin_map_, "test/" );

	TR.Debug << "floating base sampling time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_screeners( pose::Pose & pose )
{
	screeners_.clear();
	if ( residue_level_screening_ ) initialize_residue_level_screeners( pose ); // do not keep this commented out!
	initialize_pose_level_screeners( pose );
}

////////////////////////////////////////////////////////////////////////////////
// geometry checks that require base conformation (but not sugar) at moving_res
////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_residue_level_screeners( pose::Pose & pose ) {

	using namespace screener;
	using namespace core::conformation;
	using utility::tools::make_vector1;

	screeners_.push_back( new StubApplier( moving_res_base_stub_ ) ); // will pull stub out of the sampler

	screeners_.push_back( new StubDistanceScreener( moving_res_base_stub_, reference_stub_,	max_distance_squared_ ) );

	if ( base_centroid_checker_ ) screeners_.push_back( new BaseCentroidScreener( base_centroid_checker_,
																																								moving_res_base_stub_ ) );

	tag_definition_ = new TagDefinition( pose, screeners_[ screeners_.size() ] );
	screeners_.push_back( tag_definition_ );

	for ( Size n = 1; n <= five_prime_chain_breaks_.size(); n++ ) screeners_.push_back( new ChainClosableGeometryResidueBasedScreener( chain_closable_geometry_checkers_[ n ] ) );

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// clash checks
	// Some trickiness with clashes to 3'-phosphate in case of chain closure -- parin has worked out these cases, fortunately.
	ResidueCOP screening_moving_rsd_at_origin = sampler_->get_residue( moving_res_ );
	screeners_.push_back( new VDW_BinScreener( VDW_bin_checker_, *screening_pose_, moving_res_,
																						 screening_moving_rsd_at_origin, moving_res_base_stub_ ) );

	// User-input VDW: Does not work for chain_closure move and is_internal_ move yet, since the checker does not know that
	// moving residue atoms can bond to previous or next residues.
	if ( ( user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ) && ( !close_chain_to_distal_ ) && ( !is_internal_ )  ){
		screeners_.push_back( new VDW_BinScreener( user_input_VDW_bin_checker_, *screening_pose_, moving_res_,
																							 screening_moving_rsd_at_origin, moving_res_base_stub_ ) );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// 'Standard screens' that look at entire pose.
//////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_pose_level_screeners( pose::Pose & pose ) {

	using namespace screener;
	using namespace core::conformation;
	using utility::tools::make_vector1;

	//////////////////////////////////////////////////////////////////////////////////////////////
	// finally copy in the base conformation into the pose for RMSD and VDW screens
	screeners_.push_back( new SampleApplier( *screening_pose_, false /*apply_residue_alternative_sampler*/ ) );

	if ( rigid_body_sampling_ ) screeners_.push_back( new CentroidDistanceScreener( *screening_pose_, moving_res_, reference_stub_.v,	max_distance_squared_ ) );

	NativeRMSD_ScreenerOP native_rmsd_screener = new NativeRMSD_Screener( *get_native_pose(), *screening_pose_,
																																				job_parameters_, options_->sampler_native_screen_rmsd_cutoff(),
																																				native_rmsd_screen_ /*do_screen*/ );
	screeners_.push_back( native_rmsd_screener );

	if ( full_pose_level_screening_ & base_centroid_checker_ ){
		screeners_.push_back( new BaseCentroidScreener( base_centroid_checker_, screening_pose_ ) );
	}

	screeners_.push_back( new VDW_BinScreener( VDW_bin_checker_, *screening_pose_, moving_res_ ) );
	if ( ( user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ) && ( cutpoints_closed_.size() == 0 ) && ( !is_internal_ ) ){
		screeners_.push_back( new VDW_BinScreener( user_input_VDW_bin_checker_, *screening_pose_, moving_res_ ) );
	}

	AtrRepScreenerOP atr_rep_screener;
	if ( options_->VDW_atr_rep_screen() ) atr_rep_screener = new AtrRepScreener( atr_rep_checker_, *screening_pose_ );
	screeners_.push_back( atr_rep_screener );

	/////////////////////////////////////////////////////
	screeners_.push_back( new SampleApplier( *sugar_screening_pose_ ) );

	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) screeners_.push_back( new ChainClosableGeometryScreener( chain_closable_geometry_checkers_[ n ], sugar_screening_pose_ ) );

	if ( options_->VDW_atr_rep_screen() ) screeners_.push_back( new AtrRepScreener( atr_rep_checker_with_instantiated_sugar_, *sugar_screening_pose_ ) );

	for ( Size n = 1; n <= five_prime_chain_breaks_.size(); n++ ) {
		screeners_.push_back( new ChainClosableGeometryScreener( chain_closable_geometry_checkers_[ n ], sugar_screening_pose_,
																										 cutpoints_closed_.has_value( five_prime_chain_breaks_[n] ) /*strict*/ ) );
	}

	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) screeners_.push_back( new ChainClosureScreener( chain_closure_checkers_[ n ] ) );

	if ( options_->sampler_perform_phosphate_pack() ) screeners_.push_back( new PhosphateScreener( phosphate_sampler_ ) );

	if ( rigid_body_sampling_ && try_sugar_instantiation_ )	screeners_.push_back( new SugarInstantiator( *sugar_screening_pose_, moving_res_,
																																																			 o2prime_instantiation_distance_cutoff_ ) );

	if ( options_->sampler_perform_o2prime_pack() )	screeners_.push_back( new O2PrimeScreener( o2prime_packer_ ) );

	/////////////////////////////////////////////////////
	screeners_.push_back( new SampleApplier( pose ) );

	if ( !tag_definition_ ){ // may have been defined above in residue level sampling.
		tag_definition_ = new TagDefinition( pose, screeners_[1], options_->sampler_include_torsion_value_in_tag(),
																				 moving_res_, is_prepend_, extra_tag_ );
		screeners_.push_back( tag_definition_ );
	}

	screeners_.push_back( new PoseSelectionScreener( pose_selection_, pose /*const reference*/, tag_definition_,
																									 options_->verbose(), silent_file_, get_native_pose(), job_parameters_ ) );

	if ( rigid_body_sampling_ ) {
		screeners_.push_back( new BaseBinMapUpdater( base_bin_map_ ) );
		// As long as we're not instantiating any sugars in the final pose, break once we find any valid sugar rotamer [really?]
		if ( chain_closure_checkers_.size() == 0 && !try_sugar_instantiation_ ) {
			screeners_.push_back( new FastForwardToNextRigidBody );
		} else {
			screeners_.push_back( new FastForwardToNextResidueAlternative( moving_res_ ) );
		}
	}

	if ( options_->integration_test_mode() ) {
		screeners_.insert( screeners_.begin() /*right at beginning!*/,
											 new IntegrationTestBreaker( atr_rep_screener, screeners_[ screeners_.size() ], native_rmsd_screener ) );
	}


}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_poses_and_stubs_and_checkers( pose::Pose & pose  ){

	// cleanup variants for pose -- remove terminal phosphates on floating base
	// it would be more general instead to split and 'prepack' phosphates with MultiPhosphateChecker
	// that function would also hold for ClassicResidueSampler.
	remove_variant_type_from_pose_residue( pose, "THREE_PRIME_PHOSPHATE", moving_res_ );
	remove_variant_type_from_pose_residue( pose, "FIVE_PRIME_PHOSPHATE", moving_res_ );
	if ( !pose.residue( moving_res_ ).has_variant_type( "CUTPOINT_UPPER" ) ) {
		add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", moving_res_ );
	}
	if ( !pose.residue( moving_res_ ).has_variant_type( "CUTPOINT_UPPER" ) &&
			 !pose.residue( moving_res_ ).has_variant_type( "CUTPOINT_LOWER" ) ){
		add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res_ );
	}

	// sets up cutpoints_closed_, five_prime_chain_break_res_, three_prime_chain_break_res_ and chainbreak_gaps_
	figure_out_moving_chain_breaks( pose, job_parameters_->working_moving_partition_pos(),
																	cutpoints_closed_,
																	five_prime_chain_breaks_, three_prime_chain_breaks_, chain_break_gap_sizes_ );

	//move out to job_parameters assertion check.
	if ( gap_size_to_anchor_ == 0 )	runtime_assert( pose.fold_tree().is_cutpoint( moving_res_ - 1 ) );

	// hydrogens will be reinstantiated later. perhaps this should be the screening_pose?
	Pose pose_with_virtual_O2prime_hydrogen = pose;
	add_virtual_O2Prime_hydrogen( pose_with_virtual_O2prime_hydrogen );

	reference_stub_ = get_reference_stub( pose_with_virtual_O2prime_hydrogen );

	////////////////////////////////////////Checkers///////////////////////////////////////////////////
	VDW_bin_checker_ = new checker::RNA_VDW_BinChecker();
	// perhaps this should use the screening_pose_, which has terminal phosphates virtualized too.
	VDW_bin_checker_->setup_using_working_pose( pose_with_virtual_O2prime_hydrogen, job_parameters_ );
	if ( user_input_VDW_bin_checker_ ) user_input_VDW_bin_checker_->reference_xyz_consistency_check( reference_stub_.v );

	screening_pose_ = pose_with_virtual_O2prime_hydrogen.clone(); //Hard copy, used for clash checking
	phosphate::remove_terminal_phosphates( *screening_pose_ );
	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) add_variant_type_to_pose_residue( *screening_pose_, "VIRTUAL_PHOSPHATE", cutpoints_closed_[n] + 1 ); // PS May 31, 2010 -- updated to all cutpoints by rhiju, feb. 2014

	sugar_screening_pose_ = screening_pose_->clone(); //Hard copy. Used for trying out sugar at moving residue.
	pose::remove_variant_type_from_pose_residue( *sugar_screening_pose_, "VIRTUAL_RIBOSE", moving_res_ );

	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) {
		reinstantiate_backbone_and_add_constraint_at_moving_res( pose, cutpoints_closed_[ n ] );
	}


	if ( options_->sampler_perform_o2prime_pack() ) {
		remove_virtual_O2Prime_hydrogen( pose );
		o2prime_packer_ = new o2prime::O2PrimePacker( pose, scorefxn_, job_parameters_->working_moving_partition_pos() /* moving_res_*/ );
	} else {
		add_virtual_O2Prime_hydrogen( pose );
	}

	if ( options_->sampler_perform_phosphate_pack() ){
		phosphate_sampler_ = new phosphate::MultiPhosphateSampler( pose );
		utility::vector1< Size > const & working_moving_partition_pos_ = job_parameters_->working_moving_partition_pos();
		if ( working_moving_partition_pos_.size() != 1 ){
			std::cerr << pose.fold_tree() << std::endl;
			std::cerr << pose.annotated_sequence() << std::endl;
			std::cerr << "MOVING_RES " << moving_res_ << std::endl;
			std::cerr << "PARTITION_POS " << working_moving_partition_pos_ << std::endl;
		}
		runtime_assert( working_moving_partition_pos_.size() == 1 ); //generalize later.
		phosphate_sampler_->set_moving_partition_res( working_moving_partition_pos_ );
	}

	//	atr_rep_checker_ = new checker::AtrRepChecker( pose_with_virtual_O2prime_hydrogen, job_parameters_ );
	atr_rep_checker_ = new checker::AtrRepChecker( *screening_pose_, job_parameters_ );
	atr_rep_checker_with_instantiated_sugar_ = new checker::AtrRepChecker( *sugar_screening_pose_, job_parameters_ );

	for ( Size n = 1; n <= five_prime_chain_breaks_.size(); n++ ) chain_closable_geometry_checkers_.push_back( new checker::ChainClosableGeometryChecker( five_prime_chain_breaks_[n], three_prime_chain_breaks_[n], chain_break_gap_sizes_[n] ) );

	for ( Size n = 1; n <= cutpoints_closed_.size(); n++ ) {
		chain_closure_checkers_.push_back( new checker::ChainClosureChecker( pose, cutpoints_closed_[n] ) );
		chain_closure_checkers_[n]->set_reinitialize_CCD_torsions( options_->reinitialize_CCD_torsions() );
	}

	initialize_moving_residue_pose_list( pose );

	num_pose_kept_to_use_ =  options_->sampler_num_pose_kept();
	if ( is_dinucleotide_ && base_centroid_checker_->allow_base_pair_only_screen() ) num_pose_kept_to_use_ = 4 * options_->sampler_num_pose_kept();
	pose_selection_ = new StepWiseRNA_PoseSelection( job_parameters_, scorefxn_ );
	pose_selection_->set_num_pose_kept( num_pose_kept_to_use_ );
	pose_selection_->set_cluster_rmsd( options_->cluster_rmsd() );
	pose_selection_->set_PBP_clustering_at_chain_closure( options_->PBP_clustering_at_chain_closure() );
	pose_selection_->set_distinguish_pucker( options_->distinguish_pucker() );
	// Allows for accumulation and clustering of poses across multiple jobs with e.g. different conformations of sampled sugars.
	pose_selection_->set_pose_list( pose_list_ );

	native_rmsd_screen_ = options_->sampler_native_rmsd_screen(); // may get updated later.

	max_ntries_ = std::max( 100000, 1000 * int( options_->num_random_samples() ) );
	if ( chain_closure_checkers_.size() > 0 ) max_ntries_ *= 10;
	if ( native_rmsd_screen_ ) max_ntries_ *= 10;

}

/////////////////////////////////////////////////////////////////////////////////////
Size
StepWiseRNA_RigidBodyConnectionSampler::which_residue_alternative_set_is_moving_residue() const {
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
StepWiseRNA_RigidBodyConnectionSampler::initialize_moving_residue_pose_list( pose::Pose const & pose ){
	if ( which_residue_alternative_set_is_moving_residue() > 0 ) return; // already initialized

	utility::vector1< pose::PoseOP > pose_list = setup_pose_with_moving_residue_alternative_list( pose, moving_res_, options_->extra_chi(), options_->use_phenix_geo() );
	rotamer_sampler::copy_dofs::ResidueAlternativeSet residue_alternative_set( pose_list, moving_res_ );
	residue_alternative_sets_.push_back( residue_alternative_set );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::reinstantiate_backbone_and_add_constraint_at_moving_res(	pose::Pose & pose, Size const & five_prime_chain_break_res )
{ //harmonic angle and distance constraints are used ONLY by chainbreak_screening
	if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) ) pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res_ ); //May 31, 2010
	if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_O2PRIME_HYDROGEN" ) ) pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", moving_res_ );
	if ( moving_res_ == ( five_prime_chain_break_res + 1 ) && pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_PHOSPHATE" ) ) {
		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", moving_res_ ); //this virtual_phosphate was added to pose at the beginning of this function.
	}
	add_harmonic_chain_break_constraint( pose, five_prime_chain_break_res );
	runtime_assert( !pose.residue( five_prime_chain_break_res   ).has_variant_type( "VIRTUAL_RIBOSE" ) );
	runtime_assert( !pose.residue( five_prime_chain_break_res+1 ).has_variant_type( "VIRTUAL_RIBOSE" ) );
}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_euler_angle_grid_parameters(){

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// definition of euler angle grid search parameters
	// Following are set by  #define in StepWiseRNA_FloatingBase_Samper_util.hh
	euler_angle_bin_size_ = euler_angle_bin_size;
	euler_z_bin_size_     = euler_z_bin_size;
	centroid_bin_size_    = centroid_bin_size;
	if ( options_->integration_test_mode() ){ // use coarser search for speed
		euler_angle_bin_size_ *= 4;
		//euler_z_bin_size_     *= 4;
		centroid_bin_size_    *= 4;
	}

	euler_angle_bin_min_ =  - 180/euler_angle_bin_size_    ; //Should be -180/euler_angle_bin_size
	euler_angle_bin_max_ =    180/euler_angle_bin_size_ - 1;  //Should be 180/euler_angle_bin_size-1
	euler_z_bin_min_ = int(  - 1/euler_z_bin_size_ );
	euler_z_bin_max_ = int(    1/euler_z_bin_size_ );
}

//////////////////////////////////////////////////////////////////////
// max_distance is currently hard-wired... instead we should make it
// depend on gap_size_to_anchor_, right?
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_xyz_grid_parameters(){
	max_distance_ = options_->sampler_max_centroid_distance(); // if unspecified (0.0), will be replaced

	initialize_xyz_parameters( max_distance_, max_distance_squared_,
														 centroid_bin_min_, centroid_bin_max_,
														 get_moving_rsd_list(), truly_floating_base() );
	base_bin_map_.clear();
}

//////////////////////////////////////////////////////////////////////
// used in setting max_distance -- can use a tighter tether if moving res is covalently connected to the reference (anchor) residue.
Size  // Size instead of bool for historical reasons.
StepWiseRNA_RigidBodyConnectionSampler::truly_floating_base(){
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
StepWiseRNA_RigidBodyConnectionSampler::get_moving_rsd_list() const {
	runtime_assert( rigid_body_sampling_ );
	Size const n = which_residue_alternative_set_is_moving_residue();
	runtime_assert( n > 0 );
	utility::vector1< core::conformation::ResidueOP > moving_rsd_list;
	rotamer_sampler::copy_dofs::ResidueAlternativeSet const & moving_res_modeling = residue_alternative_sets_[ n ];
	for ( Size n = 1; n <= moving_res_modeling.pose_list().size(); n++ ){
		moving_rsd_list.push_back( moving_res_modeling.pose( n )->residue( moving_res_ ).clone() );
	}
	return moving_rsd_list;
}

//////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::initialize_rigid_body_sampler( pose::Pose const & pose ){

	using namespace rotamer_sampler::rigid_body;
	using namespace rotamer_sampler::copy_dofs;

	// defining this RigidBodyRotamer also produces a useful pose with the entire moving partition
	// transformed to the origin.
	RigidBodyRotamerOP rigid_body_rotamer =  define_rigid_body_rotamer( pose );

	ResidueAlternativeRotamerCombOP rsd_alternatives_rotamer = new ResidueAlternativeRotamerComb();

	// note that following will include moving_res_ for sampler, as well as any other chunks that might move...
	for ( Size n = 1; n <= residue_alternative_sets_.size(); n++ ){
		ResidueAlternativeSet const & residue_alternative_set = residue_alternative_sets_[n];
		ResidueAlternativeRotamerOP rsd_alt_rotamer = new ResidueAlternativeRotamer( residue_alternative_set,
																																								 *rigid_body_rotamer->pose_at_origin() /*take representative residues from this pose after applying copy_dofs*/);
		rsd_alternatives_rotamer->add_residue_alternative_rotamer( rsd_alt_rotamer );
	}

	sampler_ = new RigidBodyRotamerWithResidueAlternatives( rsd_alternatives_rotamer, rigid_body_rotamer );

	sampler_->set_random( options_->choose_random() );
	sampler_->init();
}

//////////////////////////////////////////////////////////////////////
rotamer_sampler::rigid_body::RigidBodyRotamerOP
StepWiseRNA_RigidBodyConnectionSampler::define_rigid_body_rotamer( pose::Pose const & pose ){

	using namespace rotamer_sampler::rigid_body;
	RigidBodyRotamerOP sampler = new rotamer_sampler::rigid_body::RigidBodyRotamer( pose, moving_res_ );

	Real max_distance_rounded = int( max_distance_/centroid_bin_size_ ) * centroid_bin_size_;
	sampler->set_x_values( -max_distance_rounded + 0.5 * Real(centroid_bin_size_),
													+max_distance_rounded - 0.5 * Real(centroid_bin_size_),
													centroid_bin_size_);
	sampler->set_y_values( -max_distance_rounded + 0.5 * Real(centroid_bin_size_),
													+max_distance_rounded - 0.5 * Real(centroid_bin_size_),
													centroid_bin_size_);

	// weird offset -- just matching Parin's original
	sampler->set_z_values( -max_distance_rounded,
													+max_distance_rounded - Real(centroid_bin_size_),
													centroid_bin_size_);

	Real const max_angle_rounded = int( 180 / euler_angle_bin_size_ ) * euler_angle_bin_size_;
	sampler->set_euler_alpha_values( -max_angle_rounded + 0.5 * Real ( euler_angle_bin_size_ ),
																		+max_angle_rounded - 0.5 * Real ( euler_angle_bin_size_ ),
																		euler_angle_bin_size_ );

	Real const max_euler_z_rounded = int( 1.0 / euler_z_bin_size_ ) * euler_z_bin_size_;
	sampler->set_euler_z_values( -max_euler_z_rounded,
																+max_euler_z_rounded,
																euler_z_bin_size_ );
	sampler->set_euler_gamma_values( -max_angle_rounded + 0.5 * Real ( euler_angle_bin_size_ ),
																		+max_angle_rounded - 0.5 * Real ( euler_angle_bin_size_ ),
																		euler_angle_bin_size_ );

	sampler->set_random( options_->choose_random() );

	return sampler;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::kinematics::Stub
StepWiseRNA_RigidBodyConnectionSampler::get_reference_stub( pose::Pose const & pose ) const{

	std::string const reference_stub_type = "base"; //"sugar"
	core::kinematics::Stub reference_stub;

	TR.Debug << "-----------------------get reference stub-----------------------" << std::endl;
	if ( reference_stub_type == "sugar" ){
		reference_stub = get_sugar_stub( pose.residue( reference_res_ ), is_prepend_, true );
	} else{ //Use the base
		reference_stub.v = core::chemical::rna::get_rna_base_centroid(  pose.residue( reference_res_ ), options_->verbose() );
		reference_stub.M = core::chemical::rna::get_rna_base_coordinate_system( pose.residue( reference_res_ ), reference_stub.v );
	}

	TR.Debug << " reference_stub.v: x = " << reference_stub.v[0] << " y = " << reference_stub.v[1] << " z = " << reference_stub.v[2] << std::endl;
	TR.Debug << "---------------------------------------------------------------------" << std::endl;

	return reference_stub;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker ){
	base_centroid_checker_ = checker;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

void
StepWiseRNA_RigidBodyConnectionSampler::add_residue_alternative_set( rotamer_sampler::copy_dofs::ResidueAlternativeSet const & residue_alternative_set ){
	residue_alternative_sets_.push_back( residue_alternative_set );
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseRNA_RigidBodyConnectionSampler::get_pose_list(){
	return pose_list_;
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_RigidBodyConnectionSampler::set_options( StepWiseRNA_ModelerOptionsCOP options ){
	options_ = options;
}





} //rigid_body
} //legacy
} //rna
} //sampling
} //stepwise
} //protocols
