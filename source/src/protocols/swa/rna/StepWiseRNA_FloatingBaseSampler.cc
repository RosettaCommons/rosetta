// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/StepWiseRNA_FloatingBaseSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSamplerUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/rna/screener/ChainClosableScreener.hh>
#include <protocols/swa/rna/screener/ChainBreakScreener.hh>
#include <protocols/swa/rna/screener/AtrRepScreener.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <protocols/swa/rna/O2PrimePacker.hh>
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
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>


using ObjexxFCL::string_of;
using ObjexxFCL::lead_zero_string_of;

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_FloatingBaseSampler" ) ;

using namespace core;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// "Floating base", a.k.a. dinucleotide or "skip bulge" sampling, developed by Parin Sripakdeevong.
//
// Intead of enumerating torsions for a residue directly attached to the fixed pose (as in StepWiseRNA_StandardResidueSampler),
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
// As with StepWiseRNA_StandardResidueSampler, this Sampler uses several poses and screeners to reduce the number of variant changes &
// energy recomputations -- faster code but somewhat complicated.
//
// Note that sugar of the floating base will be virtualized unless it needs to be instantiated to close chain to
//  anchor or distal residues. Nevertheless, code below carried out a geometric 'sanity check' that involves temporarily
//  instantiating the sugar.
//
// On the way to dinucleotide sampling, we often start with poses in which a virtual sugar in the attachment point in the fixed
//  pose was instantiated. One option would be to run this dinucleotide sampler separately for each of those instantiation cases.
//  That was the choice made in StandardResidueSampler -- and the code is simpler there. But here, that information is passed into
//  the sampler in 'anchor_sugar_modeling'; some of the first screens are based on base centroid distances and are independent of
//  sugar conformation, and so some computation that would be redundant across anchor sugar conformations can be avoided.
//
//   -- Rhiju, 2013.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace swa {
namespace rna {

//Constructor
StepWiseRNA_FloatingBaseSampler::StepWiseRNA_FloatingBaseSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
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
	num_pose_kept_( 108 ),
	num_pose_kept_to_use_( num_pose_kept_ ),
	cluster_rmsd_( 0.5001 ),
	verbose_( false ),
	native_rmsd_screen_( false ),
	native_screen_rmsd_cutoff_( 2.0 ),
	perform_o2prime_pack_( true ),
	integration_test_mode_( false ), //March 16, 2012
	centroid_screen_( true ),
	VDW_atr_rep_screen_( true ),
	distinguish_pucker_( true ),
	PBP_clustering_at_chain_closure_( false ), //New option Aug 15 2010
	reinitialize_CCD_torsions_( false ), //New option Aug 15 2010 //Reinitialize_CCD_torsion to zero before every CCD chain closure
	extra_chi_( false ),
	choose_random_( false ), // Rhiju, Jul 2013
	num_random_samples_( 1 ),
	use_phenix_geo_( false ),
	max_ntries_( 0 ), // updated below
	euler_angle_bin_min_( 0 ), // updated below
	euler_angle_bin_max_( 0 ), // updated below
	euler_z_bin_min_( 0 ), // updated below
	euler_z_bin_max_( 0 ),  // updated below
	centroid_bin_min_( 0 ),  // updated below
	centroid_bin_max_( 0 ),  // updated below
	max_distance_( 0.0 ),  // updated below
	max_distance_squared_( 0.0 ), // updated below
	try_sugar_instantiation_( false ),
	o2prime_instantiation_distance_cutoff_( 6.0 )
{
	set_native_pose( job_parameters_->working_native_pose() );
	runtime_assert ( !is_dinucleotide_ || !is_internal_ );
	TR.Debug << "GAP_SIZE_TO_ANCHOR " << gap_size_to_anchor_ << "  REFERENCE RES        " << reference_res_ << "  MOVING_RES " << moving_res_ << std::endl;
	TR.Debug << "GAP_SIZE_TO_DISTAL " << gap_size_           << "  FIVE' CHAINBREAK RES " << five_prime_chain_break_res_ << "  MOVING_RES " << moving_res_ << std::endl;
}

//Destructor
StepWiseRNA_FloatingBaseSampler::~StepWiseRNA_FloatingBaseSampler()
{}

/////////////////////
std::string
StepWiseRNA_FloatingBaseSampler::get_name() const {
	return "StepWiseRNA_FloatingBaseSampler";
}

///////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::apply( core::pose::Pose & pose ){
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace protocols::rotamer_sampler::rigid_body;

	output_title_text( "Enter StepWiseRNA_FloatingBaseSampler::floating_base_sampling", TR.Debug );

	clock_t const time_start( clock() );

	initialize_poses_and_stubs_and_screeners( pose );
	initialize_euler_angle_grid_parameters();
	initialize_xyz_grid_parameters();
	initialize_rigid_body_sampler( pose );
	initialize_other_residues_base_list( pose ); 	// places where floating base can 'dock'

	utility::vector1< PoseOP > pose_list;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MAIN LOOP
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size ntries( 0 ), num_success( 0 ); // used in choose_random mode.
	for ( sampler_->reset(); sampler_->not_end(); ++(*sampler_) ){

		if ( choose_random_ && ++ntries > max_ntries_ ) break;
		if ( break_early_for_integration_tests() ) break;

		core::kinematics::Stub moving_res_base_stub = sampler_->get_stub();
		if ( ( moving_res_base_stub.v - reference_stub_.v ).length_squared() > max_distance_squared_ ) {
			sampler_->fast_forward_to_next_translation(); continue;
		}
		if ( centroid_screen_ &&
				 ( !base_centroid_screener_->check_centroid_interaction( moving_res_base_stub, count_data_ ) ||
					 !base_centroid_screener_->check_that_terminal_res_are_unstacked() ) ) {
			sampler_->fast_forward_to_next_euler_gamma(); continue;
		}
		count_data_.tot_rotamer_count++;

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// geometry checks that require base conformation (but not sugar) at moving_res
		if ( anchor_sugar_modeling_.sample_sugar ){
			if ( !chain_closable_to_anchor_screener_->check_screen( anchor_sugar_modeling_.pose_list, moving_rsd_at_origin_list_, moving_res_base_stub, reference_res_ ) ) continue;
		} else {
			if ( !chain_closable_to_anchor_screener_->check_screen( *screening_pose_, moving_rsd_at_origin_list_, moving_res_base_stub, reference_res_ ) ) continue;
		}
		if ( close_chain_to_distal_ ){
			if ( !chain_closable_to_distal_screener_->check_screen( *screening_pose_, moving_rsd_at_origin_list_, moving_res_base_stub, chain_break_reference_res_ ) ) continue;
		}
		count_data_.chain_closable_count++;

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// clash checks
		// Some trickiness with clashes to 3'-phosphate in case of chain closure -- parin has worked out these cases, fortunately.
		Residue const & screening_moving_rsd_at_origin = ( *screening_moving_rsd_at_origin_list_[ 5 /*magic number?*/ ] );
		if ( !VDW_bin_screener_->VDW_rep_screen( *screening_pose_, moving_res_, screening_moving_rsd_at_origin, moving_res_base_stub ) ) continue;
		// User-input VDW: Does not work for chain_closure move and is_internal_ move yet, since the screener does not know that
		// moving residue atoms can bond to previous or next residues.
		if ( ( user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ) && ( !close_chain_to_distal_ ) && ( !is_internal_ )  &&
				 !user_input_VDW_bin_screener_->VDW_rep_screen( *screening_pose_, moving_res_, screening_moving_rsd_at_origin, moving_res_base_stub ) ) continue;
		count_data_.good_bin_rep_count++;

		//////////////////////////////////////////////////////////////////////////////////////////////
		// finally copy in the base conformation into the pose for RMSD and VDW  screens
		sampler_->apply( *screening_pose_, screening_moving_rsd_at_origin );

		//////////////////////////////////////////////////////
		if ( native_rmsd_screen_ && get_native_pose() ){
			//This assumes that screening_pose and native_pose are already superimposed.
			if ( suite_rmsd( *get_native_pose(), *screening_pose_, moving_res_, is_prepend_, true /*ignore_virtual_atom*/ ) > ( native_screen_rmsd_cutoff_ ) ) continue;
			if ( rmsd_over_residue_list( *get_native_pose(), *screening_pose_, job_parameters_, true /*ignore_virtual_atom*/ ) > ( native_screen_rmsd_cutoff_ ) ) continue; //Oct 14, 2010
			count_data_.rmsd_count++;
		}

		if ( VDW_atr_rep_screen_ && !atr_rep_screener_->check_screen( *screening_pose_ ) ) continue;

		//////////////////////////////////////////////////////////////////////////
		// Inner-most loop: over potential sugar/chi conformations.
		//  Note that actual sugar will end up virtual -- this
		//  is a sanity check that at least one sugar conformation is viable
		//  in terms of sterics and chain closure geometry.
		//
		// Moving_rsd_at_origin_list and anchor_sugar_modeling_ are similar in spirit,
		//  storing possible sugar conformations for moving residue and for anchor residue.
		// But they are unfortunately different data structures here.
		for ( Size n = 1; n <= moving_rsd_at_origin_list_.size(); n++ ){

			std::string tag = "U_" + lead_zero_string_of( count_data_.tot_rotamer_count, 12 ) + '_' + string_of( n );

			sampler_->apply( *sugar_screening_pose_, ( *sugar_screening_moving_rsd_at_origin_list_[n] ) );
			if ( close_chain_to_distal_ ) if ( !chain_closable_to_distal_screener_->check_screen( *sugar_screening_pose_ ) ) continue;
			if ( !screen_anchor_sugar_conformation( pose, tag ) ) continue; // this loops over anchor sugar models if there are more than one.
			count_data_.non_clash_sugar++;

			// above was all on screening poses; the actual base of the pose hasn't been moved yet!
			conformation::Residue const & moving_rsd_at_origin( *moving_rsd_at_origin_list_[n] );
			if ( close_chain_to_distal_ ) sampler_->apply( chain_break_to_distal_screener_->pose(), moving_rsd_at_origin );
			if ( close_chain_to_anchor_ ) sampler_->apply( chain_break_to_anchor_screener_->pose(), moving_rsd_at_origin );
			if ( perform_o2prime_pack_ )  sampler_->apply( o2prime_packer_->pose(), moving_rsd_at_origin );
			sampler_->apply( pose, moving_rsd_at_origin );

			//////////////////////////////////////////////////////////////////////////////////////
			// Try to close chain break from floating base to a distal connection point, if it exists.
			if ( close_chain_to_distal_ ) {
				if ( !chain_break_to_distal_screener_->check_screen() ) continue;
				if ( perform_o2prime_pack_ ) chain_break_to_distal_screener_->copy_CCD_torsions( o2prime_packer_->pose() );
				chain_break_to_distal_screener_->copy_CCD_torsions( pose );
			}

			/////////////////////////////////////////////////////////////////////////////////////
			// Try to close chain break from floating base to its anchor point, if they are directly connected.
			if ( close_chain_to_anchor_ ){
				if ( !chain_closable_to_anchor_screener_->check_screen( chain_break_to_anchor_screener_->pose(), true /*strict*/ ) ) continue;
				if ( !chain_break_to_anchor_screener_->check_screen() ) continue;
				if ( perform_o2prime_pack_ ) chain_break_to_anchor_screener_->copy_CCD_torsions( o2prime_packer_->pose() );
				chain_break_to_anchor_screener_->copy_CCD_torsions( pose );
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////
			bool instantiate_sugar = check_moving_sugar( pose );
			if ( perform_o2prime_pack_ ){
				if ( instantiate_sugar ) instantiate_moving_sugar_and_o2prime( o2prime_packer_->pose() );
				o2prime_packer_->sample_o2prime_hydrogen();
				o2prime_packer_->copy_all_o2prime_torsions( pose ); //Copy the o2prime torsions from the o2prime_pack_pose to the pose!
				if ( instantiate_sugar ) virtualize_moving_sugar_and_o2prime( o2prime_packer_->pose() );
			}

			Pose selected_pose = pose; // the reason for this copy is that we might apply a bulge variant, and that can produce thread conflicts with graphics.
			if ( instantiate_sugar ) instantiate_moving_sugar_and_o2prime( selected_pose );
			pose_selection_->pose_selection_by_full_score( selected_pose, tag );

			if ( verbose_ ) output_data( silent_file_, tag, true, pose, get_native_pose(), job_parameters_ );
			num_success++;

			// As long as we're not instantiating the sugar in the final pose, break once we found any valid sugar rotamer
			if ( !close_chain_to_distal_ && !close_chain_to_anchor_ && !try_sugar_instantiation_ ) break;

		} // floating base's sugar/chi rotamer

		update_base_bin_map( sampler_->get_values() ); // diagnostics.

		if ( choose_random_ && num_success >= num_random_samples_ ) break;
	}

	if ( choose_random_ ) TR << "Number of tries: " << ntries << ".  Passed centroid filters: " << count_data_.tot_rotamer_count++ << ". Number of successes: " << num_success << std::endl;

	pose_selection_->finalize();
	output_count_data();

	// can put this back in pretty easily...
	if ( verbose_ ) analyze_base_bin_map( base_bin_map_, "test/" );

	pose_list_ = pose_selection_->pose_list();
	TR.Debug << "floating base sampling time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::initialize_poses_and_stubs_and_screeners( pose::Pose & pose  ){

	runtime_assert( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_PHOSPHATE" ) );
	runtime_assert( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) );
	if ( gap_size_to_anchor_ == 0 ) runtime_assert( pose.fold_tree().is_cutpoint( moving_res_ - 1 ) );

	// hydrogens will be reinstantiated later. perhaps this should be the screening_pose?
	Pose pose_with_virtual_O2prime_hydrogen = pose;
	add_virtual_O2Prime_hydrogen( pose_with_virtual_O2prime_hydrogen );

	reference_stub_ = get_reference_stub( pose_with_virtual_O2prime_hydrogen );

	////////////////////////////////////////Screeners///////////////////////////////////////////////////
	VDW_bin_screener_ = new screener::StepWiseRNA_VDW_BinScreener();
	VDW_bin_screener_->setup_using_working_pose( pose_with_virtual_O2prime_hydrogen, job_parameters_ );
	if ( user_input_VDW_bin_screener_ ) user_input_VDW_bin_screener_->reference_xyz_consistency_check( reference_stub_.v );

	screening_pose_ = pose_with_virtual_O2prime_hydrogen.clone(); //Hard copy, used for clash checking
	if ( close_chain_to_distal_ ) pose::add_variant_type_to_pose_residue( *screening_pose_, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res_ + 1 ); //May 31, 2010

	sugar_screening_pose_ = screening_pose_->clone(); //Hard copy. Used for trying out sugar at moving residue.
	pose::remove_variant_type_from_pose_residue( *sugar_screening_pose_, "VIRTUAL_RIBOSE", moving_res_ );

	if ( close_chain_to_distal_ ) reinstantiate_backbone_and_add_constraint_at_moving_res( pose, five_prime_chain_break_res_ );
	if ( close_chain_to_anchor_ ) reinstantiate_backbone_and_add_constraint_at_moving_res( pose, floating_base_five_prime_chain_break_ );

	if ( perform_o2prime_pack_ ) {
		remove_virtual_O2Prime_hydrogen( pose );
		// weird -- following should not be necessary, but commenting it out changes results.
		if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) ) pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", moving_res_ );
		o2prime_packer_ = new O2PrimePacker( pose, scorefxn_, job_parameters_->working_moving_partition_pos() /* moving_res_*/ );
	} else {
		add_virtual_O2Prime_hydrogen( pose );
	}

	atr_rep_screener_ = new screener::AtrRepScreener( pose_with_virtual_O2prime_hydrogen, job_parameters_ );
	atr_rep_screener_with_instantiated_sugar_ = new screener::AtrRepScreener( *sugar_screening_pose_, job_parameters_ );

	// this seems like overkill, but anchor sugar modeling can involve minimizing, and the 'anchor' residue moves.
	for ( Size anchor_sugar_ID = 1; anchor_sugar_ID <= anchor_sugar_modeling_.pose_list.size(); anchor_sugar_ID++ ){
		pose::Pose const & anchor_sugar_modeling_pose = *anchor_sugar_modeling_.pose_list[anchor_sugar_ID];
		copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, *sugar_screening_pose_, anchor_sugar_modeling_pose );
		atr_rep_screeners_for_anchor_sugar_models_.push_back( new screener::AtrRepScreener( *sugar_screening_pose_, job_parameters_ ) );
	}

	chain_closable_to_distal_screener_ = new screener::ChainClosableScreener( five_prime_chain_break_res_, gap_size_ );
	chain_break_to_distal_screener_    = new screener::ChainBreakScreener( pose, five_prime_chain_break_res_ ); //Hard copy
	chain_break_to_distal_screener_->set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );

	chain_closable_to_anchor_screener_ = new screener::ChainClosableScreener( floating_base_five_prime_chain_break_, floating_base_three_prime_chain_break_, gap_size_to_anchor_ );
	chain_break_to_anchor_screener_    = new screener::ChainBreakScreener( pose, floating_base_five_prime_chain_break_ ); //Hard copy
	chain_break_to_anchor_screener_->set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup Residue of moving and reference of various rsd conformation (syn/anti chi, 2' and 3' endo) with
	// base at origin coordinate frame
	moving_rsd_at_origin_list_                 = setup_residue_at_origin_list(	pose,	moving_res_, extra_chi_, use_phenix_geo_ );
	screening_moving_rsd_at_origin_list_	     = setup_residue_at_origin_list( *screening_pose_, moving_res_, extra_chi_, use_phenix_geo_ );
	sugar_screening_moving_rsd_at_origin_list_ = setup_residue_at_origin_list( *sugar_screening_pose_, moving_res_,	extra_chi_, use_phenix_geo_ );
	runtime_assert ( moving_rsd_at_origin_list_.size() == screening_moving_rsd_at_origin_list_.size() );
	runtime_assert ( moving_rsd_at_origin_list_.size() == sugar_screening_moving_rsd_at_origin_list_.size() );

	num_pose_kept_to_use_ =  num_pose_kept_;
	if ( is_dinucleotide_ && base_centroid_screener_->allow_base_pair_only_screen() ) num_pose_kept_to_use_ = 4 * num_pose_kept_;
	pose_selection_ = new StepWiseRNA_PoseSelection( job_parameters_, scorefxn_ );
	pose_selection_->set_num_pose_kept( num_pose_kept_to_use_ );
	pose_selection_->set_cluster_rmsd( cluster_rmsd_ );
	pose_selection_->set_PBP_clustering_at_chain_closure( PBP_clustering_at_chain_closure_ );
	pose_selection_->set_distinguish_pucker( distinguish_pucker_ );

	max_ntries_ = std::max( 100000, 1000 * int( num_random_samples_ ) );
	if ( close_chain_to_distal_ ) max_ntries_ *= 10;
	if ( close_chain_to_anchor_ ) max_ntries_ *= 10;
	if ( native_rmsd_screen_ ) max_ntries_ *= 10;

}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::reinstantiate_backbone_and_add_constraint_at_moving_res(	pose::Pose & pose, Size const & five_prime_chain_break_res )
{ //harmonic angle and distance constraints are used ONLY by chainbreak_screening
	if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) ) pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res_ ); //May 31, 2010
	if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_O2PRIME_HYDROGEN" ) ) pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", moving_res_ );
	if ( moving_res_ == ( five_prime_chain_break_res + 1 ) && pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_PHOSPHATE" ) ) {
		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", moving_res_ ); //this virtual_phosphate was added to pose at the beginning of this function.
	}
	add_harmonic_chain_break_constraint( pose, five_prime_chain_break_res );
}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::initialize_euler_angle_grid_parameters(){

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// definition of euler angle grid search parameters
	// Following are set by  #define in StepWiseRNA_FloatingBase_Samper_util.hh
	euler_angle_bin_size_ = euler_angle_bin_size;
	euler_z_bin_size_     = euler_z_bin_size;
	centroid_bin_size_    = centroid_bin_size;
	if ( integration_test_mode_ ){ // use coarser search for speed
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
StepWiseRNA_FloatingBaseSampler::initialize_xyz_grid_parameters(){

	Distance C5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list_, " C5'" );
	Distance O5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list_, " O3'" );
	Distance const Max_O3_to_C5_DIST = ( gap_size_to_anchor_ == 0 ) ? O3I_C5I_PLUS_ONE_MAX_DIST : O3I_C5I_PLUS_TWO_MAX_DIST;

	//Theoretical maximum dist between the two base's centroid, +1 is to be lenient
	max_distance_ = Max_O3_to_C5_DIST + C5_centroid_dist + O5_centroid_dist + 1.0;
	max_distance_squared_ = max_distance_ * max_distance_;

	centroid_bin_min_ = int(  - max_distance_/centroid_bin_size );
	centroid_bin_max_ = int( max_distance_/centroid_bin_size ) - 1;

	base_bin_map_.clear();
}

//////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::initialize_rigid_body_sampler( pose::Pose const & pose ){

	sampler_ = new rotamer_sampler::rigid_body::RigidBodyRotamer( moving_res_, pose.residue( moving_res_ ), reference_stub_ );

	Real max_distance_rounded = int( max_distance_/centroid_bin_size_ ) * centroid_bin_size_;
	sampler_->set_x_values( -max_distance_rounded + 0.5 * Real(centroid_bin_size_),
													+max_distance_rounded - 0.5 * Real(centroid_bin_size_),
													centroid_bin_size_);
	sampler_->set_y_values( -max_distance_rounded + 0.5 * Real(centroid_bin_size_),
													+max_distance_rounded - 0.5 * Real(centroid_bin_size_),
													centroid_bin_size_);

	// weird offset -- just matching Parin's original
	sampler_->set_z_values( -max_distance_rounded,
													+max_distance_rounded - Real(centroid_bin_size_),
													centroid_bin_size_);

	Real const max_angle_rounded = int( 180 / euler_angle_bin_size_ ) * euler_angle_bin_size_;
	sampler_->set_euler_alpha_values( -max_angle_rounded + 0.5 * Real ( euler_angle_bin_size_ ),
																		+max_angle_rounded - 0.5 * Real ( euler_angle_bin_size_ ),
																		euler_angle_bin_size_ );

	Real const max_euler_z_rounded = int( 1.0 / euler_z_bin_size_ ) * euler_z_bin_size_;
	sampler_->set_euler_z_values( -max_euler_z_rounded,
																+max_euler_z_rounded,
																euler_z_bin_size_ );
	sampler_->set_euler_gamma_values( -max_angle_rounded + 0.5 * Real ( euler_angle_bin_size_ ),
																		+max_angle_rounded - 0.5 * Real ( euler_angle_bin_size_ ),
																		euler_angle_bin_size_ );

	sampler_->set_random( choose_random_ );

	sampler_->init();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_FloatingBaseSampler::break_early_for_integration_tests() {
	if ( integration_test_mode_ && count_data_.both_count > 10000 )     return true;
	if ( integration_test_mode_ && count_data_.full_score_count >= 10 ) native_rmsd_screen_ = true;
	if ( integration_test_mode_ && count_data_.rmsd_count >= 10 )       return true;
	return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// following function might be dramatically simplifiable, in which case we may want to return
// it to main loop.
bool
StepWiseRNA_FloatingBaseSampler::screen_anchor_sugar_conformation( pose::Pose & pose, std::string & tag ){

	//OK check that with this sugar, the chain can be theoretically closed..
	std::string const moving_atom_name    = ( is_prepend_ ) ? " O3'" : " C5'";
	std::string const reference_atom_name = ( is_prepend_ ) ? " C5'" : " O3'";

	if ( !anchor_sugar_modeling_.sample_sugar ){
		if ( !chain_closable_to_anchor_screener_->check_screen( *sugar_screening_pose_ ) ) return false;
		if ( VDW_atr_rep_screen_ && !atr_rep_screener_with_instantiated_sugar_->check_screen( *sugar_screening_pose_ ) ) return false; // wait a minute... why is this in here? oh, because base can move in different sampled sugar modeling conformations.
		return true;
	} else {

		// The point of this section is to look for *any* conformation of sugar in anchor residue that passes
		// screens, going from the lowest energy option on up.
		//Ok, since anchor_sugar_modeling_.pose_list is sorted by SCORE, the lower energy conformations are tried first!

		for ( Size n = 1; n <= anchor_sugar_modeling_.pose_list.size(); n++ ){

			pose::Pose const & anchor_sugar_modeling_pose = *anchor_sugar_modeling_.pose_list[n];

			// is_prepend --> sugar_screening_pose [moving] = 5' pose [need O3'], anchor_sugar = 3' pose [need C5']
			if ( !chain_closable_to_anchor_screener_->check_screen( *sugar_screening_pose_, anchor_sugar_modeling_pose, is_prepend_) ) continue;

			copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, *sugar_screening_pose_, anchor_sugar_modeling_pose );
			// This is in here because the anchor sugar models can have slightly shifted bases (not just riboses!) due to a minimization step that
			// can occur in VirtualRiboseSampler [see the option: do_minimize].
			if ( VDW_atr_rep_screen_ && !atr_rep_screeners_for_anchor_sugar_models_[ n ]->check_screen( *sugar_screening_pose_ ) ) continue;

			copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, pose, anchor_sugar_modeling_pose );
			if ( perform_o2prime_pack_ ) copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, o2prime_packer_->pose(), anchor_sugar_modeling_pose );
			tag += tag_from_pose( anchor_sugar_modeling_pose );
			return true;
		}
		return false;
	}

	runtime_assert( false /*should not get here*/ )
	return true; // shouldn't actually get here.
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::instantiate_moving_sugar_and_o2prime( pose::Pose & pose ){
	if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) )           pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res_ ); //May 31, 2010
	if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_O2PRIME_HYDROGEN" ) ) pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", moving_res_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::virtualize_moving_sugar_and_o2prime( pose::Pose & pose ){
	if ( !pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) ) pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res_ ); //May 31, 2010
	if ( !pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_O2PRIME_HYDROGEN" ) ) pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", moving_res_ );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_FloatingBaseSampler::check_moving_sugar( pose::Pose & pose ){
	if ( !try_sugar_instantiation_ ) return false;
	if ( !pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) ) return false; // nothing to do.

	bool instantiate_sugar( false );
	Vector const & moving_O2prime_xyz = pose.residue( moving_res_ ).xyz( " O2'" );
	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		if ( i == moving_res_ ) continue;
		if ( pose.residue( i ).is_virtual_residue() ) continue;
		for ( Size j = 1; j <= pose.residue( i ).nheavyatoms(); j++ ){
			if ( pose.residue_type( i ).is_virtual( j ) ) continue;
			if ( pose.residue_type( i ).heavyatom_is_an_acceptor( j ) ||
					 pose.residue_type( i ).heavyatom_has_polar_hydrogens( j )  ){
				Distance dist = ( pose.residue(i).xyz( j ) - moving_O2prime_xyz ).length();
				//std::cout << "Distance to rsd " << i << "  atom " << pose.residue_type( i ).atom_name( j ) << ": " << dist << std::endl;
				if ( dist < o2prime_instantiation_distance_cutoff_ ) {
					instantiate_sugar = true; break;
				}
			}
		} // atoms j
		if ( instantiate_sugar ) break;
	} // residues i

	return instantiate_sugar;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::kinematics::Stub
StepWiseRNA_FloatingBaseSampler::get_reference_stub( pose::Pose const & pose ) const{

	std::string const reference_stub_type = "base"; //"sugar"
	core::kinematics::Stub reference_stub;

	TR.Debug << "-----------------------get reference stub-----------------------" << std::endl;
	if ( reference_stub_type == "sugar" ){
		reference_stub = get_sugar_stub( pose.residue( reference_res_ ), is_prepend_, true );
	} else{ //Use the base
		reference_stub.v = core::chemical::rna::get_rna_base_centroid(  pose.residue( reference_res_ ), verbose_ );
		reference_stub.M = core::chemical::rna::get_rna_base_coordinate_system( pose.residue( reference_res_ ), reference_stub.v );
	}

	TR.Debug << " reference_stub.v: x = " << reference_stub.v[0] << " y = " << reference_stub.v[1] << " z = " << reference_stub.v[2] << std::endl;
	TR.Debug << "---------------------------------------------------------------------" << std::endl;

	return reference_stub;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::initialize_other_residues_base_list( pose::Pose const & pose ) {

	runtime_assert( max_distance_ > 0.0 );
	other_residues_base_list_.clear();
	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

		conformation::Residue const & residue_object = pose.residue( seq_num );
		if ( residue_object.aa() == core::chemical::aa_vrt ) continue;
		if ( !residue_object.is_RNA() ) continue;
		if ( residue_object.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
		if ( seq_num == job_parameters_->working_moving_res() ) continue;
		if ( ( job_parameters_->working_terminal_res() ).has_value( seq_num ) ) continue;

		core::kinematics::Stub base_info;
		base_info.v = core::chemical::rna::get_rna_base_centroid( residue_object, verbose_ );
		base_info.M = core::chemical::rna::get_rna_base_coordinate_system( residue_object, base_info.v );

		Real const distance = ( base_info.v - core::chemical::rna::get_rna_base_centroid( pose.residue( reference_res_ ), false ) ).length();

		if (  distance > ( max_distance_ + 6.364 ) ) continue; //6.364 is max centroid to centroid distance for base-stacking centroid_screening.

		TR.Debug << " Added to other_residues_base_list: " << seq_num << std::endl;
		other_residues_base_list_.push_back( base_info );
	}

}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener ){
	base_centroid_screener_ = screener;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

void
StepWiseRNA_FloatingBaseSampler::set_anchor_sugar_modeling( SugarModeling const & anchor_sugar_modeling ){
	anchor_sugar_modeling_ = anchor_sugar_modeling;
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseRNA_FloatingBaseSampler::get_pose_list(){
	return pose_list_;
}

/////////////////////////////////////////////////////////////////////////////////////
// diagnostics
void
StepWiseRNA_FloatingBaseSampler::update_base_bin_map( Base_bin const & base_bin ){

	std::map< Base_bin, int, compare_base_bin > ::const_iterator it = base_bin_map_.find( base_bin );
	if ( it == base_bin_map_.end() )	base_bin_map_[base_bin] = 0;
	base_bin_map_[base_bin] ++;
}

/////////////////////////////////////////////////////////////////////////////////////
// diagnostics
void
StepWiseRNA_FloatingBaseSampler::update_base_bin_map( utility::vector1< Real > const & rigid_body_values ){
	Base_bin base_bin;
	base_bin.centroid_x  = rigid_body_values[6];
	base_bin.centroid_y  = rigid_body_values[5];
	base_bin.centroid_z  = rigid_body_values[4];
	base_bin.euler_alpha = rigid_body_values[3];
	base_bin.euler_z     = rigid_body_values[2];
	base_bin.euler_gamma = rigid_body_values[1];
	update_base_bin_map( base_bin );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::output_count_data(){
	Size const gap_size_( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */

	if ( gap_size_ == 0 ){
		TR.Debug << " angle_n = " << count_data_.good_angle_count << " dist_n = " << count_data_.good_distance_count;
		TR.Debug << " chain_break_screening = " << count_data_.chain_break_screening_count << std::endl;
	}

	TR.Debug << " stack_n = " << count_data_.base_stack_count << " pair_n = " << count_data_.base_pairing_count;
	TR.Debug << " strict_pair_n = " << count_data_.strict_base_pairing_count << " centroid_n = " << count_data_.pass_base_centroid_screen;
	TR.Debug << " bin_rep = " << count_data_.good_bin_rep_count << " atr = " << count_data_.good_atr_rotamer_count << " rep = " << count_data_.good_rep_rotamer_count;
	TR.Debug << " both = " << count_data_.both_count << " total_bin = " << count_data_.tot_rotamer_count << " total_screen_bin = " << base_bin_map_.size();
	TR.Debug << "  closable = " << count_data_.chain_closable_count << "  non_clash_sugar = " << count_data_.non_clash_sugar << std::endl;
	TR.Debug << " WARNING centroid_n count is severely UNDERESTIMATED...need to be multiply by ( euler_angle_bin_max_ - euler_angle_bin_min_ + 1 ): ";
	TR.Debug << ( euler_angle_bin_max_ - euler_angle_bin_min_ + 1 ) << std::endl;
}

void
StepWiseRNA_FloatingBaseSampler::set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){ user_input_VDW_bin_screener_ = user_input_VDW_bin_screener; }

} //rna
} //swa
} //protocols
