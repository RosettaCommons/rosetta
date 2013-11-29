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
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/rna/screener/ChainBreakScreener.hh>
#include <protocols/swa/rna/screener/AtrRepScreener.hh>
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

Real const RADS_PER_DEG = numeric::NumericTraits < Real > ::pi() / 180.;

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
	five_prime_chain_break_res_( job_parameters_->five_prime_chain_break_res() ),
	num_nucleotides_( job_parameters_->working_moving_res_list().size() ),
	reference_res_( job_parameters_->working_reference_res() ), //the last static_residues that this attach to the moving residues
	floating_base_five_prime_chain_break_( ( is_prepend_ ) ? moving_res_ : moving_res_ - 1 ), //for floating base chain closure when num_nucleotides=1
	is_dinucleotide_( num_nucleotides_ == 2 ),
	close_chain_to_distal_( gap_size_ == 0 ),
	close_chain_to_anchor_( num_nucleotides_ == 1 ),
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
	allow_base_pair_only_centroid_screen_( false ), //allow for possibility of conformation that base_pair but does not base_stack
	VDW_atr_rep_screen_( true ),
	distinguish_pucker_( true ),
	PBP_clustering_at_chain_closure_( false ), //New option Aug 15 2010
	reinitialize_CCD_torsions_( false ), //New option Aug 15 2010 //Reinitialize_CCD_torsion to zero before every CCD chain closure
	extra_chi_( false ),
	choose_random_( false ), // Rhiju, Jul 2013
	num_random_samples_( 1 ),
	use_phenix_geo_( false ),
	euler_angle_bin_min_( 0 ), // updated below
	euler_angle_bin_max_( 0 ), // updated below
	euler_z_bin_min_( 0 ), // updated below
	euler_z_bin_max_( 0 ),  // updated below
	centroid_bin_min_( 0 ),  // updated below
	centroid_bin_max_( 0 ),  // updated below
	max_distance_( 0.0 ),  // updated below
	max_distance_squared_( 0.0 ) // updated below
{
	set_native_pose( job_parameters_->working_native_pose() );
	if ( is_dinucleotide_ && is_internal_ ) utility_exit_with_message( "is_dinucleotide_ == true && is_internal_ == true )!!" );
	if ( num_nucleotides_ != 1 && num_nucleotides_ != 2 ) utility_exit_with_message( "num_nucleotides != 1 and num_nucleotides != 2" );
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

	output_title_text( "Enter StepWiseRNA_FloatingBaseSampler::floating_base_sampling", TR.Debug );

	clock_t const time_start( clock() );

	initialize_poses_and_stubs_and_screeners( pose );
	initialize_euler_angle_grid_parameters(); // later will move into a rotamer sampler.
	initialize_xyz_grid_parameters();
	initialize_other_residues_base_list( pose ); 	// places where floating base can 'dock'

	core::kinematics::Stub moving_res_base_stub;
	Euler_angles euler_angles;
	Base_bin base_bin;
	utility::vector1< PoseOP > pose_data_list;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MAIN LOOP
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( base_bin.euler_alpha = euler_angle_bin_min_; base_bin.euler_alpha <= euler_angle_bin_max_; base_bin.euler_alpha++ ){
		for ( base_bin.euler_z = euler_z_bin_min_; base_bin.euler_z <= euler_z_bin_max_; base_bin.euler_z++ ){

			if ( break_early_for_integration_tests() ) break;

			// Note: most of the initial screening occurs at the level of residues, not in poses.
			Matrix O_frame_rotation;
			euler_angles.alpha = ( base_bin.euler_alpha + 0.5 )*euler_angle_bin_size*( RADS_PER_DEG ); //convert to radians
			euler_angles.z = ( base_bin.euler_z )*euler_z_bin_size;
			euler_angles.beta = acos( euler_angles.z );
			euler_angles.gamma = 0;

			convert_euler_to_coordinate_matrix( euler_angles, O_frame_rotation );
			moving_res_base_stub.M = reference_stub_.M * O_frame_rotation;

			////////////////////////////////////////////////////////////////////////////////////////////////////////
			for ( base_bin.centroid_z = centroid_bin_min_; base_bin.centroid_z <= centroid_bin_max_; base_bin.centroid_z++ ){
				for ( base_bin.centroid_x = centroid_bin_min_; base_bin.centroid_x <= centroid_bin_max_; base_bin.centroid_x++ ){
					for ( base_bin.centroid_y = centroid_bin_min_; base_bin.centroid_y <= centroid_bin_max_; base_bin.centroid_y++ ){

						Vector O_frame_centroid;
						O_frame_centroid[0] = ( base_bin.centroid_x + 0.5 ) * centroid_bin_size;
						O_frame_centroid[1] = ( base_bin.centroid_y + 0.5 ) * centroid_bin_size;
						O_frame_centroid[2] = ( base_bin.centroid_z ) * centroid_bin_size;
						moving_res_base_stub.v = ( reference_stub_.M * O_frame_centroid ) + reference_stub_.v;

						////////////////////Not dependent on euler gamma value//////////////////////////////////////////////////////////////
						if ( ( moving_res_base_stub.v - reference_stub_.v ).length_squared() > max_distance_squared_ ) continue;
						if ( centroid_screen_  &&  !floating_base_centroid_screening( moving_res_base_stub, other_residues_base_list_, num_nucleotides_, count_data_, allow_base_pair_only_centroid_screen_ ) ) continue;

						//////////////////////Update the moving_res_base_stub -- euler angle is last /////////////////////////////////////////
						for ( base_bin.euler_gamma = euler_angle_bin_min_; base_bin.euler_gamma <= euler_angle_bin_max_; base_bin.euler_gamma++ ){

							count_data_.tot_rotamer_count++;

							euler_angles.gamma = ( base_bin.euler_gamma + 0.5 )*euler_angle_bin_size*( RADS_PER_DEG ); //convert to radians
							convert_euler_to_coordinate_matrix( euler_angles, O_frame_rotation );
							moving_res_base_stub.M = reference_stub_.M * O_frame_rotation;

							//////////////////////////////////////////////////////////////////////////////////////////////////////
							// geometry checks that do not require sugar instantiated at moving_res -- can we close chains?
							if ( anchor_sugar_modeling_.sample_sugar ){
								if (	!check_floating_base_chain_closable( reference_res_, anchor_sugar_modeling_.PDL, moving_rsd_at_origin_list_, moving_res_base_stub, is_prepend_, ( num_nucleotides_ - 1 ) ) ) continue;
							} else{
								if (	!check_floating_base_chain_closable( reference_res_, *screening_pose_, moving_rsd_at_origin_list_, moving_res_base_stub, is_prepend_, ( num_nucleotides_ - 1 ) ) ) continue;
							}
							if ( close_chain_to_distal_ ){
								Size const chain_break_reference_res_ = ( is_prepend_ ) ? five_prime_chain_break_res_ : five_prime_chain_break_res_ + 1;
								if ( !check_floating_base_chain_closable( chain_break_reference_res_, *screening_pose_, moving_rsd_at_origin_list_, moving_res_base_stub, !is_prepend_, 0 /*gap_size_*/ ) ) continue;
							}
							count_data_.chain_closable_count++;

							//////////////////////////////////////////////////////////////////////////////////////////////////////
							// clash checks
							Residue const & screening_moving_rsd_at_origin = ( *screening_moving_rsd_at_origin_list_[ 5 /*magic number?*/ ] );
							// Some trickiness with clashes to 3'-phosphate in case of chain closure -- parin has worked out these cases, fortunately.
							if ( !VDW_bin_screener_->VDW_rep_screen( *screening_pose_, moving_res_, screening_moving_rsd_at_origin, moving_res_base_stub ) ) continue;
							// User-input VDW: Does not work for chain_closure move and is_internal_ move yet, since the screener does not know that moving residue atoms can bond
							//  to previous or next residues.
							if ( ( user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ) && ( !close_chain_to_distal_ ) && ( !is_internal_ )  &&
									 !user_input_VDW_bin_screener_->VDW_rep_screen( *screening_pose_, moving_res_, screening_moving_rsd_at_origin, moving_res_base_stub ) ) continue;
							count_data_.good_bin_rep_count++;

							//////////////////////////////////////////////////////////////////////////////////////////////
							// finally copy in the base conformation into the pose for RMSD and VDW  screens
							set_base_coordinate_frame( *screening_pose_, moving_res_, screening_moving_rsd_at_origin, moving_res_base_stub );

							//////////////////////////////////////////////////////
							if ( native_rmsd_screen_ && get_native_pose() ){
								//This assumes that screening_pose and native_pose are already superimposed.
								if ( suite_rmsd( *get_native_pose(), *screening_pose_, moving_res_, is_prepend_, true /*ignore_virtual_atom*/ ) > ( native_screen_rmsd_cutoff_ ) ) continue;
								if ( rmsd_over_residue_list( *get_native_pose(), *screening_pose_, job_parameters_, true /*ignore_virtual_atom*/ ) > ( native_screen_rmsd_cutoff_ ) ) continue; //Oct 14, 2010
								count_data_.rmsd_count++;
							}

							if ( VDW_atr_rep_screen_ && !atr_rep_screener_->check_screen( *screening_pose_, gap_size_, is_internal_, false /*kic_sampling_*/ ) ) continue;
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

								set_base_coordinate_frame( *sugar_screening_pose_, moving_res_, ( *sugar_screening_moving_rsd_at_origin_list_[n] ), moving_res_base_stub );
								if ( gap_size_ == 0 ) if ( !check_chain_closable( *sugar_screening_pose_, five_prime_chain_break_res_, 0 ) ) continue;
								if ( !screen_anchor_sugar_conformation( pose, tag ) ) continue; // this loops over anchor sugar models if there are more than one.
								count_data_.non_clash_sugar++;

								// above was all on screening poses; the actual base of the pose hasn't been moved yet!
								conformation::Residue const & moving_rsd_at_origin( *moving_rsd_at_origin_list_[n] );
								if ( close_chain_to_distal_ ) set_base_coordinate_frame( chain_break_screener_->pose(), moving_res_, moving_rsd_at_origin, moving_res_base_stub );
								if ( close_chain_to_anchor_ ) set_base_coordinate_frame( chain_break_screener_to_anchor_->pose(), moving_res_, moving_rsd_at_origin, moving_res_base_stub );
								if ( perform_o2prime_pack_ ) set_base_coordinate_frame( o2prime_packer_->pose(), moving_res_, moving_rsd_at_origin, moving_res_base_stub );
								set_base_coordinate_frame( pose, moving_res_, moving_rsd_at_origin, moving_res_base_stub );

								//////////////////////////////////////////////////////////////////////////////////////
								// Try to close chain break from floating base to a connection point, if it exists.
								if ( close_chain_to_distal_ ) {
									if ( !chain_break_screener_->check_screen() ) continue;
									if ( perform_o2prime_pack_ ) chain_break_screener_->copy_CCD_torsions( o2prime_packer_->pose() );
									chain_break_screener_->copy_CCD_torsions( pose );
								}

								/////////////////////////////////////////////////////////////////////////////////////
								// Try to close chain break from floating base to its anchor point, if they are directly connected.
								if ( close_chain_to_anchor_ ){
									if ( !check_chain_closable_floating_base( chain_break_screener_to_anchor_->pose(), chain_break_screener_to_anchor_->pose(), floating_base_five_prime_chain_break_, 0 ) ) continue; //strict version of the check chain closable.
									if ( !chain_break_screener_to_anchor_->check_screen() ) continue;
									if ( perform_o2prime_pack_ ) chain_break_screener_to_anchor_->copy_CCD_torsions( o2prime_packer_->pose() );
									chain_break_screener_to_anchor_->copy_CCD_torsions( pose );
								}

								////////////////////////////////////////////////////////////////////////////////////////////////////////
								if ( perform_o2prime_pack_ ){
									o2prime_packer_->sample_o2prime_hydrogen();
									o2prime_packer_->copy_all_o2prime_torsions( pose ); //Copy the o2prime torsions from the o2prime_pack_pose to the pose!
								}

								pose_selection_->pose_selection_by_full_score( pose, tag );
								if ( verbose_ ) output_data( silent_file_, tag, true, pose, get_native_pose(), job_parameters_ );

								// As long as we're not instantiating the sugar in the final pose, break once we found any valid sugar rotamer
								if ( !close_chain_to_distal_ && !close_chain_to_anchor_ ) break;

							} // floating base's sugar/chi rotamer

							update_base_bin_map( base_bin); // diagnostics.

						} // base centroid x
					} // base centroid y
				} // base centroid z
			} // euler gamma
		} // euler z
	} // euler alpha

	pose_selection_->finalize();
	output_count_data();

	if ( verbose_ ) analyze_base_bin_map( base_bin_map_, "test/" );

	pose_data_list_ = pose_selection_->pose_data_list();
	TR.Debug << "floating base sampling time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_FloatingBaseSampler::initialize_poses_and_stubs_and_screeners( pose::Pose & pose  ){

	runtime_assert( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_PHOSPHATE" ) );
	runtime_assert( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) );
	if ( num_nucleotides_ == 1 ) runtime_assert( pose.fold_tree().is_cutpoint( moving_res_ - 1 ) );

	// hydrogens will be reinstantiated later. perhaps this should be the screening_pose?
	Pose pose_with_virtual_O2prime_hydrogen = pose;
	add_virtual_O2Star_hydrogen( pose_with_virtual_O2prime_hydrogen );

	atr_rep_screener_ = new screener::AtrRepScreener( pose_with_virtual_O2prime_hydrogen, job_parameters_ );

	reference_stub_ = get_reference_stub( pose_with_virtual_O2prime_hydrogen );

	////////////////////////////////////////Screeners///////////////////////////////////////////////////
	VDW_bin_screener_ = new screener::StepWiseRNA_VDW_BinScreener();
	VDW_bin_screener_->setup_using_working_pose( pose_with_virtual_O2prime_hydrogen, job_parameters_ );
	if ( user_input_VDW_bin_screener_ ) user_input_VDW_bin_screener_->reference_xyz_consistency_check( reference_stub_.v );

	screening_pose_ = pose.clone(); //Hard copy, used for clash checking
	if ( close_chain_to_distal_ ) pose::add_variant_type_to_pose_residue( *screening_pose_, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res_ + 1 ); //May 31, 2010

	sugar_screening_pose_ = screening_pose_->clone(); //Hard copy. Used for trying out sugar at moving residue.
	pose::remove_variant_type_from_pose_residue( *sugar_screening_pose_, "VIRTUAL_RIBOSE", moving_res_ );

	if ( close_chain_to_distal_ ) reinstantiate_backbone_and_add_constraint_at_moving_res( pose, five_prime_chain_break_res_ );
	if ( close_chain_to_anchor_ ) reinstantiate_backbone_and_add_constraint_at_moving_res( pose, floating_base_five_prime_chain_break_ );

	if ( perform_o2prime_pack_ ) {
		remove_virtual_O2Star_hydrogen( pose );
		// weird -- following should not be necessary, but commenting it out changes results.
		if ( pose.residue( moving_res_ ).has_variant_type( "VIRTUAL_RIBOSE" ) ) pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", moving_res_ );
		o2prime_packer_ = new O2PrimePacker( pose, scorefxn_, job_parameters_->working_moving_partition_pos() /* moving_res_*/ );
	} else {
		add_virtual_O2Star_hydrogen( pose );
	}

	chain_break_screener_ = new screener::ChainBreakScreener( pose, five_prime_chain_break_res_ ); //Hard copy
	chain_break_screener_->set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );

	chain_break_screener_to_anchor_ = new screener::ChainBreakScreener( pose, floating_base_five_prime_chain_break_ ); //Hard copy
	chain_break_screener_to_anchor_->set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup Residue of moving and reference of various rsd conformation (syn/anti chi, 2' and 3' endo) with
	// base at origin coordinate frame
	moving_rsd_at_origin_list_                 = setup_residue_at_origin_list(	pose,	moving_res_, extra_chi_, use_phenix_geo_ );
	screening_moving_rsd_at_origin_list_	     = setup_residue_at_origin_list( *screening_pose_, moving_res_, extra_chi_, use_phenix_geo_ );
	sugar_screening_moving_rsd_at_origin_list_ = setup_residue_at_origin_list( *sugar_screening_pose_, moving_res_,	extra_chi_, use_phenix_geo_ );
	runtime_assert ( moving_rsd_at_origin_list_.size() == screening_moving_rsd_at_origin_list_.size() );
	runtime_assert ( moving_rsd_at_origin_list_.size() == sugar_screening_moving_rsd_at_origin_list_.size() );

	num_pose_kept_to_use_ =  num_pose_kept_;
	if ( is_dinucleotide_ && allow_base_pair_only_centroid_screen_ ) num_pose_kept_to_use_ = 4 * num_pose_kept_;
	pose_selection_ = new StepWiseRNA_PoseSelection( job_parameters_, scorefxn_ );
	pose_selection_->set_num_pose_kept( num_pose_kept_to_use_ );
	pose_selection_->set_cluster_rmsd( cluster_rmsd_ );
	pose_selection_->set_PBP_clustering_at_chain_closure( PBP_clustering_at_chain_closure_ );
	pose_selection_->set_distinguish_pucker( distinguish_pucker_ );
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
		// instead this should be its own class, and these should be private variables.
		Real euler_angle_bin_size_ = euler_angle_bin_size;
		Real euler_z_bin_size_     = euler_z_bin_size;
		Real centroid_bin_size_    = centroid_bin_size;
		if ( integration_test_mode_ ){ // use coarser search for speed
			euler_angle_bin_size_ *= 4;
			//euler_z_bin_size_     *= 4;
			centroid_bin_size_    *= 4;
		}

		euler_angle_bin_min_ =  - 180/euler_angle_bin_size_; //Should be -180/euler_angle_bin_size
		euler_angle_bin_max_ = 180/euler_angle_bin_size_ - 1;  //Should be 180/euler_angle_bin_size-1
		euler_z_bin_min_ = int(  - 1/euler_z_bin_size_ );
		euler_z_bin_max_ = int( 1/euler_z_bin_size_ );
}


/////////////////////////////////////////////////////////////////////////////////////////
// Following defines a grid search centered around 'takeoff' position for floating base.
void
StepWiseRNA_FloatingBaseSampler::initialize_xyz_grid_parameters(){
	Distance C5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list_, " C5'" );
	Distance O5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list_, " O3'" );
	Distance const Max_O3_to_C5_DIST = ( num_nucleotides_ == 1 ) ? O3I_C5I_PLUS_ONE_MAX_DIST : O3I_C5IPLUS2_MAX_DIST;

	max_distance_ = Max_O3_to_C5_DIST + C5_centroid_dist + O5_centroid_dist + 1.0; //Theoretical maximum dist between the two base's centroid, +1 is to be lenient
	max_distance_squared_ = max_distance_ * max_distance_;

	centroid_bin_min_ = int(  - max_distance_/centroid_bin_size );
	centroid_bin_max_ = int( max_distance_/centroid_bin_size ) - 1;
	base_bin_map_.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_FloatingBaseSampler::break_early_for_integration_tests() {
	if ( integration_test_mode_ && count_data_.both_count > 10000 ) return true;
	if ( integration_test_mode_ && count_data_.full_score_count >= 10 ) native_rmsd_screen_ = true;
	if ( integration_test_mode_ && count_data_.rmsd_count >= 10 ) return true;
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
		if ( !check_chain_closable( sugar_screening_pose_->residue( moving_res_ ).xyz( moving_atom_name ),
																sugar_screening_pose_->residue( reference_res_ ).xyz( reference_atom_name ), ( num_nucleotides_ - 1 ) ) ) return false;
		if ( VDW_atr_rep_screen_ && !atr_rep_screener_->check_screen( *sugar_screening_pose_, gap_size_, is_internal_, false /*kic_sampling_*/ ) ) return false;
		return true;
	} else {
		// The point of this section is to look for *any* conformation of sugar in attached residue that passes
		// screens, going from the lowest energy option on up.
		//Ok, since anchor_sugar_modeling_.PDL is sorted by SCORE, the lower energy conformations are tried first!
		for ( Size anchor_sugar_ID = 1; anchor_sugar_ID <= anchor_sugar_modeling_.PDL.size(); anchor_sugar_ID++ ){

			if ( !check_chain_closable( sugar_screening_pose_->residue( moving_res_ ).xyz( moving_atom_name ),
																	( *anchor_sugar_modeling_.PDL[anchor_sugar_ID] ).residue( reference_res_ ).xyz( reference_atom_name ), ( num_nucleotides_ - 1 ) ) )	 continue;

			copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, *sugar_screening_pose_, ( *anchor_sugar_modeling_.PDL[anchor_sugar_ID] ) );
			if ( VDW_atr_rep_screen_ && !atr_rep_screener_->check_screen( *sugar_screening_pose_, gap_size_, is_internal_, false /*kic_sampling_*/ ) ) continue;

			copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, pose, ( *anchor_sugar_modeling_.PDL[anchor_sugar_ID] ) );
			if ( perform_o2prime_pack_ ) copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, o2prime_packer_->pose(), ( *anchor_sugar_modeling_.PDL[anchor_sugar_ID] ) );
			tag += tag_from_pose( *anchor_sugar_modeling_.PDL[ anchor_sugar_ID ] );
			return true;
		}
		return false;
	}

	runtime_assert( false /*should not get here*/ )
	return true; // shouldn't actually get here.
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
StepWiseRNA_FloatingBaseSampler::get_pose_data_list(){
	return pose_data_list_;
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
