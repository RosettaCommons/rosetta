// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/suite/StepWiseRNA_ClassicResidueSampler.cc
/// @brief Sampler of a single residue, the heart of StepWise Assembly of RNA
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Parin Sripakdeevong, sripakpa@stanford.edu


#include <protocols/stepwise/enumerate/rna/suite/StepWiseRNA_ClassicResidueSampler.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/enumerate/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/enumerate/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/enumerate/rna/phosphate/PhosphateUtil.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/stepwise/enumerate/rna/screener/AtrRepScreener.hh>
#include <protocols/stepwise/enumerate/rna/screener/ChainClosableScreener.hh>
#include <protocols/stepwise/enumerate/rna/screener/ChainBreakScreener.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/rotamer_sampler/rna/RNA_KicSampler.hh>
#include <protocols/rotamer_sampler/rna/RNA_SuiteRotamer.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

using ObjexxFCL::string_of;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Enumerates over suite backbone, sugar pucker, and base chi angle; and packs 2' hydroxyls, for a single RNA residue
// at the terminus of an RNA pose.
//
//                 epsilon,zeta,alpha,
//                    beta,gamma      pucker
//
//     5'   --Sugar --  Phosphate -- Sugar  ... -- Phosphate -- Sugar -- Phosphate -- 3'
//              |                      | chi                       |
//          Reference               Moving                       Distal
//             Base                 Residue                      Base
//                                                        |
//                                       If no residues to distal, close_chain_to_distal
//
//                                     |<------ gap_size + 1 ----->|
//
//
// Refactored out of StepWiseRNA_ResidueSampler.cc
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_ClassicResidueSampler" ) ;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace suite {

//Constructor
StepWiseRNA_ClassicResidueSampler::StepWiseRNA_ClassicResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
		job_parameters_( job_parameters ),
		moving_res_(  job_parameters_->working_moving_res() ), // Might not corresponds to user input.
		moving_suite_(  job_parameters_->working_moving_suite() ), // dofs betweeen this value and value+1 actually move.
		is_prepend_(  job_parameters_->is_prepend() ),
		is_internal_(  job_parameters_->is_internal() ), // no cutpoints before or after moving_res_.
		actually_moving_res_( job_parameters_->actually_moving_res() ), //Now same as moving_res_
		gap_size_( job_parameters_->gap_size() ), /* If this is zero or one, need to screen or closable chain break */
		working_moving_partition_pos_( job_parameters_->working_moving_partition_pos() ),
		num_nucleotides_(  job_parameters_->working_moving_res_list().size() ),
		is_dinucleotide_( num_nucleotides_ == 2 ),
		close_chain_to_distal_( gap_size_ == 0 ),
		five_prime_chain_break_res_( job_parameters_->five_prime_chain_break_res() ),
		last_append_res_( ( is_prepend_ ) ? moving_res_ - 1: moving_res_ ),
		last_prepend_res_( ( is_prepend_ ) ? moving_res_: moving_res_ + 1 ),
		atom_atom_overlap_dist_cutoff_(-1.0 ),
		extra_tag_( "" ),
		build_pose_from_scratch_( job_parameters_->working_sequence().length() == ( num_nucleotides_ + 1 ) ), // somewhat hacky, used for rna puzzle
		kic_sampling_( false ), // will be updated below.
		rebuild_bulge_mode_( job_parameters_->rebuild_bulge_mode() ),
		scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ), // can be replaced from the outside
		silent_file_( "silent_file.txt" ),
		native_rmsd_screen_( false ), // will be updated below
		bin_size_( 20 )
{
	set_native_pose( job_parameters_->working_native_pose() );
	if ( num_nucleotides_ != 1 && num_nucleotides_ != 2 )		utility_exit_with_message( "num_nucleotides_ != 1 and num_nucleotides_ != 2" );
}

//Destructor
StepWiseRNA_ClassicResidueSampler::~StepWiseRNA_ClassicResidueSampler()
{}

	/////////////////////
std::string
StepWiseRNA_ClassicResidueSampler::get_name() const {
	return "StepWiseRNA_ClassicResidueSampler";
}

////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::apply( core::pose::Pose & pose ) {

	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;

	output_title_text( "Enter StepWiseRNA_ClassicResidueSampler::standard_sampling", TR.Debug );

	clock_t const time_start( clock() );

	initialize_poses_and_screeners( pose );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MAIN LOOP --> rotamer sampling.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size ntries( 0 ), num_success( 0 ); // used in choose_random mode.
	Size max_ntries = std::max( 10000, 100 * int( options_->num_random_samples() ) );
	if ( kic_sampling_ ) max_ntries = 5 * options_->num_random_samples(); // some chains just aren't closable.

	for ( rotamer_sampler_->reset(); rotamer_sampler_->not_end(); ++( *rotamer_sampler_ ) ) {

		if ( options_->choose_random() && ++ntries > max_ntries ) break;

		rotamer_sampler_->apply( *screening_pose_ );
		count_data_.tot_rotamer_count++;

		if ( break_early_for_integration_tests() ) break;

		std::string tag = create_tag( "U" + extra_tag_, ntries );

		// For KIC closure, immediate check that a closed loop solution was actually found.
		if ( kic_sampling_ && !chain_break_screener_->check_loop_closed( *screening_pose_ ) ) continue;

		if ( native_rmsd_screen_ && get_native_pose() ){
			if ( suite_rmsd( *get_native_pose(), *screening_pose_, actually_moving_res_, is_prepend_ ) > ( options_->sampler_native_screen_rmsd_cutoff() ) ) continue;
			if ( rmsd_over_residue_list( *get_native_pose(), *screening_pose_, job_parameters_, false ) > ( options_->sampler_native_screen_rmsd_cutoff() ) ) continue; //Oct 14, 2010
			count_data_.rmsd_count++;
		}

		if ( options_->combine_long_loop_mode() && !close_chain_to_distal_ ){ //residue-residue contact screen;
			//Nov 18,2010
			if ( !is_residues_in_contact( last_append_res_, *screening_pose_, last_prepend_res_, *screening_pose_, atom_atom_overlap_dist_cutoff_, 1 /*num_atom_contacts_cutoff*/ ) ) continue;
			count_data_.residues_contact_screen++; //mistakenly put this inside the if loop, fix on Sept 22, 2010
		}

		bool is_possible_bulge = false;

		if ( base_centroid_screener_ ){
			//Reminder note of independency: Important that base_stub_list is updated even in the case where gap_size_ == 0 (and bulge is allowed) ////
			// since base_stub_list is used later below in the chain_break_screening section Jan 28, 2010 Parin S. ///////////////////////////////////
			bool found_a_centroid_interaction_partner( false );
			found_a_centroid_interaction_partner = base_centroid_screener_->update_base_stub_list_and_check_centroid_interaction( *screening_pose_, count_data_ );

			if ( close_chain_to_distal_ && !found_a_centroid_interaction_partner ){ //does not stack or base_pair
				if ( working_moving_partition_pos_.size() == 1 ) is_possible_bulge = true;
			}
			//Essentially this doesn't screen for centroid interaction at chain_break.
			if ( ( !close_chain_to_distal_ || options_->force_centroid_interaction() ) && !found_a_centroid_interaction_partner ) continue;
			if ( num_nucleotides_ > 1 && is_possible_bulge ) utility_exit_with_message( "num_nucleotides_ > 1 but is_possible_bulge == true!" );

			// Note that is does not update base_stub_list. To do that, use update_base_stub_list_and_check_that_terminal_res_are_unstacked
			if ( !base_centroid_screener_->check_that_terminal_res_are_unstacked() ) continue;
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		///////////////Chain_break_screening -- distance cut                     /////////////////
		//////////////////////////////////////////////////////////////////////////////////////////
		if ( close_chain_to_distal_ && !kic_sampling_ ){
			if ( !chain_closable_screener_->check_screen( *screening_pose_, options_->finer_sampling_at_chain_closure() ) ) continue;
			count_data_.chain_closable_count++;
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		/////////////// Van_der_Waals_screening                                  /////////////////
		//////////////////////////////////////////////////////////////////////////////////////////
		if ( options_->VDW_atr_rep_screen() && !atr_rep_screener_->check_screen( *screening_pose_ ) ) continue;

		if ( ( user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ) && ( gap_size_ != 0 ) && ( !is_internal_ ) ){
			//Does not work for chain_closure move and is_internal_ move yet...
			//Residue at 3' of building region have a phosphate that is NOT VIRTUALIZED. This Residue should not be excluded in VDW_bin_screen_pose! Feb 21, 2011.
			//Residue at 5' of building region also have O3' atom that is covalently bond to the phosphate atom of loop res next to it. VDW_rep doesn't realize this and this lead to clash. Hence This Residue should not be excluded in the VDW_bin_screen_pose as well. Feb 21, 2011.

			if ( !user_input_VDW_bin_screener_->VDW_rep_screen( *screening_pose_, moving_res_ ) ) continue;
			count_data_.good_bin_rep_count++;
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		// Almost ready to actually score pose.
		//////////////////////////////////////////////////////////////////////////////////////////
		rotamer_sampler_->apply( pose );
		if ( options_->sampler_perform_o2prime_pack() ) rotamer_sampler_->apply( o2prime_packer_->pose() );

		//////////////////////////////////////////////////////////////////////////////////////////
		///////////////Chain_break_screening -- CCD closure /////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////
		bool add_bulge( false );
		if ( close_chain_to_distal_ ) {
			if ( !kic_sampling_ ){ // kic sampling already should have closed chain.
				rotamer_sampler_->apply( chain_break_screener_->pose() );
				if ( ! chain_break_screener_->check_screen() ) continue;
				chain_break_screener_->copy_CCD_torsions( pose );
				if ( options_->sampler_perform_o2prime_pack() ) chain_break_screener_->copy_CCD_torsions( o2prime_packer_->pose() );
				if ( is_possible_bulge ){
					add_bulge = bulge_variant_decision( pose, atr_rep_screener_->delta_atr_score() ); /*further cut on atr, inside*/
				}
			}
		}

		if ( options_->sampler_perform_phosphate_pack() ){
			rotamer_sampler_->apply( phosphate_sampler_->pose() );
			if ( close_chain_to_distal_ ) chain_break_screener_->copy_CCD_torsions( phosphate_sampler_->pose() );
			phosphate_sampler_->sample_phosphates();
			phosphate_sampler_->copy_phosphates( o2prime_packer_->pose() );
			phosphate_sampler_->copy_phosphates( pose );
		}

		if ( options_->sampler_perform_o2prime_pack() ){
			if ( add_bulge ) apply_bulge_variant( o2prime_packer_->pose() );
			o2prime_packer_->sample_o2prime_hydrogen();
			o2prime_packer_->copy_all_o2prime_torsions( pose ); //Copy the o2prime torsions from the o2prime_pack_pose to the pose!
			if ( add_bulge ) remove_virtual_rna_residue_variant_type( o2prime_packer_->pose(), moving_res_ );
		}

		if ( options_->sampler_include_torsion_value_in_tag() ) tag += create_rotamer_string( pose, moving_res_, is_prepend_ );

		///////Add pose to pose_list if pose has good score///////////
		Pose selected_pose = pose; // the reason for this copy is that we might apply a variant, which can produce thread conflicts with graphics.

		if ( add_bulge ) apply_bulge_variant( selected_pose );
		pose_selection_->pose_selection_by_full_score( selected_pose, tag );

		TR.Debug << tag <<  std::endl;
		if (options_->verbose() ) output_data( silent_file_, tag, true, selected_pose, get_native_pose(), job_parameters_ );

		num_success++;
		if ( options_->choose_random() && num_success >= options_->num_random_samples() ) break;
	}

	if ( options_->choose_random() ) {
		TR << "Number of tries: " << count_data_.tot_rotamer_count++ << ". Number of successes: " << num_success <<  std::endl;
		TR.Debug << "Was shooting for max_tries: " << max_ntries << ". num_random_samples: " << options_->num_random_samples() << std::endl;
	}

	output_title_text( "Final sort and clustering", TR.Debug );
	pose_selection_->finalize( !build_pose_from_scratch_ /*do_clustering*/ );
	pose_list_ = pose_selection_->pose_list();

	output_count_data();

	TR.Debug << "Total time in StepWiseRNA_ClassicResidueSampler: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::initialize_poses_and_screeners( pose::Pose & pose ){

	if ( options_->combine_long_loop_mode() && !close_chain_to_distal_ ){ //residue-residue contact screen;
		TR.Debug << "options_->combine_long_loop_mode() && gap_size_ == 0" << std::endl;
		TR.Debug << "Enforcing contact between LAST_APPEND_RES: " << last_append_res_ << " and LAST_PREPEND_RES: " << last_prepend_res_  << std::endl;
		TR.Debug << "atom_atom_overlap_dist_cutoff " << atom_atom_overlap_dist_cutoff_ << std::endl;
	}

	if ( options_->sampler_perform_phosphate_pack() ){
		phosphate_sampler_ = new phosphate::MultiPhosphateSampler( pose );
		phosphate_sampler_->set_moving_partition_res( working_moving_partition_pos_ );
	}

	/////////////////////////////// O2prime sampling/virtualization //////////////////////////
	Pose pose_with_virtual_O2prime_hydrogen = pose;
	add_virtual_O2Prime_hydrogen( pose_with_virtual_O2prime_hydrogen );

	if ( options_->sampler_perform_o2prime_pack() ) {
		o2prime_packer_ = new o2prime::O2PrimePacker( pose, scorefxn_, working_moving_partition_pos_ /* moving_res_*/ );
		o2prime_packer_->set_use_green_packer( options_->use_green_packer() ); // green_packer not on by default
		o2prime_packer_->set_partition_definition( job_parameters_->partition_definition() ); // would be needed for green packer
	} else {
		// Otherwise, virtualize the 2-OH.
		pose = pose_with_virtual_O2prime_hydrogen;
	}

	kic_sampling_ = ( options_->kic_sampling_if_relevant() && close_chain_to_distal_ );

	screening_pose_ = pose_with_virtual_O2prime_hydrogen.clone(); //Hard copy
	phosphate::remove_terminal_phosphates( *screening_pose_ );

	// wait shouldn't this happen below? after all virtualization is complete?
	//atr_rep_screener_ = new screener::AtrRepScreener( *screening_pose_, job_parameters_ );
	//		atr_rep_screener_->set_kic_sampling( kic_sampling_ );

	//Necessary for the case where gap_size_ == 0. In this case, Pose_setup does not automatically create a VIRTUAL_PHOSPHATE.///
	//However since screening_pose is scored before the CCD correctly positions the chain_break phosphate atoms, ///////////////
	// the VIRTUAL_PHOSPHATE is needed to prevent artificial clashes. Parin Jan 28, 2010////////////////////////////////////
	if ( close_chain_to_distal_ )	pose::add_variant_type_to_pose_residue( *screening_pose_,      "VIRTUAL_PHOSPHATE",    five_prime_chain_break_res_ + 1 );

	atr_rep_screener_ = new screener::AtrRepScreener( *screening_pose_, job_parameters_ );
	atr_rep_screener_->set_kic_sampling( kic_sampling_ );

	chain_closable_screener_ = new screener::ChainClosableScreener( five_prime_chain_break_res_, gap_size_ );

	chain_break_screener_ = new screener::ChainBreakScreener( pose_with_virtual_O2prime_hydrogen, five_prime_chain_break_res_ ); //Hard copy
	chain_break_screener_->set_reinitialize_CCD_torsions( options_->reinitialize_CCD_torsions() );

	if ( close_chain_to_distal_ ) add_harmonic_chain_break_constraint( pose, five_prime_chain_break_res_ );

	pose_selection_ = new StepWiseRNA_PoseSelection( job_parameters_, scorefxn_ );
	pose_selection_->set_num_pose_kept( get_num_pose_kept() );
	pose_selection_->set_cluster_rmsd( options_->cluster_rmsd() );
	pose_selection_->set_PBP_clustering_at_chain_closure( options_->PBP_clustering_at_chain_closure() );
	pose_selection_->set_distinguish_pucker( options_->distinguish_pucker() );
	// Allows for accumulation and clustering of poses across multiple jobs with e.g. different conformations of sampled sugars.
	pose_selection_->set_pose_list( pose_list_ );

	native_rmsd_screen_ = options_->sampler_native_rmsd_screen();

	/////Get the Rotamer Sampler/////
	runtime_assert( base_centroid_screener_ );
	rotamer_sampler_ = setup_rotamer_sampler( *screening_pose_ );

	set_extra_tag( tag_from_pose( pose ) );

}


/////////////////////////////////////////////////////////////////////////////////
Size
StepWiseRNA_ClassicResidueSampler::get_num_pose_kept()
{
	Size num_pose_kept = options_->sampler_num_pose_kept();

	if ( build_pose_from_scratch_ ){
		TR.Debug << "Since build_pose_from_scratch, choose to increase NUM_POSE_KEPT by 36 fold. ";
		TR.Debug << " Somewhat hacky..since sample both sugar..want to make sure that we keep are good energy score states " << std::endl;
		TR.Debug << "Old_num_pose_kept = " << num_pose_kept  << std::endl;
		num_pose_kept = 36* num_pose_kept;
		TR.Debug << "New_num_pose_kept = " << num_pose_kept  << std::endl;
	}

	if ( job_parameters_->sample_both_sugar_base_rotamer() ){
		TR.Debug << "Since build_pose_from_scratch, choose to increase NUM_POSE_KEPT by 12 fold. ";
		TR.Debug << " Somewhat hacky..since sample both sugar..want to make sure that we keep are good energy score states " << std::endl;
		TR.Debug << "Old_num_pose_kept = " << num_pose_kept  << std::endl;
		num_pose_kept = 12* num_pose_kept;
		TR.Debug << "New_num_pose_kept = " << num_pose_kept  << std::endl;
	}

	return num_pose_kept;
}

//////////////////////////////////////////////////////////////////
// Fang-Chieh Chou nicely refactored the rotamer sampling
//////////////////////////////////////////////////////////////////
rotamer_sampler::RotamerBaseOP
StepWiseRNA_ClassicResidueSampler::setup_rotamer_sampler( pose::Pose const & pose ) const {

	using namespace rotamer_sampler::rna;
	using namespace chemical::rna;
	using namespace pose::rna;

	/////Load in constants being used/////
	utility::vector1<Size> const & working_moving_suite_list(
		job_parameters_->working_moving_suite_list() );
	utility::vector1<Size> const & syn_chi_res(
		job_parameters_->working_force_syn_chi_res_list() );
	utility::vector1<Size> const & north_puckers(
		job_parameters_->working_force_north_sugar_list() );
	utility::vector1<Size> const & south_puckers(
		job_parameters_->working_force_south_sugar_list() );

	runtime_assert( working_moving_suite_list.size() == 1 );
	Size const moving_suite_( working_moving_suite_list[1] );

	/////Get the base and pucker state/////
	utility::vector1<bool> sample_sugar( 2, false );
	utility::vector1<Size> base_state( 2, WHATEVER );
	utility::vector1<Size> pucker_state( 2, WHATEVER );

	if ( build_pose_from_scratch_ || job_parameters_->sample_both_sugar_base_rotamer()) {
		sample_sugar[1] = true;
		sample_sugar[2] = true;
	} else if ( !is_internal_  ) {
		if ( is_prepend_ ) {
			sample_sugar[1] = true;
		} else {
			sample_sugar[2] = true;
		}
	} else {
		runtime_assert( is_internal_ );
	}

	for ( Size i = 1; i <= 2; ++i ) {
		Size const curr_rsd( moving_suite_ + i - 1 );
		if ( sample_sugar[i] ) {
			bool is_north ( north_puckers.has_value( curr_rsd ) );
			bool is_south ( south_puckers.has_value( curr_rsd ) );
			bool is_syn ( syn_chi_res.has_value( curr_rsd ) );
			runtime_assert( !is_north || !is_south );
			if ( is_north ) pucker_state[i] = NORTH;
			if ( is_south ) pucker_state[i] = SOUTH;
			if ( !options_->allow_syn_pyrimidine() && !is_purine( pose.residue( curr_rsd ) ) ) {
				runtime_assert( !is_syn );
				base_state[i] = ANTI;
			} else if ( is_syn ) {
				base_state[i] = SYN;
			}
		} else {
			pucker_state[i] = assign_pucker( pose, curr_rsd );
			base_state[i] = NONE;
		}
	}

	/////Set up the sampler/////
	if ( kic_sampling_ ) {
		runtime_assert( close_chain_to_distal_ );

		Size const chainbreak_suite( five_prime_chain_break_res_ );
		runtime_assert( five_prime_chain_break_res_ > 0 );

		pose::PoseOP new_pose = new pose::Pose( pose ); //hard copy
		RNA_KicSamplerOP sampler = new RNA_KicSampler(
				new_pose, moving_suite_, chainbreak_suite );
		if ( !sample_sugar[2] ) {
			sampler->set_base_state( NONE );
			sampler->set_pucker_state( NONE );
		} else {
			sampler->set_base_state( base_state[2] );
			sampler->set_pucker_state( pucker_state[2] );
		}
		sampler->set_verbose( options_->verbose() );
		sampler->set_skip_same_pucker( options_->use_phenix_geo() );
		sampler->set_idealize_coord( options_->use_phenix_geo() );
		sampler->set_extra_epsilon( options_->sampler_extra_epsilon_rotamer() );
		sampler->set_extra_chi(	options_->extra_chi() );
		sampler->set_random( options_->choose_random() );
		sampler->set_fast( options_->integration_test_mode() ); // overrules extra_chi, extra_epsilon; and sets bin size to 40!
		if ( options_->finer_sampling_at_chain_closure() ) sampler->set_bin_size( 10 );
		sampler->init();
		return sampler;
	}

	RNA_SuiteRotamerOP sampler = new RNA_SuiteRotamer( moving_suite_,
			pucker_state[1], pucker_state[2], base_state[1], base_state[2] );
	sampler->set_skip_same_pucker( options_->use_phenix_geo() );
	sampler->set_idealize_coord( options_->use_phenix_geo() );
	sampler->set_sample_nucleoside_lower( sample_sugar[1] );
	sampler->set_sample_nucleoside_upper( sample_sugar[2] );
	sampler->set_fast( options_->integration_test_mode() );
	sampler->set_extra_epsilon( options_->sampler_extra_epsilon_rotamer() );
	sampler->set_extra_beta( options_->sampler_extra_beta_rotamer() );
	sampler->set_extra_chi(	options_->extra_chi() );
	sampler->set_random( options_->choose_random() );
	if ( close_chain_to_distal_ && options_->finer_sampling_at_chain_closure()  )	sampler->set_bin_size( 10 );
	sampler->init();

	return sampler;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener ){
	base_centroid_screener_ = screener;
}

////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_ClassicResidueSampler::bulge_variant_decision( core::pose::Pose & pose, Real const & delta_atr_score ){
	using namespace ObjexxFCL;

	if ( rebuild_bulge_mode_ ) return false; //Hacky want to output sample diverse bulge conformation
	static Real const atr_cutoff_for_bulge( -999999.0 ); //Feb 02, 2012

	runtime_assert ( delta_atr_score <= (  + 0.01 ) );
	runtime_assert ( !is_virtual_base( pose.residue( moving_res_ ) ) );

	bool add_bulge = false;
	if ( options_->allow_bulge_at_chainbreak() ) {
		if ( delta_atr_score >= atr_cutoff_for_bulge ) {
			add_bulge = true;
			count_data_.bulge_at_chain_closure_count++;
			if ( options_->verbose() ){
				TR.Debug << "delta_atr " << delta_atr_score << " passes cutoff for bulge. " << atr_cutoff_for_bulge;
				TR.Debug << "  bulge = " << count_data_.bulge_at_chain_closure_count << "  both = " << count_data_.both_count << " tot = " << count_data_.tot_rotamer_count << std::endl;
			}
		} else {
			add_bulge = false;
			if ( options_->verbose() ) TR.Debug << "delta_atr " << delta_atr_score << " DOES NOT PASS cutoff for bulge " << atr_cutoff_for_bulge << std::endl;
		}
	}

	return add_bulge;

}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::apply_bulge_variant( core::pose::Pose & pose ) const{
	//Note that there is problem in that even after applying virtual_rna_residue, the chain break torsion potential is still scored for the chain_break torsions.
	//The should_score_torsion function in RNA_torsional_potential returns true (indicating that the score should be scored) if it finds a chain_break torsion,
	//even if this torsion contain virtual atoms.. May 4, 2010
	runtime_assert( !is_virtual_base( pose.residue( moving_res_ ) ) );
	apply_virtual_rna_residue_variant_type( pose, moving_res_, true );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::output_count_data(){

	if ( gap_size_ <= 1 ) TR.Debug << " chain_closable_count = " << count_data_.chain_closable_count << std::endl;
	if ( gap_size_ == 0 ){
		TR.Debug << " angle_n = " << count_data_.good_angle_count << " dist_n = " << count_data_.good_distance_count;
		TR.Debug << " chain_break_screening = " << count_data_.chain_break_screening_count << std::endl;
	}
	if ( options_->combine_long_loop_mode() && !close_chain_to_distal_ ) TR.Debug << "res_contact = " << count_data_.residues_contact_screen << " ";

	TR.Debug << "stack = " << count_data_.base_stack_count << " pair = " << count_data_.base_pairing_count;
	TR.Debug << " strict_pair_n = " << count_data_.strict_base_pairing_count;
	TR.Debug << " atr = " << count_data_.good_atr_rotamer_count;
	TR.Debug << " rep = " << count_data_.good_rep_rotamer_count;
	TR.Debug << " both = " << count_data_.both_count;
	TR.Debug << " bulge = " << count_data_.bulge_at_chain_closure_count;
	TR.Debug << " rmsd = " << count_data_.rmsd_count << " tot = " << count_data_.tot_rotamer_count << std::endl;

}


//////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_ClassicResidueSampler::break_early_for_integration_tests() {
	if ( options_->integration_test_mode() && count_data_.both_count >= 100 ) return true;
	if ( options_->integration_test_mode() && count_data_.full_score_count >= 10 ) native_rmsd_screen_ = true;
	if ( options_->integration_test_mode() && count_data_.rmsd_count >= 10 ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){ user_input_VDW_bin_screener_ = user_input_VDW_bin_screener; }

//////////////////////////////////////////////////////////////////
utility::vector1< pose::PoseOP > &
StepWiseRNA_ClassicResidueSampler::pose_list(){
	return pose_list_;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::set_pose_list( utility::vector1< pose::PoseOP > &	pose_list ){
	pose_list_ = pose_list;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ClassicResidueSampler::set_options( StepWiseRNA_ModelerOptionsCOP options ){
	options_ = options;
}

} //suite
} //rna
} //enumerate
} //stepwise
} //protocols
