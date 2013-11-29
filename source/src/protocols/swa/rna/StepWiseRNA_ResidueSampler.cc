// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_ResidueSampler
/// @brief Master sampler for StepWise Assembly of RNA.
/// @detailed
/// @author Parin Sripakdeevong
/// @author Rhiju Das

//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_StandardResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSamplerUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <protocols/swa/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <core/io/pdb/pose_io.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <numeric/NumericTraits.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>
#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>


using namespace core;
using io::pdb::dump_pdb;
using ObjexxFCL::string_of;

//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of RNA
//
//  This carries out standard enumerative sampling of a nucleotide and its
//   connections to the next base OR floating base sampling, which skips
//   bulges (a.k.a., 'dinucleotide move')
//  Prior to these moves, any virtualized sugars at the residue to be sampled
//   or at residues involved in chain closure are instantiated.
//
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_ResidueSampler" ) ;

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_ResidueSampler::StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
	job_parameters_( job_parameters ),
	sfd_( new core::io::silent::SilentFileData ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ), // can be replaced from the outside
	silent_file_( "silent_file.txt" ),
	output_filename_( "data.txt" ),
	rep_cutoff_( 4.0 ),
	bin_size_( 20 ),
	num_pose_kept_( 108 ),
	cluster_rmsd_( 0.5001 ),
	verbose_( false ),
	native_rmsd_screen_( false ),
	native_screen_rmsd_cutoff_( 2.0 ),
	perform_o2prime_pack_( true ),
	use_green_packer_( false ),
	allow_bulge_at_chainbreak_( false ),
	integration_test_mode_( false ), //March 16, 2012
	centroid_screen_( true ),
	allow_base_pair_only_centroid_screen_( false ), //allow for possibility of conformation that base_pair but does not base_stack
	VDW_atr_rep_screen_( true ),
	allow_syn_pyrimidine_( false ), //New option Nov 15, 2010
	distinguish_pucker_( true ),
	finer_sampling_at_chain_closure_( false ), //New option Jun 10 2010
	PBP_clustering_at_chain_closure_( false ), //New option Aug 15 2010
	reinitialize_CCD_torsions_( false ), //New option Aug 15 2010 //Reinitialize_CCD_torsion to zero before every CCD chain closure
	extra_epsilon_rotamer_( false ), //New option Aug 30, 2010
	extra_beta_rotamer_( false ), //New option Aug 30, 2010
	extra_chi_( false ),
	sample_both_sugar_base_rotamer_( false ), //New option Nov 12, 2010 (mainly for square_RNA)
	include_torsion_value_in_tag_( false ), //For checking if the extra rotamer are important
	combine_long_loop_mode_( false ), //in this mode, the moving_residues must contact the last residue built from the other side.
	do_not_sample_multiple_virtual_sugar_( false ), //Nov 13, 2010, optimize the chain closure step speed
	sample_ONLY_multiple_virtual_sugar_( false ), //Nov 13, 2010, optimize the chain closure step speed
	assert_no_virt_sugar_sampling_( false ), //July 28 2011
	output_pdb_( false ), //Sept 24, 2011
	choose_random_( false ), // Rhiju, Jul 2013
	num_random_samples_( 1 ),
	force_centroid_interaction_( false ),  // Rhiju, Jul 2013
	use_phenix_geo_( false ),
	kic_sampling_( false )
{
	set_native_pose( job_parameters_->working_native_pose() );
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseRNA_ResidueSampler::~StepWiseRNA_ResidueSampler()
{}

/////////////////////
std::string
StepWiseRNA_ResidueSampler::get_name() const {
	return "StepWiseRNA_ResidueSampler";
}

////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::apply( core::pose::Pose & pose ) {
	using namespace ObjexxFCL;

	output_title_text( "Enter StepWiseRNA_ResidueSampler::apply", TR.Debug );
	clock_t const time_start( clock() );
	output_options();
	Pose const pose_save = pose; pose = pose_save; //this recopy is actually useful for triggering graphics.

	instantiate_any_virtual_sugars( pose ); // saves work in "anchor_sugar_modeling_", etc. FloatingBaseJobParameter.

	if ( job_parameters_->floating_base() ){
		floating_base_sampling( pose );
	} else{
		standard_sampling_WRAPPER( pose );
	}

	TR.Debug << "Total time in StepWiseRNA_ResidueSampler::apply " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
	output_title_text( "Exit StepWiseRNA_ResidueSampler::apply", TR.Debug );
	pose = pose_save;
}


// This can/should move to its own class...

///////////////////////////////////////////////////////////////////////////////////
//////////////////////Build previously virtualize sugar/////////////////////
// A virtualized sugar occurs when a previous move was a 'floating base'
// step, which only samples euler angles of a base, but virtualized the
// attached sugar and any residues connecting that nucleotide to the
//'instantiated' body of the RNA.

// There are potentially four different virtualized sugar positions.
// (1) Anchor_sugar (most common case) :
//      This is the sugar of the nucleotide that was built immediately
//      before the current/moving nucleotide along the chain. Corresponds
//      to sugar of residue (moving_res - 1) if current step is built in
//      the forward (3') direction. Likewise, corresponds to sugar of
//      residue (moving_res + 1) if current step built in the backward
//      (5') direction.
//
// (2) Current_sugar (rare) :
//      The sugar of the current/moving nucleotide. This sugar can
//      be virtual in the situation where the current step is a step to
//      combine two moving_elements that were previously built with SWA. It is
//      possible that the current nucleotide (in the moving moving_element) was
//      previously built with a 'floating base' step and therefore will
//      contain a virtual sugar.
//
// (3) Five_prime_chain_break_sugar :
//      Occurs only in the context where the current step is a
//      chain-closure (closing loop) step. This is the sugar of the
//      nucleotide immediately 5' of the chain-closure phosphate group.
//      Five_prime_chain_break_sugar can be identical to current_sugar
//      in certain situations. The code below contains check statements to
//      prevent double-counting in these situations.
//
// (4) Three_prime_chain_break_sugar :
//      Occurs only in the context where the current step is a
//      chain-closure (closing loop) step. This is the sugar of the
//      nucleotide immediately 3' of the chain-closure phosphate group.
//      Three_prime_chain_break_sugar can be identical to current_sugar
//      in certain situations. The code below contains check statements to
//      prevent double-counting in these situations.
void
StepWiseRNA_ResidueSampler::instantiate_any_virtual_sugars( pose::Pose & pose ){

	output_title_text( "Build previously virtualize sugar", TR.Debug );
	bool const is_anchor_sugar_virt = is_anchor_sugar_virtual( pose );
	bool const is_curr_sugar_virt = !job_parameters_->floating_base() && is_current_sugar_virtual( pose ); // occurs less frequently -- typically only when connecting two pieces.
	bool const is_five_prime_CB_sugar_virt = is_five_prime_chain_break_sugar_virtual( pose ); // may occur when closing loop
	bool const is_three_prime_CB_sugar_virt = is_three_prime_chain_break_sugar_virtual( pose ); // may occur when closing loop
	Size const num_nucleotides =  job_parameters_->working_moving_res_list().size();
	Size const gap_size = job_parameters_->gap_size();

	output_boolean( " is_anchor_sugar_virt = ", is_anchor_sugar_virt, TR.Debug );
	output_boolean( " is_current_sugar_virt = ", is_curr_sugar_virt, TR.Debug );
	output_boolean( " is_five_prime_chain_break_sugar_virt = ", is_five_prime_CB_sugar_virt, TR.Debug );
	output_boolean( " is_three_prime_chain_break_sugar_virt = ", is_three_prime_CB_sugar_virt, TR.Debug );
	TR.Debug << std::endl;

	initialize_scorefunctions();

	///Nov 13, 2010
	Size num_virtual_sugar = 0;
	if ( is_anchor_sugar_virt ) num_virtual_sugar++;
	if ( is_curr_sugar_virt ) num_virtual_sugar++;
	if ( is_five_prime_CB_sugar_virt ) num_virtual_sugar++;
	if ( is_three_prime_CB_sugar_virt ) num_virtual_sugar++;

	TR.Debug << "num_virtual_sugar = " << num_virtual_sugar << std::endl;

	bool const floating_base = job_parameters_->floating_base();

	// rd2013 -- need comments explaining these checks
	if ( assert_no_virt_sugar_sampling_ /*run checks*/ ){
		if ( floating_base && num_nucleotides == 2 ){ //Hacky..ok the only acception right now is in floating_base + dinucleotide mode.
			if ( num_virtual_sugar > 1 ) utility_exit_with_message( "assert_no_virt_sugar_sampling_ == true and floating_base, but num_virtual_sugar > 1" );
			if ( num_virtual_sugar == 1 ){
				if ( is_anchor_sugar_virt == false ) {
					utility_exit_with_message( "assert_no_virt_sugar_sampling_ == true and floating_base and num_virtual_sugar == 1 BUT is_anchor_sugar_virt == false! )" );
				}
			}
		} else{
			if ( num_virtual_sugar != 0 ) utility_exit_with_message( "assert_no_virt_sugar_sampling_ == true but num_virtual_sugar != 0" );
		}
	}

	// rd2013 -- I don't yet understand what logic guarantees that these assertions are met.
	if ( is_curr_sugar_virt ){ //Consistency test.
		//Right now, the only possiblility for is_curr_sugar_virt==true is when combining two silent_files moving_element at chain-break
		if ( gap_size != 0 ) utility_exit_with_message( "is_curr_sugar_virt == true but gap_size != 0 !!" );
		utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();
		if ( working_moving_partition_pos.size() <= 1 ) utility_exit_with_message( "is_curr_sugar_virt == true but working_moving_partition_pos.size() <= 1" );
	}
	if ( gap_size != 0 && num_virtual_sugar > 1 ) utility_exit_with_message( "gap_size != 0 but num_virtual_sugar > 1" ); //Obsolete!

	if ( floating_base ){
		if ( is_five_prime_CB_sugar_virt ){ //This is rare since floating_base sampling is not often used at chain-closure step!
			TR.Debug << "WARNING: floating_base and is_five_prime_CB_sugar_virt case. Code not implemented yet, early return!" << std::endl;
			return;
		}

		if ( is_three_prime_CB_sugar_virt ){ //This is rare since floating_base sampling is not often used at chain-closure step!
			TR.Debug << "WARNING: floating_base and is_three_prime_CB_sugar_virt case. Code not implemented yet, early return!" << std::endl;
			return;
		}

		if ( is_curr_sugar_virt ){
			utility_exit_with_message( "floating_base and is_curr_sugar_virt case. Code not implemented yet!" );
		}
	}

	if ( do_not_sample_multiple_virtual_sugar_ && sample_ONLY_multiple_virtual_sugar_ ){
		utility_exit_with_message( "do_not_sample_multiple_virtual_sugar_ == true && sample_ONLY_multiple_virtual_sugar_ == true" );
	}

	// this was just for testing.
	if ( do_not_sample_multiple_virtual_sugar_ ){
		if ( num_virtual_sugar > 1 ) return;
	}

	if ( sample_ONLY_multiple_virtual_sugar_ ){
		if ( gap_size != 0 ) utility_exit_with_message( "sample_ONLY_multiple_virtual_sugar_ == true but gap_size != 0" );
		if ( num_virtual_sugar <= 1 ) return;
	}

	bool const is_prepend(  job_parameters_->is_prepend() );
	Size const moving_res(  job_parameters_->working_moving_res() );
	Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
	Size const three_prime_chain_break_res = five_prime_chain_break_res + 1;

	// These parameters will be used for sampling virtual sugars and associated bulge. Reset them here.
	anchor_sugar_modeling_ = SugarModeling();
	curr_sugar_modeling_ = SugarModeling();
	five_prime_CB_sugar_modeling_ = SugarModeling();
	three_prime_CB_sugar_modeling_ = SugarModeling();

	// In following, PDL means pose_data_list.
	if ( is_anchor_sugar_virt ){
		TR.Debug << "anchor_sugar floating_base_chain_closure" << std::endl;

		Size const anchor_moving_res = ( is_prepend ) ? ( moving_res + num_nucleotides ) : ( moving_res - num_nucleotides );
		Size const anchor_ref_res    = ( is_prepend ) ? ( moving_res + ( num_nucleotides + 2 ) ) : ( moving_res - ( num_nucleotides + 2 ) );

		anchor_sugar_modeling_ = SugarModeling( anchor_moving_res, anchor_ref_res );
		anchor_sugar_modeling_.set_base_and_pucker_state( pose, job_parameters_ );

		// do the sampling!
		anchor_sugar_modeling_.PDL = anchor_floating_base_chain_closure( pose, anchor_sugar_modeling_, "anchor" );

		if ( anchor_sugar_modeling_.PDL.size() == 0 ){
			TR.Debug << "is_anchor_sugar_virt == True but anchor_sugar_modeling_.PDL.size() == 0!" << std::endl;
			return;
		}
	}

	if ( is_curr_sugar_virt ){ //June 12, 2011
		TR.Debug << "current_sugar floating_base_chain_closure" << std::endl;

		Size const curr_ref_res = ( is_prepend ) ? ( moving_res - 2 ) : ( moving_res + 2 );

		curr_sugar_modeling_ = SugarModeling( moving_res, curr_ref_res );
		curr_sugar_modeling_.set_base_and_pucker_state( pose, job_parameters_ );

		// do the sampling!
		curr_sugar_modeling_.PDL = anchor_floating_base_chain_closure( pose, curr_sugar_modeling_, "current" );

		if ( curr_sugar_modeling_.PDL.size() == 0 ){
			TR.Debug << "is_curr_sugar_virt == True but curr_sugar_modeling_.PDL.size() == 0!" << std::endl;
			return;
		}
	}

	if ( is_five_prime_CB_sugar_virt ){
		TR.Debug << "five_prime_CB_sugar floating_base_chain_closure" << std::endl;

		five_prime_CB_sugar_modeling_ = SugarModeling( five_prime_chain_break_res, five_prime_chain_break_res - 2 );
		five_prime_CB_sugar_modeling_.set_base_and_pucker_state( pose, job_parameters_ );

		// do the sampling!
		five_prime_CB_sugar_modeling_.PDL = anchor_floating_base_chain_closure( pose, five_prime_CB_sugar_modeling_, "five_prime_CB" );

		if ( five_prime_CB_sugar_modeling_.PDL.size() == 0 ) {
			TR.Debug << "is_five_prime_CB_sugar_virt == True but five_prime_CB_sugar_modeling_.PDL.size() == 0!" << std::endl;
			return;
		}
	}

	if ( is_three_prime_CB_sugar_virt ){
		TR.Debug << "three_prime_CB_sugar floating_base_chain_closure:" << std::endl;

		three_prime_CB_sugar_modeling_ = SugarModeling( three_prime_chain_break_res, three_prime_chain_break_res + 2 );
		three_prime_CB_sugar_modeling_.set_base_and_pucker_state( pose, job_parameters_ );

		// do the sampling!
		three_prime_CB_sugar_modeling_.PDL = anchor_floating_base_chain_closure( pose, three_prime_CB_sugar_modeling_, "three_prime_CB" );

		if ( three_prime_CB_sugar_modeling_.PDL.size() == 0 ){
			TR.Debug << "is_three_prime_CB_sugar_virt == True but three_prime_CB_sugar_modeling_.PDL.size() == 0!" << std::endl;
			return;
		}
	}

	if (	num_virtual_sugar > 0 ){
		TR.Debug << "Contain_virtual_sugar == true" << std::endl;
	} else{
		TR.Debug << "Contain_virtual_sugar == false" << std::endl;
	}

	/////////////Sort the pose_data_list by score..should be determined according to torsional potential score.	/////////////////////////
	std::sort( anchor_sugar_modeling_.PDL.begin(), anchor_sugar_modeling_.PDL.end(), sort_pose_by_score );
	std::sort( curr_sugar_modeling_.PDL.begin(), curr_sugar_modeling_.PDL.end(), sort_pose_by_score );
	std::sort( five_prime_CB_sugar_modeling_.PDL.begin(),  five_prime_CB_sugar_modeling_.PDL.end(), sort_pose_by_score );
	std::sort( three_prime_CB_sugar_modeling_.PDL.begin(), three_prime_CB_sugar_modeling_.PDL.end(), sort_pose_by_score );
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::floating_base_sampling( pose::Pose & pose ){

	StepWiseRNA_FloatingBaseSampler floating_base_sampler( job_parameters_ );
	floating_base_sampler.set_base_centroid_screener( base_centroid_screener_ );
	floating_base_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );
	floating_base_sampler.set_silent_file ( silent_file_ );
	floating_base_sampler.set_scorefxn ( scorefxn_ );
	floating_base_sampler.set_num_pose_kept ( num_pose_kept_ );
	floating_base_sampler.set_native_rmsd_screen ( native_rmsd_screen_ );
	floating_base_sampler.set_native_screen_rmsd_cutoff ( native_screen_rmsd_cutoff_ );
	floating_base_sampler.set_perform_o2prime_pack ( perform_o2prime_pack_ );
	floating_base_sampler.set_verbose ( verbose_ );
	floating_base_sampler.set_cluster_rmsd (	cluster_rmsd_	);
	floating_base_sampler.set_distinguish_pucker ( distinguish_pucker_ );
	floating_base_sampler.set_PBP_clustering_at_chain_closure ( PBP_clustering_at_chain_closure_ );
	floating_base_sampler.set_use_phenix_geo ( use_phenix_geo_  );
	floating_base_sampler.set_centroid_screen ( centroid_screen_ );
	floating_base_sampler.set_VDW_atr_rep_screen ( VDW_atr_rep_screen_ );
	floating_base_sampler.set_integration_test_mode( integration_test_mode_ ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	floating_base_sampler.set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );
	floating_base_sampler.set_allow_base_pair_only_centroid_screen( allow_base_pair_only_centroid_screen_ );
	runtime_assert( !use_green_packer_ );
	runtime_assert( !combine_long_loop_mode_ );

	floating_base_sampler.apply( pose );
	pose_data_list_ = floating_base_sampler.get_pose_data_list();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
StepWiseRNA_ResidueSampler::standard_sampling_WRAPPER( core::pose::Pose & pose ){

	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::scoring;

	// extraneous?
	pose::Pose const pose_save = pose;

	StepWiseRNA_StandardResidueSampler standard_residue_sampler( job_parameters_ );
	standard_residue_sampler.set_pose_data_list( pose_data_list_ ); // allows for accumulation of poses, if desired.
	standard_residue_sampler.set_base_centroid_screener( base_centroid_screener_ );
	standard_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );
	standard_residue_sampler.set_silent_file ( silent_file_ );
	standard_residue_sampler.set_scorefxn ( scorefxn_ );
	standard_residue_sampler.set_num_pose_kept ( num_pose_kept_ );
	standard_residue_sampler.set_native_rmsd_screen ( native_rmsd_screen_ );
	standard_residue_sampler.set_native_screen_rmsd_cutoff ( native_screen_rmsd_cutoff_ );
	standard_residue_sampler.set_perform_o2prime_pack ( perform_o2prime_pack_ );
	standard_residue_sampler.set_use_green_packer ( use_green_packer_ );
	standard_residue_sampler.set_verbose ( verbose_ );
	standard_residue_sampler.set_cluster_rmsd (	cluster_rmsd_	);
	standard_residue_sampler.set_distinguish_pucker ( distinguish_pucker_ );
	standard_residue_sampler.set_finer_sampling_at_chain_closure ( finer_sampling_at_chain_closure_ );
	standard_residue_sampler.set_PBP_clustering_at_chain_closure ( PBP_clustering_at_chain_closure_ );
	standard_residue_sampler.set_allow_syn_pyrimidine( allow_syn_pyrimidine_ );
	standard_residue_sampler.set_extra_chi( extra_chi_ );
	standard_residue_sampler.set_use_phenix_geo ( use_phenix_geo_  );
	standard_residue_sampler.set_kic_sampling( kic_sampling_ );
	standard_residue_sampler.set_centroid_screen ( centroid_screen_ );
	standard_residue_sampler.set_VDW_atr_rep_screen ( VDW_atr_rep_screen_ );
	standard_residue_sampler.set_force_centroid_interaction ( force_centroid_interaction_ );
	standard_residue_sampler.set_integration_test_mode( integration_test_mode_ ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	standard_residue_sampler.set_allow_bulge_at_chainbreak( allow_bulge_at_chainbreak_ );
	standard_residue_sampler.set_parin_favorite_output( parin_favorite_output_ );
	standard_residue_sampler.set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );
	standard_residue_sampler.set_extra_epsilon_rotamer( extra_epsilon_rotamer_ );
	standard_residue_sampler.set_extra_beta_rotamer( extra_beta_rotamer_ );
	standard_residue_sampler.set_sample_both_sugar_base_rotamer( sample_both_sugar_base_rotamer_ ); //Nov 12, 2010
	standard_residue_sampler.set_include_torsion_value_in_tag( include_torsion_value_in_tag_ );
	standard_residue_sampler.set_combine_long_loop_mode( combine_long_loop_mode_ );
	standard_residue_sampler.set_allow_base_pair_only_centroid_screen( allow_base_pair_only_centroid_screen_ );

	if ( !sampling_sugar() ){
		standard_residue_sampler.apply( pose );
	} else{ //Case where have to sample virtual sugar...
		utility::vector1< PoseOP > starting_pose_data_list;
		prepare_from_prior_sampled_sugar_jobs( pose, starting_pose_data_list );
		for ( Size n = 1; n <= starting_pose_data_list.size(); n++ ){
			pose = ( *starting_pose_data_list[n] ); //set viewer_pose;
			standard_residue_sampler.set_extra_tag( tag_from_pose( *starting_pose_data_list[n] )  );
			standard_residue_sampler.apply( pose );
		}
	}

	pose_data_list_ = standard_residue_sampler.pose_data_list();
	pose = pose_save;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
StepWiseRNA_ResidueSampler::set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){ user_input_VDW_bin_screener_ = user_input_VDW_bin_screener; }



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP >
StepWiseRNA_ResidueSampler::anchor_floating_base_chain_closure( pose::Pose & viewer_pose, SugarModeling const & sugar_modeling, std::string const name ){

	return sample_virtual_sugar_and_bulge_and_close_chain(
			viewer_pose, sugar_modeling, name,	scorefxn_, sampling_scorefxn_,
			atr_rep_screening_scorefxn_, chainbreak_scorefxn_, job_parameters_,
			integration_test_mode_, use_phenix_geo_ );
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
StepWiseRNA_ResidueSampler::is_anchor_sugar_virtual( core::pose::Pose const & pose ) const {

	//Check if anchor sugar is virtual, if virtual then need to sample it.
	bool const is_prepend(  job_parameters_->is_prepend() ); // If true, moving_suite+1 is fixed. Otherwise, moving_suite is fixed.
	Size const moving_res(  job_parameters_->working_moving_res() ); // Corresponds to user input.
	Size const num_nucleotides(  job_parameters_->working_moving_res_list().size() );
	Size const anchor_moving_res = ( is_prepend ) ? ( moving_res + num_nucleotides ) : ( moving_res - num_nucleotides );
	Size const anchor_bulge_res = ( is_prepend ) ? ( moving_res + ( num_nucleotides + 1 ) ) : ( moving_res - ( num_nucleotides + 1 ) );

	return is_sugar_virtual(  pose, anchor_moving_res, anchor_bulge_res );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////New June 12, 2011/////
bool
StepWiseRNA_ResidueSampler::is_current_sugar_virtual( core::pose::Pose const & pose ) const {

	// Check if curr sugar is virtual. If virtual then need to sample it.
	// This occur when combining two moving_element and the moving_res
	// in the moving_element was built with a dinucleotide move.
	bool const is_prepend( job_parameters_->is_prepend() );
	Size const moving_res( job_parameters_->working_moving_res() );
	Size const virtual_sugar_res = moving_res;
	Size const bulge_res = ( is_prepend ) ? ( moving_res - 1 ) : ( moving_res + 1 );

	return is_sugar_virtual(  pose, virtual_sugar_res, bulge_res );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
StepWiseRNA_ResidueSampler::is_five_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) const {

	Size const moving_res( job_parameters_->working_moving_res() );
	Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
	Size const gap_size( job_parameters_->gap_size() );

	if ( gap_size != 0 ) return false;

	// Make sure to not over count number of virtual_sugar-virtual_bulge
	// pairs to be build!
	if ( moving_res == five_prime_chain_break_res ) return false;

	Size const five_prime_CB_bulge_res = ( five_prime_chain_break_res - 1 );

	bool sugar_is_virtual = is_sugar_virtual( pose, five_prime_chain_break_res, five_prime_CB_bulge_res );

	return sugar_is_virtual;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
StepWiseRNA_ResidueSampler::is_three_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) const {

	Size const moving_res(  job_parameters_->working_moving_res() );
	Size const three_prime_chain_break_res = job_parameters_->five_prime_chain_break_res() + 1;
	Size const gap_size( job_parameters_->gap_size() );

	if ( gap_size != 0 ) return false;

	// Make sure to not over count number of virtual_sugar-virtual_bulge
	// pairs to be build!
	if ( moving_res == three_prime_chain_break_res ) return false;

	Size const three_prime_CB_bulge_res = ( three_prime_chain_break_res + 1 );

	bool sugar_is_virtual = is_sugar_virtual(  pose, three_prime_chain_break_res, three_prime_CB_bulge_res );

	return sugar_is_virtual;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::initialize_scorefunctions(){

	initialize_common_scorefxns( scorefxn_, sampling_scorefxn_, atr_rep_screening_scorefxn_, chainbreak_scorefxn_, o2prime_pack_scorefxn_ );

}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseRNA_ResidueSampler::get_pose_data_list(){
	return pose_data_list_;
}

///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_silent_file( std::string const & silent_file ){
	silent_file_ = silent_file;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_num_pose_kept( core::Size const & num_pose_kept ){
	num_pose_kept_ = num_pose_kept ;
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_native_rmsd_screen( bool const & setting ){
	native_rmsd_screen_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_native_screen_rmsd_cutoff( core::Real const & setting ){
	native_screen_rmsd_cutoff_ = setting;
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_integration_test_mode( bool const & setting ){
	integration_test_mode_ = setting;
	if ( integration_test_mode_ ){
		num_pose_kept_ = 2;
		native_rmsd_screen_ = false; // will be switched to true mid-way through sampling.
		native_screen_rmsd_cutoff_ = 1.0;
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_verbose( bool const & setting ){
	verbose_ = setting;
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_perform_o2prime_pack( bool const & setting ){
	perform_o2prime_pack_ = setting;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_allow_bulge_at_chainbreak( bool const & setting ){
	allow_bulge_at_chainbreak_ = setting;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_output_filename( std::string const & output_filename ){
	output_filename_ = output_filename;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
core::io::silent::SilentFileDataOP &
StepWiseRNA_ResidueSampler::silent_file_data(){
	return sfd_;
}


//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::output_pose_data_list( std::string const final_sampler_output_silent_file ) const{
	using namespace core::io::silent;

	if ( verbose_ == false ){ //consistency check Apr 3, 2010
		utility_exit_with_message( "verbose_ == false, but StepWiseRNA_ResidueSampler::output_pose_data_list is still called?!" );
	}

	SilentFileData silent_file_data;

	for ( Size n = 1; n <= pose_data_list_.size(); n++ ) {
		output_data( silent_file_data, final_sampler_output_silent_file, tag_from_pose( *pose_data_list_[n] ), false, *( pose_data_list_[n] ), get_native_pose(), job_parameters_ );
	}

}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener ){
	base_centroid_screener_ = screener;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_cluster_rmsd( Real const & setting ){
	cluster_rmsd_ = setting;
	TR.Debug << "Set cluster_rmsd to " << cluster_rmsd_ << std::endl;
}


//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::output_options(){
	//output screen options
	TR.Debug << "--------SCREEN OPTIONS---------- " << std::endl;
	output_boolean( "integration_test_mode_ = ", integration_test_mode_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "native_rmsd_screen = ", native_rmsd_screen_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << "native_screen_rmsd_cutoff = " << native_screen_rmsd_cutoff_ << std::endl;
	output_boolean( "perform_o2prime_pack = ", perform_o2prime_pack_, TR.Debug ); TR.Debug << std::endl;
	output_seq_num_list( "working_moving_partition_pos = ", job_parameters_->working_moving_partition_pos(), TR.Debug );
	output_boolean( "centroid_screen = ", centroid_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "allow_base_pair_only_centroid_screen = ", allow_base_pair_only_centroid_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "VDW_atr_rep_screen = ", VDW_atr_rep_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "sample_both_sugar_base_rotamer_ = ", sample_both_sugar_base_rotamer_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "do_not_sample_multiple_virtual_sugar_ = ", do_not_sample_multiple_virtual_sugar_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "sample_ONLY_multiple_virtual_sugar_ = ", sample_ONLY_multiple_virtual_sugar_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "assert_no_virt_sugar_sampling_ = ", assert_no_virt_sugar_sampling_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "distinguish_pucker_ ", distinguish_pucker_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << "--------------------------------" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::check_res_not_bulged(){
	for ( Size n = 1; n <= pose_data_list_.size(); n++ ){
		if ( ( *pose_data_list_[n] ).residue( job_parameters_->working_moving_res() ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
			utility_exit_with_message( "working_moving_res: " + string_of( job_parameters_->working_moving_res() ) + " of pose " + string_of( n ) + " is a virtual res!"  );
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_ResidueSampler::sampling_sugar() const{
	return ( anchor_sugar_modeling_.sample_sugar || curr_sugar_modeling_.sample_sugar || five_prime_CB_sugar_modeling_.sample_sugar || three_prime_CB_sugar_modeling_.sample_sugar );
}


////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_ResidueSampler::prepare_from_prior_sampled_sugar_jobs( pose::Pose const & pose,
																																	 utility::vector1< PoseOP > & starting_pose_data_list ) {

	if ( anchor_sugar_modeling_.PDL.size() == 0 && curr_sugar_modeling_.PDL.size() == 0 && five_prime_CB_sugar_modeling_.PDL.size() == 0 && three_prime_CB_sugar_modeling_.PDL.size() == 0 ){
		utility_exit_with_message( "pose_data_list is empty for all 4 possible virtual sugar!" );
	}

	Size count = 0;

	for ( Size anchor_sugar_ID = 1; anchor_sugar_ID <= anchor_sugar_modeling_.PDL.size() || anchor_sugar_ID == 1; anchor_sugar_ID++ ){
		for ( Size curr_sugar_ID = 1; curr_sugar_ID <= curr_sugar_modeling_.PDL.size() || curr_sugar_ID == 1; curr_sugar_ID++ ){
			for ( Size five_prime_CB_sugar_ID = 1; five_prime_CB_sugar_ID <= five_prime_CB_sugar_modeling_.PDL.size() || five_prime_CB_sugar_ID == 1; five_prime_CB_sugar_ID++ ){
				for ( Size three_prime_CB_sugar_ID = 1; three_prime_CB_sugar_ID <= three_prime_CB_sugar_modeling_.PDL.size() || three_prime_CB_sugar_ID == 1; three_prime_CB_sugar_ID++ ){

					count++;

					PoseOP start_pose_data;

					start_pose_data = new pose::Pose;
					( *start_pose_data ) = pose;
					pose::Pose & start_pose = ( *start_pose_data );
					tag_into_pose( start_pose, "" );

					if ( anchor_sugar_modeling_.PDL.size() > 0 ) {
						tag_into_pose( start_pose,   tag_from_pose(start_pose) + tag_from_pose( *anchor_sugar_modeling_.PDL[anchor_sugar_ID] ) );
						copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, start_pose, ( *anchor_sugar_modeling_.PDL[anchor_sugar_ID] ) );
					} else{
						tag_into_pose( start_pose, tag_from_pose( start_pose ) + "_null" );
					}

					if ( curr_sugar_modeling_.PDL.size() > 0 ) {
						tag_into_pose( start_pose,   tag_from_pose(start_pose) + tag_from_pose( *curr_sugar_modeling_.PDL[curr_sugar_ID] ) );
						copy_bulge_res_and_sugar_torsion( curr_sugar_modeling_, start_pose, ( *curr_sugar_modeling_.PDL[curr_sugar_ID] ) );
					} else{
						tag_into_pose( start_pose, tag_from_pose( start_pose ) + "_null" );
					}

					if ( five_prime_CB_sugar_modeling_.PDL.size() > 0 ){
						tag_into_pose( start_pose,   tag_from_pose(start_pose) + tag_from_pose( *five_prime_CB_sugar_modeling_.PDL[five_prime_CB_sugar_ID] ) );
						copy_bulge_res_and_sugar_torsion( five_prime_CB_sugar_modeling_, start_pose, ( *five_prime_CB_sugar_modeling_.PDL[five_prime_CB_sugar_ID] ) );
					} else{
						tag_into_pose( start_pose, tag_from_pose( start_pose ) + "_null" );
					}

					if ( three_prime_CB_sugar_modeling_.PDL.size() > 0 ){
						tag_into_pose( start_pose,   tag_from_pose(start_pose) + tag_from_pose( *three_prime_CB_sugar_modeling_.PDL[three_prime_CB_sugar_ID] ) );
						copy_bulge_res_and_sugar_torsion( three_prime_CB_sugar_modeling_, start_pose, ( *three_prime_CB_sugar_modeling_.PDL[three_prime_CB_sugar_ID] ) );
					} else{
						tag_into_pose( start_pose, tag_from_pose( start_pose ) + "_null" );
					}

					starting_pose_data_list.push_back( start_pose_data );
				}
			}
		}
	}

	utility::vector1< SugarModeling > sampled_sugar_modeling_list;
	if ( anchor_sugar_modeling_.PDL.size() > 0 ) sampled_sugar_modeling_list.push_back( anchor_sugar_modeling_ );
	if ( curr_sugar_modeling_.PDL.size() > 0 ) sampled_sugar_modeling_list.push_back( curr_sugar_modeling_ );
	if ( five_prime_CB_sugar_modeling_.PDL.size() > 0 ) sampled_sugar_modeling_list.push_back( five_prime_CB_sugar_modeling_ );
	if ( three_prime_CB_sugar_modeling_.PDL.size() > 0 ) sampled_sugar_modeling_list.push_back( three_prime_CB_sugar_modeling_ );


	///////Ok, finally have to remove clashes that may arise due to the fact that the floating base sugar sampling and minimization were done individually of each other///
	Pose pose_copy = pose;
	minimize_all_sampled_floating_bases( pose_copy, sampled_sugar_modeling_list, starting_pose_data_list, sampling_scorefxn_, job_parameters_, true /*virtual_sugar_is_from_prior_step*/ );

	TR.Debug << "starting_pose_data_list.size() = " << starting_pose_data_list.size() << std::endl;

}

}
}
}
