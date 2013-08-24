// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_ResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Parin Sripakdeevong
/// @author Rhiju Das

//////////////////////////////////
//Test_comment
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSamplerUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.hh>

#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <core/chemical/rna/RNA_Util.hh>
#include <core/pose/rna/RNA_Util.hh>


//#include <protocols/swa/rna/StepWiseRNA_Dinucleotide_Sampler_Util.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
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

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>

#include <protocols/rotamer_sampler/rna/RNA_SuiteRotamer.hh>

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

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of RNA
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.rna.stepwise_rna_residue_sampler" ) ;

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
		bin_size_( 20 ),
		rep_cutoff_( 4.0 ),
		num_pose_kept_( 108 ),
		multiplier_( 2 ), //Sort and cluster poses when the number of pose is pose_data_list exceed multiplier*num_pose_kept,
		cluster_rmsd_( 0.5001 ),
		verbose_( false ),
		native_rmsd_screen_( false ),
		native_screen_rmsd_cutoff_( 2.0 ),
		perform_o2star_pack_( true ),
		use_green_packer_( false ),
		allow_bulge_at_chainbreak_( false ),
		fast_( false ),
		medium_fast_( false ),
		integration_test_mode_( false ), //March 16, 2012
		floating_base_( false ), //July 02, 2012
		centroid_screen_( true ),
		allow_base_pair_only_centroid_screen_( false ), //allow for possibility of conformation that base_pair but does not base_stack
		VDW_atr_rep_screen_( true ),
		include_syn_chi_( true ),
		allow_syn_pyrimidine_( false ), //New option Nov 15, 2010
		distinguish_pucker_( true ),
		build_pose_from_scratch_( false ), //July 03, 2012. This was previously uninitialized in floating_base_sampling mode! lead to early return from cluster_pose_data_list() function on some compiler/machine! [Specifically on the Biox2-cluster RUNS after converted to qsub (~April 2012)
		current_score_cutoff_( 999999.9 ), //Feb 02, 2012
		//current_score_cutoff_(99999999999.9999), //New option May 12, 2010, Feb 02, 2012; This might lead to server-test error at R47200
		finer_sampling_at_chain_closure_( false ), //New option Jun 10 2010
		PBP_clustering_at_chain_closure_( false ), //New option Aug 15 2010
		reinitialize_CCD_torsions_( false ), //New option Aug 15 2010 //Reinitialize_CCD_torsion to zero before every CCD chain closure
		extra_epsilon_rotamer_( false ), //New option Aug 30, 2010
		extra_beta_rotamer_( false ), //New option Aug 30, 2010
		extra_chi_( false ),
		sample_both_sugar_base_rotamer_( false ), //New option Nov 12, 2010 (mainly for square_RNA)
		include_torsion_value_in_tag_( false ), //For checking if the extra rotamer are important
		rebuild_bulge_mode_( false ),
		debug_epsilon_south_sugar_mode_( false ),
		exclude_alpha_beta_gamma_sampling_( false ),
		combine_long_loop_mode_( false ), //in this mode, the moving_residues must contact the last residue built from the other side.
		do_not_sample_multiple_virtual_sugar_( false ), //Nov 13, 2010, optimize the chain closure step speed
		sample_ONLY_multiple_virtual_sugar_( false ), //Nov 13, 2010, optimize the chain closure step speed
		assert_no_virt_ribose_sampling_( false ), //July 28 2011
		output_pdb_( false ), //Sept 24, 2011
		choose_random_( false ), // Rhiju, Jul 2013
		num_random_samples_( 1 ),
		force_centroid_interaction_( false )  // Rhiju, Jul 2013
  {
		set_native_pose( job_parameters_->working_native_pose() );

		////////////////Parin Feb 28, 2010////////////////////////////////////////////////
		utility::vector1 < core::Size > const & rmsd_res_list = job_parameters_->rmsd_res_list();
		working_rmsd_res_ = apply_full_to_sub_mapping( rmsd_res_list, job_parameters );

		std::map< core::Size, bool > const & Is_prepend_map = job_parameters_->Is_prepend_map();

		Output_is_prepend_map( "Is_prepend_map = ", Is_prepend_map, job_parameters_->full_sequence().size(), TR.Debug, 30 );
		Output_seq_num_list( "rmsd_res = ", rmsd_res_list, TR.Debug, 30 );
		Output_seq_num_list( "working_rmsd_res = ", working_rmsd_res_, TR.Debug, 30 );
		////////////////////////////////////////////////////////////////////////////////

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


	//////////////////////////////////////////////////////////////////////////
	bool
	sort_criteria( pose_data_struct2  pose_data_1, pose_data_struct2 pose_data_2 ) {  //This function used to be call sort_criteria2
		return ( pose_data_1.score < pose_data_2.score );
	}
	////////////////////////////////////////////////////////////////////////////


	void
	StepWiseRNA_ResidueSampler::apply( core::pose::Pose & pose ) {

		using namespace ObjexxFCL;

		Output_title_text( "Enter StepWiseRNA_ResidueSampler::apply", TR.Debug );

		clock_t const time_start( clock() );

		//output screen options
		TR.Debug << "--------SCREEN OPTIONS---------- " << std::endl;
		Output_boolean( "fast_ = ", fast_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "medium_fast_ = ", medium_fast_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "integration_test_mode_ = ", integration_test_mode_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "native_rmsd_screen = ", native_rmsd_screen_, TR.Debug ); TR.Debug << std::endl;
		TR.Debug << "native_screen_rmsd_cutoff = " << native_screen_rmsd_cutoff_ << std::endl;
		Output_boolean( "perform_o2star_pack = ", perform_o2star_pack_, TR.Debug ); TR.Debug << std::endl;
		Output_seq_num_list( "working_moving_partition_pos = ", job_parameters_->working_moving_partition_pos(), TR.Debug );
		Output_boolean( "centroid_screen = ", centroid_screen_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "allow_base_pair_only_centroid_screen = ", allow_base_pair_only_centroid_screen_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "VDW_atr_rep_screen = ", VDW_atr_rep_screen_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "sample_both_sugar_base_rotamer_ = ", sample_both_sugar_base_rotamer_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "do_not_sample_multiple_virtual_sugar_ = ", do_not_sample_multiple_virtual_sugar_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "sample_ONLY_multiple_virtual_sugar_ = ", sample_ONLY_multiple_virtual_sugar_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "assert_no_virt_ribose_sampling_ = ", assert_no_virt_ribose_sampling_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( "distinguish_pucker_ ", distinguish_pucker_, TR.Debug ); TR.Debug << std::endl;
		TR.Debug << "--------------------------------" << std::endl;


		Pose const pose_save = pose;
		pose = pose_save; //this recopy is useful for triggering graphics.

		// Sets up scorefunctions for a bunch of different screens
		initialize_scorefunctions();

		//Sept 2, 2010
		if ( rebuild_bulge_mode_ ) remove_virtual_rna_residue_variant_type( pose, job_parameters_->working_moving_res() );

		//////////////////////Build previously virtualize sugar/////////////////////
		// A virtualized sugar occurs when a previous move was a 'floating base'
		// step, which only samples euler angles of a base, but virtualized the
		// attached sugar and any residues connecting that nucleotide to the
		//'instantiated' body of the RNA.

		// There are potentially four different virtualized sugar positions.
		// (1) Prev_sugar (most common case) :
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
		//      combine two chunks that were previously built with SWA. It is
		//      possible that the current nucleotide (in the moving chunk) was
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

		Output_title_text( "Build previously virtualize sugar", TR.Debug );
		bool const Is_prev_sugar_virt = Is_previous_sugar_virtual( pose );
		bool const Is_curr_sugar_virt = Is_current_sugar_virtual( pose ); // occurs less frequently -- typically only when connecting two pieces.
		bool const Is_five_prime_CB_sugar_virt = Is_five_prime_chain_break_sugar_virtual( pose ); // may occur when closing loop
		bool const Is_three_prime_CB_sugar_virt = Is_three_prime_chain_break_sugar_virtual( pose ); // may occur when closing loop
		Size const num_nucleotides =  job_parameters_->working_moving_res_list().size();
		Size const gap_size = job_parameters_->gap_size();

		Output_boolean( " Is_previous_sugar_virt = ", Is_prev_sugar_virt, TR.Debug );
		Output_boolean( " Is_current_sugar_virt = ", Is_curr_sugar_virt, TR.Debug );
		Output_boolean( " Is_five_prime_chain_break_sugar_virt = ", Is_five_prime_CB_sugar_virt, TR.Debug );
		Output_boolean( " Is_three_prime_chain_break_sugar_virt = ", Is_three_prime_CB_sugar_virt, TR.Debug );
		TR.Debug << std::endl;

		///Nov 13, 2010

		Size num_virtual_sugar = 0;
		if ( Is_prev_sugar_virt ) num_virtual_sugar++;
		if ( Is_curr_sugar_virt ) num_virtual_sugar++;
		if ( Is_five_prime_CB_sugar_virt ) num_virtual_sugar++;
		if ( Is_three_prime_CB_sugar_virt ) num_virtual_sugar++;

		TR.Debug << "num_virtual_sugar = " << num_virtual_sugar << std::endl;

		// rd2013 -- Parin, please put in comments explaining these checks
		if ( assert_no_virt_ribose_sampling_ /*run checks*/ ){
			if ( floating_base_ && num_nucleotides == 2 ){ //Hacky..ok the only acception right now is in floating_base_ + dinucleotide mode.
				if ( num_virtual_sugar > 1 ) utility_exit_with_message( "assert_no_virt_ribose_sampling_ == true and floating_base, but num_virtual_sugar > 1" );
				if ( num_virtual_sugar == 1 ){
					if ( Is_prev_sugar_virt == false ) {
						utility_exit_with_message( "assert_no_virt_ribose_sampling_ == true and floating_base and num_virtual_sugar == 1 BUT Is_prev_sugar_virt == false! )" );
					}
				}
			} else{
				if ( num_virtual_sugar != 0 ) utility_exit_with_message( "assert_no_virt_ribose_sampling_ == true but num_virtual_sugar != 0" );
			}
		}

		// rd2013 -- Parin, please put in comments explaining these checks. I don't understand what logic guarantees that these assertions are met.
		if ( Is_curr_sugar_virt ){ //Consistency test.
			//Right now, the only possiblility for Is_curr_sugar_virt==true is when combining two silent_files chunk at chain-break
			if ( gap_size != 0 ) utility_exit_with_message( "Is_curr_sugar_virt == true but gap_size != 0 !!" );
			utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();
			if ( working_moving_partition_pos.size() <= 1 ) utility_exit_with_message( "Is_curr_sugar_virt == true but working_moving_partition_pos.size() <= 1" );
		}
		if ( gap_size != 0 && num_virtual_sugar > 1 ) utility_exit_with_message( "gap_size != 0 but num_virtual_sugar > 1" ); //Obsolete!

		if ( floating_base_ ){
			if ( Is_five_prime_CB_sugar_virt ){ //This is rare since floating_base sampling is not often used at chain-closure step!
				TR.Debug << "WARNING: floating_base_ and Is_five_prime_CB_sugar_virt case. Code not implemented yet, early return!" << std::endl;
				return;
			}

			if ( Is_three_prime_CB_sugar_virt ){ //This is rare since floating_base sampling is not often used at chain-closure step!
				TR.Debug << "WARNING: floating_base_ and Is_three_prime_CB_sugar_virt case. Code not implemented yet, early return!" << std::endl;
				return;
			}

			if ( Is_curr_sugar_virt ){
				utility_exit_with_message( "floating_base_ and Is_curr_sugar_virt case. Code not implemented yet!" );
			}
		}

		if ( do_not_sample_multiple_virtual_sugar_ == true && sample_ONLY_multiple_virtual_sugar_ == true ){
			utility_exit_with_message( "do_not_sample_multiple_virtual_sugar_ == true && sample_ONLY_multiple_virtual_sugar_ == true" );
		}

		// this was just for testing.
		if ( do_not_sample_multiple_virtual_sugar_ == true ){
			if ( num_virtual_sugar > 1 ) return;
		}

		if ( sample_ONLY_multiple_virtual_sugar_ == true ){
			if ( gap_size != 0 ) utility_exit_with_message( "sample_ONLY_multiple_virtual_sugar_ == true but gap_size != 0" );
			if ( num_virtual_sugar <= 1 ) return;
		}



		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const moving_res(  job_parameters_->working_moving_res() );
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
		Size const three_prime_chain_break_res = five_prime_chain_break_res + 1;

		//utility::vector1< pose_data_struct2 > prev_sugar_PDL, five_prime_CB_sugar_PDL, three_prime_CB_sugar_PDL;

		// These parameters will be used for sampling virtual sugars and associated bulge.
		FloatingBaseChainClosureJobParameter prev_sugar_FB_JP = FloatingBaseChainClosureJobParameter();
		FloatingBaseChainClosureJobParameter curr_sugar_FB_JP = FloatingBaseChainClosureJobParameter();
		FloatingBaseChainClosureJobParameter five_prime_CB_sugar_FB_JP = FloatingBaseChainClosureJobParameter();
		FloatingBaseChainClosureJobParameter three_prime_CB_sugar_FB_JP = FloatingBaseChainClosureJobParameter();

		// In following, PDL means pose_data_list.
		if ( Is_prev_sugar_virt ){
			TR.Debug << "previous_sugar floating_base_chain_closure" << std::endl;

			Size const prev_moving_res = ( Is_prepend ) ? ( moving_res + num_nucleotides ) : ( moving_res - num_nucleotides );
			Size const prev_ref_res = ( Is_prepend ) ? ( moving_res + ( num_nucleotides + 2 ) ) : ( moving_res - ( num_nucleotides + 2 ) );

			prev_sugar_FB_JP = FloatingBaseChainClosureJobParameter( prev_moving_res, prev_ref_res );
			prev_sugar_FB_JP.set_base_and_pucker_state( pose, job_parameters_ );

			// do the sampling!
			prev_sugar_FB_JP.PDL = previous_floating_base_chain_closure( pose, prev_sugar_FB_JP, "previous" );

			if ( prev_sugar_FB_JP.PDL.size() == 0 ){
				TR.Debug << "Is_prev_sugar_virt == True but prev_sugar_FB_JP.PDL.size() == 0!" << std::endl;
				return;
			}
		}

		if ( Is_curr_sugar_virt ){ //June 12, 2011
			TR.Debug << "current_sugar floating_base_chain_closure" << std::endl;

			Size const curr_ref_res = ( Is_prepend ) ? ( moving_res - 2 ) : ( moving_res + 2 );

			curr_sugar_FB_JP = FloatingBaseChainClosureJobParameter( moving_res, curr_ref_res );
			curr_sugar_FB_JP.set_base_and_pucker_state( pose, job_parameters_ );

			// do the sampling!
			curr_sugar_FB_JP.PDL = previous_floating_base_chain_closure( pose, curr_sugar_FB_JP, "current" );

			if ( curr_sugar_FB_JP.PDL.size() == 0 ){
				TR.Debug << "Is_curr_sugar_virt == True but curr_sugar_FB_JP.PDL.size() == 0!" << std::endl;
				return;
			}
		}

		if ( Is_five_prime_CB_sugar_virt ){
			TR.Debug << "five_prime_CB_sugar floating_base_chain_closure" << std::endl;

			five_prime_CB_sugar_FB_JP = FloatingBaseChainClosureJobParameter( five_prime_chain_break_res, five_prime_chain_break_res - 2 );
			five_prime_CB_sugar_FB_JP.set_base_and_pucker_state( pose, job_parameters_ );

			// do the sampling!
			five_prime_CB_sugar_FB_JP.PDL = previous_floating_base_chain_closure( pose, five_prime_CB_sugar_FB_JP, "five_prime_CB" );

			if ( five_prime_CB_sugar_FB_JP.PDL.size() == 0 ) {
				TR.Debug << "Is_five_prime_CB_sugar_virt == True but five_prime_CB_sugar_FB_JP.PDL.size() == 0!" << std::endl;
				return;
			}
		}

		if ( Is_three_prime_CB_sugar_virt ){
			TR.Debug << "three_prime_CB_sugar floating_base_chain_closure:" << std::endl;

			three_prime_CB_sugar_FB_JP = FloatingBaseChainClosureJobParameter( three_prime_chain_break_res, three_prime_chain_break_res + 2 );
			three_prime_CB_sugar_FB_JP.set_base_and_pucker_state( pose, job_parameters_ );

			// do the sampling!
			three_prime_CB_sugar_FB_JP.PDL = previous_floating_base_chain_closure( pose, three_prime_CB_sugar_FB_JP, "three_prime_CB" );

			if ( three_prime_CB_sugar_FB_JP.PDL.size() == 0 ){
				TR.Debug << "Is_three_prime_CB_sugar_virt == True but three_prime_CB_sugar_FB_JP.PDL.size() == 0!" << std::endl;
				return;
			}
		}

		if (	num_virtual_sugar > 0 ){
			TR.Debug << "Contain_virtual_sugar == true" << std::endl;
		} else{
			TR.Debug << "Contain_virtual_sugar == false" << std::endl;
		}

		/////////////Sort the pose_data_list by score..should be determined according to torsional potential score.	/////////////////////////
		std::sort( prev_sugar_FB_JP.PDL.begin(), prev_sugar_FB_JP.PDL.end(), sort_criteria );
		std::sort( curr_sugar_FB_JP.PDL.begin(), curr_sugar_FB_JP.PDL.end(), sort_criteria );
		std::sort( five_prime_CB_sugar_FB_JP.PDL.begin(),  five_prime_CB_sugar_FB_JP.PDL.end(), sort_criteria );
		std::sort( three_prime_CB_sugar_FB_JP.PDL.begin(), three_prime_CB_sugar_FB_JP.PDL.end(), sort_criteria );
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if ( floating_base_ ){

			floating_base_sampling( pose, prev_sugar_FB_JP );

		} else{

			standard_sampling_WRAPPER( pose, prev_sugar_FB_JP, curr_sugar_FB_JP, five_prime_CB_sugar_FB_JP, three_prime_CB_sugar_FB_JP );
		}

		pose = pose_save;


		if ( rebuild_bulge_mode_ ){ //Ensure that bulge_res is not virtualized in final output.

			for ( Size n = 1; n <= pose_data_list_.size(); n++ ){
				if ( ( *pose_data_list_[n].pose_OP ).residue( job_parameters_->working_moving_res() ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
					utility_exit_with_message( "working_moving_res: " + string_of( job_parameters_->working_moving_res() ) + " of pose " + string_of( n ) + " is a virtual res!"  );
				}
			}
		}

		TR.Debug << "Total time in StepWiseRNA_ResidueSampler::apply " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		Output_title_text( "Exit StepWiseRNA_ResidueSampler::apply", TR.Debug );

	}


	void
	StepWiseRNA_ResidueSampler::floating_base_sampling( pose::Pose & pose, FloatingBaseChainClosureJobParameter const & prev_sugar_FB_JP ){

		Output_title_text( "Enter StepWiseRNA_ResidueSampler::floating_base_sampling", TR.Debug );

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace core::kinematics;

		clock_t const time_start( clock() );

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                const Real RADS_PER_DEG = numeric::NumericTraits < Real > ::pi() / 180.;

		SilentFileData silent_file_data;

		Size const moving_res(  job_parameters_->working_moving_res() ); // Might not corresponds to user input.
		Size const moving_suite(  job_parameters_->working_moving_suite() ); // dofs betweeen this value and value+1 actually move.
		bool const Is_prepend(  job_parameters_->Is_prepend() );
		bool const Is_internal(  job_parameters_->Is_internal() ); // no cutpoints before or after moving_res.
		Size const gap_size( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
		Size const num_nucleotides(  job_parameters_->working_moving_res_list().size() );
		Size const reference_res( job_parameters_->working_reference_res() ); //the last static_residues that this attach to the moving residues
		Size const floating_base_five_prime_chain_break = ( Is_prepend ) ? moving_res : moving_res - 1; //for floating base chain closure when num_nucleotides=1

		if ( combine_long_loop_mode_ ) utility_exit_with_message( "combine_long_loop_mode_ have not been implement for floating base sampling yet!!" );

		// note, in principle can do single nucleotide too (and was used for testing), but typical use case is dinucleotide.
		bool const Is_dinucleotide = ( num_nucleotides == 2 );

		TR.Debug << " NUM_NUCLEOTIDES = " <<  num_nucleotides << std::endl;
		Output_boolean( " IS_DINUCLEOTIDE = ", Is_dinucleotide, TR.Debug ); TR.Debug << std::endl;
		TR.Debug << " GAP SIZE " << gap_size << std::endl;
		TR.Debug << " MOVING RES " << moving_res << std::endl;
		TR.Debug << " MOVING SUITE " << moving_suite << std::endl;
		Output_boolean( " PREPEND ", Is_prepend, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( " INTERNAL ", Is_internal, TR.Debug ); TR.Debug << std::endl;
		TR.Debug << " REFERENCE_RES " << reference_res << std::endl;
		TR.Debug << " FLOATING_BASE_FIVE_PRIME_CHAIN_BREAK " << floating_base_five_prime_chain_break << std::endl;
		Output_boolean( " build_pose_from_scratch_ = ", build_pose_from_scratch_, TR.Debug ); TR.Debug << std::endl;

		if ( Is_dinucleotide == true && Is_internal == true ) utility_exit_with_message( "Is_dinucleotide == true && Is_internal == true )!!" );
		if ( num_nucleotides != 1 && num_nucleotides != 2 ) utility_exit_with_message( "num_nucleotides != 1 and num_nucleotides != 2" );

		if ( Is_dinucleotide == true && allow_base_pair_only_centroid_screen_ == true ){ //Feb 09, 2012: FIXED BUG. Used to be "and" instead of "&&"

			Size const user_input_num_pose_kept = num_pose_kept_;
			num_pose_kept_ = 4*num_pose_kept_;

			//TR.Debug << "Accessible conformational space is larger when sampling a dinucleotide " << std::endl;
			TR.Debug << "allow_base_pair_only_centroid_screen_ == true for floating base + dinucleotide sampling mode " << std::endl;
			TR.Debug << "Increase num_pose_kept by 4 folds" << std::endl;

			TR.Debug << " user_input_num_pose_kept = " << user_input_num_pose_kept << " num_pose_kept_ " << num_pose_kept_ << std::endl;

		}

		/////////////////////////////////////////////Virtual moving res (This will be moved to PoseSetup later)/////////////////////////////
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", moving_res ); //This is unique to floating_base_mode
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res ); //This is unique to floating_base_mode
		Add_virtual_O2Star_hydrogen( pose ); //Apr 3, 2010

	  //////////////////////////////////////////Setup Atr_rep_screening/////////////////////////////////////////////////
		//This is getting annoying...need to move this up here since it create a chainbreak conflict with setup_chain_break_jump_point) below..

		//Real base_rep_score(-9999999999), base_atr_score(-9999999999); //Feb 02, 2012 This might lead to server-test error at R47200
		Real base_rep_score(  - 999999 ), base_atr_score(  - 999999 ); //Feb 02, 2012
		get_base_atr_rep_score( pose, base_atr_score, base_rep_score );

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if ( num_nucleotides == 1 ){ //Sept 16, 2010
			//PURPOSE OF THIS IS FOR FLOATING BASE CHAIN CLOSURE
			//NEED TO ADD THIS BEFORE calculating moving_rsd_at_origin_list since adding OVL1, OVL2 and OVU1 will change atoms number and positions in the residue

			TR.Debug << "setup_chain_break_jump_point for floating_base_five_prime_chain_break = " <<  floating_base_five_prime_chain_break << std::endl;
			//pose.dump_pdb( "pose_before_setup_floating_base_chain_break.pdb" );

			setup_chain_break_jump_point( pose, moving_res, reference_res, floating_base_five_prime_chain_break, true );
			//OK THIS IS HACKY...still have to find a way to remove the cutpoint (perhap after MINIMIZER?)
			//pose.dump_pdb( "pose_after_setup_floating_base_chain_break.pdb" );
		}

		/////////////////////////////////////////////Setup reference stub///////////////////////////////////////////////////////////////////

		core::kinematics::Stub const reference_stub = get_reference_stub( reference_res, pose );

		StepWiseRNA_VDW_BinScreenerOP VDW_bin_screener = new StepWiseRNA_VDW_BinScreener();
		VDW_bin_screener->setup_using_working_pose( pose, job_parameters_ );

		user_input_VDW_bin_screener_->reference_xyz_consistency_check( reference_stub.v );

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		print_atom_info( pose, moving_res, "pose setup_residue_at_origin_list" );


		////////////////////////////////////////Screening poses///////////////////////////////////////////////////


		pose::Pose screening_pose = pose; //Hard copy
		if ( gap_size == 0 ) pose::add_variant_type_to_pose_residue( screening_pose, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res + 1 ); //May 31, 2010

		print_atom_info( screening_pose, moving_res, "screening_pose" );


		pose::Pose ribose_screening_pose = pose; //Hard copy
		if ( gap_size == 0 ) pose::add_variant_type_to_pose_residue( ribose_screening_pose, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res + 1 ); //May 31, 2010
		pose::remove_variant_type_from_pose_residue( ribose_screening_pose, "VIRTUAL_RIBOSE", moving_res );

		print_atom_info( ribose_screening_pose, moving_res, "ribose_screening_pose" );


		if ( gap_size == 0 )	{ //harmonic angle and distnace constraints are used ONLY by chainbreak_screening

			pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res ); //May 31, 2010

			if ( moving_res == ( five_prime_chain_break_res + 1 ) ){ //prepend.
				pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", moving_res ); //this virtual_phosphate was added to pose at the beginning of this function.
			}

			if ( pose.residue_type( five_prime_chain_break_res + 1 ).has_variant_type( "VIRTUAL_PHOSPHATE" ) ){
				utility_exit_with_message( "pose have VIRTUAL_PHOSPHATE AT five_prime_chain_break_res + 1!" );
			}

			TR.Debug << "Adding harmonic chainbreak to standard job_params five_prime_chain_break_res = " << five_prime_chain_break_res << std::endl;
		 	Add_harmonic_chainbreak_constraint( pose, five_prime_chain_break_res );

		}

		if ( num_nucleotides == 1 ){ //Sept 16, 2010

			pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", moving_res ); //May 31, 2010

			if ( moving_res == ( floating_base_five_prime_chain_break + 1 ) ){
				pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", moving_res ); //this virtual_phosphate was added to pose at the beginning of this function.
			}

			if ( pose.residue_type( floating_base_five_prime_chain_break + 1 ).has_variant_type( "VIRTUAL_PHOSPHATE" ) ){
				utility_exit_with_message( "pose have VIRTUAL_PHOSPHATE AT floating_base_five_prime_chain_break + 1!" );
			}

			TR.Debug << "Adding harmonic chainbreak to floating_base_five_prime_chain_break = " <<  floating_base_five_prime_chain_break << std::endl;
		 	Add_harmonic_chainbreak_constraint( pose, floating_base_five_prime_chain_break );

		}

		print_atom_info( pose, moving_res, "pose after modify variant_types" );


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		pose::Pose pose_with_original_HO2star_torsion;
		pose::Pose o2star_pack_pose;
		if ( perform_o2star_pack_ ) {
			pose_with_original_HO2star_torsion = pose;

			if ( use_green_packer_ ) utility_exit_with_message( "green packer mode have not been tested for floating base sampling!" );

			//need to un-virtualize the hydrogen (except for the moving base hydrogen, since set_base_coordinate_frame() requires the atom numbering to be consistent
			//ACTUALLY COULD ALSO UN-VIRTUALIZE THE MOVING_BASE O2STAR_HYDROGEN PROVIDED THAT WE CREATE A SEPERATE moving_rsd_at_origin.
			//BUT NO POINT IN DURING SO, SINCE THE FULL RIBOSE IS VIRTUAL ANYWAYS!

			for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

				//SML PHENIX conference
				if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
					if ( !pose.residue( seq_num ).is_RNA() ) continue;
				}
				if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code.

				if ( seq_num == moving_res ) continue; //moving_res is actually the working_moving_res

				if ( pose.residue_type( seq_num ).has_variant_type( "VIRTUAL_O2STAR_HYDROGEN" ) ){
					pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", seq_num );
				}
			}

			o2star_pack_pose = pose; //NEED THE VARIANT TYPE OF CB_screening_pose and pose to be the same!

		}

		pose::Pose CB_screening_pose = pose; //NEED THE VARIANT TYPE OF CB_screening_pose and pose to be the same!

		print_atom_info( pose, moving_res, "pose after remove o2star variant type" );


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		///////Setup Residue of moving and reference of various rsd conformation (syn/anti chi, 2' and 3' endo) with base at origin coordinate frame////////////////////

		utility::vector1 < core::conformation::ResidueOP > const moving_rsd_at_origin_list
																												= setup_residue_at_origin_list( pose, moving_res, extra_chi_ );
		//TR.Debug << "Change to Residue const & moving_rsd_at_origin=(*moving_rsd_at_origin_list[5]) " << std::endl;
		//Residue const & moving_rsd_at_origin=(*moving_rsd_at_origin_list[5]);

		//March 26, 2011.. after move create screening poses and modify pose variant to ABOVE
		utility::vector1 < core::conformation::ResidueOP > const screening_moving_rsd_at_origin_list
																												= setup_residue_at_origin_list( screening_pose, moving_res, extra_chi_ );

		TR.Debug << "Change to Residue const & screening_moving_rsd_at_origin = ( *screening_moving_rsd_at_origin_list[5] ) " << std::endl;
		Residue const & screening_moving_rsd_at_origin = ( *screening_moving_rsd_at_origin_list[5] );
		//Doesn't really matter which sugar/base conformation the rsd has...any member of the list is fine. (NEED TO CHECK THAT THIS WORKS!)
		//UMM, in score vs. rmsd (Trail 7(old) vs. Trail 15(new))....rmsd and disrimination slight worst....could be an artifact of new sugar position.....
		//Should calculate rmsd...excluding the sugar to see if this is the case..MAKE SURE IT is not something deeper that needs to be fix...Apr 10,2010

		utility::vector1 < core::conformation::ResidueOP > const ribose_screening_moving_rsd_at_origin_list
																												= setup_residue_at_origin_list( ribose_screening_pose, moving_res, extra_chi_ );

		if ( moving_rsd_at_origin_list.size() != screening_moving_rsd_at_origin_list.size() ){
			TR.Debug << "moving_rsd_at_origin_list.size() = " << moving_rsd_at_origin_list.size() << std::endl;
			TR.Debug << "screening_moving_rsd_at_origin_list.size() = " << screening_moving_rsd_at_origin_list.size() << std::endl;
			utility_exit_with_message( "moving_rsd_at_origin_list.size() != screening_moving_rsd_at_origin_list.size()" );
		}

		if ( moving_rsd_at_origin_list.size() != ribose_screening_moving_rsd_at_origin_list.size() ){
			TR.Debug << "moving_rsd_at_origin_list.size() = " << moving_rsd_at_origin_list.size() << std::endl;
			TR.Debug << "ribose_screening_moving_rsd_at_origin_list.size() = " << ribose_screening_moving_rsd_at_origin_list.size() << std::endl;
			utility_exit_with_message( "moving_rsd_at_origin_list.size() != ribose_screening_moving_rsd_at_origin_list.size()" );
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//int const euler_angle_bin_size_MOCK=180;
		//int const euler_angle_bin_min_MOCK=-180/(euler_angle_bin_size_MOCK);
		//int const euler_angle_bin_max_MOCK=180/(euler_angle_bin_size_MOCK)-1;

		// definition of euler angle grid search parameters
		// Following are set by  #define in StepWiseRNA_FloatingBase_Samper_util.hh
		// instead this should be its own class, and these should be private variables.
		Real euler_angle_bin_size_ = euler_angle_bin_size;
		Real euler_z_bin_size_     = euler_z_bin_size;
		Real centroid_bin_size_    = centroid_bin_size;

		if ( fast_ || medium_fast_ ){ // use coarser search for speed
			euler_angle_bin_size_ *= 4;
			euler_z_bin_size_     *= 4;
			centroid_bin_size_    *= 4;
		}

		int const euler_angle_bin_min =  - 180/euler_angle_bin_size_; //Should be -180/euler_angle_bin_size
		int const euler_angle_bin_max = 180/euler_angle_bin_size_ - 1;  //Should be 180/euler_angle_bin_size-1

		int const euler_z_bin_min = int(  - 1/euler_z_bin_size_ );
		int const euler_z_bin_max = int( 1/euler_z_bin_size_ );

		Real C5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list, " C5'" );
		Real O5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list, " O3'" );

		Real const Max_O3_to_C5_DIST = ( num_nucleotides == 1 ) ? O3I_C5I_PLUS_ONE_MAX_DIST : O3I_C5IPLUS2_MAX_DIST;

		Real const max_distance = Max_O3_to_C5_DIST + C5_centroid_dist + O5_centroid_dist + 1; //Theoretical maximum dist between the two base's centroid, +1 is to be lenient

		TR.Debug << "max centroid to centroid distance: " << max_distance << std::endl; ;

		int const centroid_bin_min = int(  - max_distance/centroid_bin_size_ );
		int const centroid_bin_max = int( max_distance/centroid_bin_size_ ) - 1;

		TR.Debug << "euler_angle_bin min = " << euler_angle_bin_min << " max " << euler_angle_bin_max << std::endl;
		TR.Debug << "euler_z_bin_min min = " << euler_z_bin_min << " max " << euler_z_bin_max << std::endl;
		TR.Debug << "centroid_bin min = " << centroid_bin_min << " max " << centroid_bin_max << std::endl;

		std::map< Base_bin, int, compare_base_bin > base_bin_map;
		std::map< Base_bin, int, compare_base_bin > ::const_iterator it;


//////////////////////Probably sure convert to use the class soon...////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////
		// rd2013 -- why is this not using the BaseCentroidScreener class?
		///////////////////////////////////////////////////////////////////////

		// places where floating base can 'dock'
		utility::vector1 < core::kinematics::Stub > other_residues_base_list;


		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ){ //Fang's electron density code
				TR.Debug << "Residue.aa() " << seq_num << " is core::chemical::aa_vrt!" << std::endl;
				continue;
			}
			//SML PHENIX conference
			if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
				if ( !pose.residue( seq_num ).is_RNA() ) continue;
			}

			conformation::Residue const & residue_object = pose.residue( seq_num );
			if ( residue_object.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				TR.Debug << "Residue " << seq_num << " is a VIRTUAL_RNA_RESIDUE!" << std::endl;
				continue;
			}

			if ( seq_num == moving_res ){
				TR.Debug << "Residue " << seq_num << " is a MOVING RESIDUE!" << std::endl;
				continue;
			}

			if ( ( job_parameters_->working_terminal_res() ).has_value( seq_num ) ){
				TR.Debug << "Residue " << seq_num << " is a TERMINAL_RESIDUE!" << std::endl;
				continue;
			}

			core::kinematics::Stub base_info;
			base_info.v = core::chemical::rna::get_rna_base_centroid( residue_object, true );
			base_info.M = core::chemical::rna::get_rna_base_coordinate_system( residue_object, base_info.v );


			//	Real distance_square=( base_info.v - pose.residue(reference_res).xyz(" C5'") ).length_squared();
			//Fix on Jan 15, 2011..how did I let allow this error to exist for so long!
			Real const distance_square = ( base_info.v - core::chemical::rna::get_rna_base_centroid( pose.residue( reference_res ), false ) ).length_squared();


			Real const distance = std::sqrt( distance_square );

			TR.Debug << "distance = " << distance;

			if (  distance > ( max_distance + 6.364 ) ) { //6.364 is max centroid to centroid distance for base-stacking centroid_screening.
				TR.Debug << " Residue " << seq_num << " is too far from the reference res: " <<  reference_res << std::endl;
				continue;
			}

			TR.Debug << " Add to other_residues_base_list: " << seq_num << std::endl;

			other_residues_base_list.push_back( base_info );
		}
		//create screening pose and modify pose variant used to be AFTER THIS POINT.



		core::kinematics::Stub moving_res_base_stub;

		Euler_angles euler_angles;

		Base_bin base_bin;

		Size total_screen_bin( 0 );

		Real /*current_score( 0.0 ), */ delta_rep_score( 0.0 ), delta_atr_score( 0.0 );
		utility::vector1< pose_data_struct2 > pose_data_list;

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Output_title_text( "START FLOATING BASE SAMPLING", TR.Debug );


		for ( base_bin.euler_alpha = euler_angle_bin_min; base_bin.euler_alpha <= euler_angle_bin_max; base_bin.euler_alpha++ ){
		for ( base_bin.euler_z = euler_z_bin_min; base_bin.euler_z <= euler_z_bin_max; base_bin.euler_z++ ){

			if ( fast_ && count_data_.both_count > 1000 ) break;
			if ( medium_fast_ && count_data_.both_count > 1000 && pose_data_list.size() > 5 ) break;
			if ( integration_test_mode_ && count_data_.full_score_count >= 10 ) native_rmsd_screen_ = true;
			if ( integration_test_mode_ && count_data_.rmsd_count >= 10 ) break;

			Matrix O_frame_rotation;
			euler_angles.alpha = ( base_bin.euler_alpha + 0.5 )*euler_angle_bin_size*( RADS_PER_DEG ); //convert to radians
			euler_angles.z = ( base_bin.euler_z )*euler_z_bin_size;
			euler_angles.beta = acos( euler_angles.z );
			euler_angles.gamma = 0;
//			euler_angles.gamma=(base_bin.euler_gamma+0.5)*euler_angle_bin_size*(RADS_PER_DEG); //convert to radians

			convert_euler_to_coordinate_matrix( euler_angles, O_frame_rotation );

			moving_res_base_stub.M = reference_stub.M * O_frame_rotation;

	  ////////////////////////////////////////////////////////////////////////////////////////////////////////
		for ( base_bin.centroid_z = centroid_bin_min; base_bin.centroid_z <= centroid_bin_max; base_bin.centroid_z++ ){
		for ( base_bin.centroid_x = centroid_bin_min; base_bin.centroid_x <= centroid_bin_max; base_bin.centroid_x++ ){
		for ( base_bin.centroid_y = centroid_bin_min; base_bin.centroid_y <= centroid_bin_max; base_bin.centroid_y++ ){

			Vector O_frame_centroid;
			O_frame_centroid[0] = ( base_bin.centroid_x + 0.5 )*centroid_bin_size;
			O_frame_centroid[1] = ( base_bin.centroid_y + 0.5 )*centroid_bin_size;
			O_frame_centroid[2] = ( base_bin.centroid_z )*centroid_bin_size;
			moving_res_base_stub.v = ( reference_stub.M * O_frame_centroid ) + reference_stub.v;

			//count_data_.test_count_one++;
			////////////////////Not dependent on euler gamma value//////////////////////////////////////////////////////////////
			if ( ( moving_res_base_stub.v - reference_stub.v ).length_squared() > max_distance*max_distance ) continue;

			//count_data_.test_count_two++;

	    ///////////////////Current implementation of Base_centroid_screen is not dependent on euler_gamma//////////////////////////////////////
			if ( centroid_screen_ ){
				if ( Base_centroid_screening( moving_res_base_stub, other_residues_base_list, num_nucleotides, count_data_, allow_base_pair_only_centroid_screen_ ) == false ) continue;
			}

//////////////////////Update the moving_res_base_stub/////////////////////////////////////////
			for ( base_bin.euler_gamma = euler_angle_bin_min; base_bin.euler_gamma <= euler_angle_bin_max; base_bin.euler_gamma++ ){


			count_data_.tot_rotamer_count++;

			euler_angles.gamma = ( base_bin.euler_gamma + 0.5 )*euler_angle_bin_size*( RADS_PER_DEG ); //convert to radians


			convert_euler_to_coordinate_matrix( euler_angles, O_frame_rotation );
			moving_res_base_stub.M = reference_stub.M * O_frame_rotation;

//////////////////////////////////////////////////////////////////////////////////////////////////////

			//WHY IS THIS SLOW, CAN IT BE MADE FASTER???
			// rd2013  The two functions (check_floating_base_chain_closable) copy code, should be integrated.
			if ( prev_sugar_FB_JP.sample_sugar ){
				if (	check_floating_base_chain_closable( reference_res, prev_sugar_FB_JP.PDL, moving_rsd_at_origin_list, moving_res_base_stub, Is_prepend, ( num_nucleotides - 1 ) ) == false ) continue;
			} else{
				if (	check_floating_base_chain_closable( reference_res, ribose_screening_pose, moving_rsd_at_origin_list, moving_res_base_stub, Is_prepend, ( num_nucleotides - 1 ) ) == false ) continue;
			}

			//Right now this only work for the case where the moving base is a single residue element...
			if ( gap_size == 0 ){ //WAIT this needs to be coupled to the line above!!!...OK since there is a stricter check below...
				Size const chain_break_reference_res = ( Is_prepend ) ? five_prime_chain_break_res : five_prime_chain_break_res + 1;
				if ( check_floating_base_chain_closable( chain_break_reference_res, ribose_screening_pose, moving_rsd_at_origin_list, moving_res_base_stub, !Is_prepend, 0 /*gap_size*/ ) == false ) continue;
			}

			count_data_.chain_closable_count++;
			//////////////////////////////////////////////////////
			//Feb 21, 2011:
			//When gap_size!=0, there is a always a virtual res serving as buffer between moving_res and surrounding_bin. Virtual res is ignored in create_VDW_screen_bin()
			//Also the phosphate at 3' prime of working_moving_res_list is ignored.
			//Potential error for gap_size==0:
			//VDW_rep doesn't realize that Phosphate and O3' are chain break should be covalently bonded (bond length < Sum VDW)
			//This beg the question whether normal Rosetta FA_REP realize this or does it score the P-O3' clash?
			//This is aside from the fact that the 3' phosphate pos is not yet defined.
			//The code works as long as the 3' phosphate is ignored
			//For append:  3' phosphate is in surrounding_bin res which is next to res in working_moving_res_list and is ignored in create_VDW_screen_bin()
			//For prepend: 3' phosphate is in sampling_res is virtualized and hence is ignored in the function VDW_bin_screener->VDW_rep_screen()
			if ( VDW_bin_screener->VDW_rep_screen( screening_pose, moving_res, screening_moving_rsd_at_origin, moving_res_base_stub ) == false ) continue;

			if ( ( user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ) && ( gap_size != 0 ) && ( Is_internal == false ) ){
				//Does not work for chain_closure move and Is_internal move yet...
				//Residue at 3' of building region have a phosphate that is NOT VIRTUALIZED. This Residue should not be excluded in VDW_bin_screen_pose! Feb 21, 2011.
				//Residue at 5' of building region also have O3' atom that is covalently bond to the phosphate atom of loop res next to it. VDW_rep doesn't realize this and this lead to clash. Hence This Residue should not be excluded in the VDW_bin_screen_pose as well. Feb 21, 2011.

				if ( user_input_VDW_bin_screener_->VDW_rep_screen( screening_pose, moving_res, screening_moving_rsd_at_origin, moving_res_base_stub ) == false ) continue;
			}

			count_data_.good_bin_rep_count++;

			set_base_coordinate_frame( screening_pose, moving_res, screening_moving_rsd_at_origin, moving_res_base_stub );

			//screening_pose.dump_pdb( "FIXED_screening_pose.pdb" );
			//////////////////////////////////////////////////////

			if ( native_rmsd_screen_ && get_native_pose() ){
				//This assumes that screening_pose and native_pose are already superimposed.
				if ( suite_rmsd( *get_native_pose(), screening_pose, moving_res, Is_prepend, true /*ignore_virtual_atom*/ ) > ( native_screen_rmsd_cutoff_ ) ) continue;
				if ( rmsd_over_residue_list( *get_native_pose(), screening_pose, job_parameters_, true /*ignore_virtual_atom*/ ) > ( native_screen_rmsd_cutoff_ ) ) continue; //Oct 14, 2010

				count_data_.rmsd_count++;
				if ( verbose_ ) TR.Debug << "rmsd_count = " << count_data_.rmsd_count << " total count = " << count_data_.tot_rotamer_count << std::endl;
			}

			if ( !Full_atom_van_der_Waals_screening( screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, gap_size, Is_internal ) ) continue;


			//van_der_Waal_screening, for ribose clash and CCD loop closure if gap_size==0////////////////////////////////////////////////////////////////////////////

			for ( Size n = 1; n <= moving_rsd_at_origin_list.size(); n++ ){

				std::string tag = "U_" + lead_zero_string_of( count_data_.tot_rotamer_count, 12 ) + '_' + string_of( n );

				set_base_coordinate_frame( ribose_screening_pose, moving_res, ( *ribose_screening_moving_rsd_at_origin_list[n] ), moving_res_base_stub );

				//ribose_screening_pose.dump_pdb( "FIXED_ribose_screening_pose_"+tag+".pdb" );

				//OK check that with this sugar, the chain can be theoretically closed..
				std::string const moving_atom_name = ( Is_prepend ) ? "O3'" : " C5'";
				std::string const reference_atom_name = ( Is_prepend ) ? " C5'" : "O3'";

				if ( gap_size == 0 ) if ( Check_chain_closable( ribose_screening_pose, five_prime_chain_break_res, 0 ) == false ) continue;

				if ( prev_sugar_FB_JP.sample_sugar ){

					bool pass_for_loop_screen_1 = false;

					for ( Size prev_sugar_ID = 1; prev_sugar_ID <= prev_sugar_FB_JP.PDL.size(); prev_sugar_ID++ ){
						if ( !Check_chain_closable( ribose_screening_pose.residue( moving_res ).xyz( moving_atom_name ),
																	  ( *prev_sugar_FB_JP.PDL[prev_sugar_ID].pose_OP ).residue( reference_res ).xyz( reference_atom_name ), ( num_nucleotides - 1 ) ) ) continue;
						pass_for_loop_screen_1 = true;
						break;
					}

					if ( pass_for_loop_screen_1 == false ) continue;

					//OK if fail van_der_Waal_screening without previous_moving_res sugar...then should fail WITH the previous_moving_res sugar as well..
					if ( !Full_atom_van_der_Waals_screening_REPLICATE( ribose_screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, gap_size, Is_internal ) ) continue;

					bool pass_for_loop_screen_2 = false;

					//Ok, since prev_sugar_FB_JP.PDL is sorted by SCORE, the lower energy conformations are tried first!
					for ( Size prev_sugar_ID = 1; prev_sugar_ID <= prev_sugar_FB_JP.PDL.size(); prev_sugar_ID++ ){

						//Add this statement on Dec 8, 2010//////////////////////
						if ( !Check_chain_closable( ribose_screening_pose.residue( moving_res ).xyz( moving_atom_name ),
											 						  ( *prev_sugar_FB_JP.PDL[prev_sugar_ID].pose_OP ).residue( reference_res ).xyz( reference_atom_name ), ( num_nucleotides - 1 ) ) ) continue;
						/////////////////////////////////////////////////////////

						copy_bulge_res_and_ribose_torsion( prev_sugar_FB_JP, ribose_screening_pose, ( *prev_sugar_FB_JP.PDL[prev_sugar_ID].pose_OP ) );

						// rd2013 why is this REPLICATE?
						if ( !Full_atom_van_der_Waals_screening_REPLICATE( ribose_screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, gap_size, Is_internal ) ) continue;

						pose::add_variant_type_to_pose_residue( ribose_screening_pose, "VIRTUAL_RIBOSE", prev_sugar_FB_JP.moving_res ); // copy_bulge_res_and_ribose_torsion removed the variant type

						copy_bulge_res_and_ribose_torsion( prev_sugar_FB_JP, pose, ( *prev_sugar_FB_JP.PDL[prev_sugar_ID].pose_OP ) );
						pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", prev_sugar_FB_JP.moving_res ); // copy_bulge_res_and_ribose_torsion removed the variant type

						if ( perform_o2star_pack_ ){ //Is this really necessary given that these DOF are virtual anyways!
							copy_bulge_res_and_ribose_torsion( prev_sugar_FB_JP, o2star_pack_pose, ( *prev_sugar_FB_JP.PDL[prev_sugar_ID].pose_OP ) );
							pose::add_variant_type_to_pose_residue( o2star_pack_pose, "VIRTUAL_RIBOSE", prev_sugar_FB_JP.moving_res ); // copy_bulge_res_and_ribose_torsion removed the variant type
						}

						tag += prev_sugar_FB_JP.PDL[prev_sugar_ID].tag;
						pass_for_loop_screen_2 = true;
						break;
					}

					if ( pass_for_loop_screen_2 == false ) continue;


				} else{
					if ( !Check_chain_closable( ribose_screening_pose.residue( moving_res ).xyz( moving_atom_name ), ribose_screening_pose.residue( reference_res ).xyz( reference_atom_name ), ( num_nucleotides - 1 ) ) ) continue;
					// rd2013 why is this REPLICATE?
					if ( !Full_atom_van_der_Waals_screening_REPLICATE( ribose_screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, gap_size, Is_internal ) ) continue;
				}

				count_data_.non_clash_ribose++;   //OK CLEARLY base_rep_score and base_atr_score is not correct!

				// above was all on screening poses; the actual base of the pose hasn't been moved yet!

				set_base_coordinate_frame( CB_screening_pose, moving_res, ( *moving_rsd_at_origin_list[n] ), moving_res_base_stub );

				set_base_coordinate_frame( pose, moving_res, ( *moving_rsd_at_origin_list[n] ), moving_res_base_stub );

				if ( perform_o2star_pack_ ) set_base_coordinate_frame( o2star_pack_pose, moving_res, ( *moving_rsd_at_origin_list[n] ), moving_res_base_stub );

				////////////////////////////////////////////////////////////////////////////////////////////////////////

				if ( gap_size == 0 ){
					if ( !Chain_break_screening( pose, chainbreak_scorefxn_ ) ) continue;

					Copy_CCD_torsions( pose, CB_screening_pose );
					if ( perform_o2star_pack_ ) Copy_CCD_torsions( o2star_pack_pose, CB_screening_pose );
				}

				if ( num_nucleotides == 1 ){
					if ( !Check_chain_closable_floating_base( CB_screening_pose, CB_screening_pose, floating_base_five_prime_chain_break, 0 ) ) continue; //strict version of the check chain closable.
					if ( !Chain_break_screening_general( CB_screening_pose, chainbreak_scorefxn_, floating_base_five_prime_chain_break ) ) continue;

					Copy_CCD_torsions_general( pose, CB_screening_pose, floating_base_five_prime_chain_break, floating_base_five_prime_chain_break + 1 );
					if ( perform_o2star_pack_ ) Copy_CCD_torsions_general( o2star_pack_pose, CB_screening_pose, floating_base_five_prime_chain_break, floating_base_five_prime_chain_break + 1 );
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////

				if ( perform_o2star_pack_ ){
					sample_o2star_hydrogen( o2star_pack_pose, pose_with_original_HO2star_torsion );
					copy_all_o2star_torsions( pose, o2star_pack_pose ); //Copy the o2star torsions from the o2star_pack_pose to the pose!
				}


				Pose_selection_by_full_score( pose_data_list, pose, tag );

				if ( verbose_ ){
					TR.Debug << tag <<  std::endl;
					Output_data( silent_file_data, silent_file_, tag, true, pose, get_native_pose(), job_parameters_ );
				}

				if ( gap_size != 0 && ( num_nucleotides != 1 ) ) break; //Break once found a valid sugar rotamer, at chain_break(gap_size) keep multiple poses..

			}


			// diagnostics?
			it = base_bin_map.find( base_bin );

			if ( it == base_bin_map.end() ){
				base_bin_map[base_bin] = 1;
				total_screen_bin++;
			} else{
				base_bin_map[base_bin] = base_bin_map[base_bin] + 1;
			}

		}
		}
		}
		}
		}
		}


		Output_title_text( "Final sort and clustering: BEFORE floating base_chain_closure", TR.Debug );

		std::sort( pose_data_list.begin(), pose_data_list.end(), sort_criteria );
		cluster_pose_data_list( pose_data_list );
		if ( pose_data_list.size() > num_pose_kept_ ) pose_data_list.erase( pose_data_list.begin() + num_pose_kept_, pose_data_list.end() );
		TR.Debug << "after erasing.. pose_data_list = " << pose_data_list.size() << std::endl;

		TR.Debug << "floating base sampling time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		if ( gap_size == 0 ){
		 	TR.Debug << " angle_n = " << count_data_.good_angle_count << " dist_n = " << count_data_.good_distance_count;
			TR.Debug << " chain_break_screening = " << count_data_.chain_break_screening_count << std::endl;
		}

		TR.Debug << " stack_n = " << count_data_.base_stack_count << " pair_n = " << count_data_.base_pairing_count;
		TR.Debug << " strict_pair_n = " << count_data_.strict_base_pairing_count << " centroid_n = " << count_data_.pass_base_centroid_screen;
		TR.Debug << " bin_rep = " << count_data_.good_bin_rep_count << " atr = " << count_data_.good_atr_rotamer_count << " rep = " << count_data_.good_rep_rotamer_count;
		TR.Debug << " both = " << count_data_.both_count << " total_bin = " << count_data_.tot_rotamer_count << " total_screen_bin = " << total_screen_bin;
		TR.Debug << "  closable = " << count_data_.chain_closable_count << "  non_clash_ribose = " << count_data_.non_clash_ribose << std::endl;
		TR.Debug << " WARNING centroid_n count is severely UNDERSTIMATED...need to be multiply by ( euler_angle_bin_max - euler_angle_bin_min + 1 ): ";
		TR.Debug << ( euler_angle_bin_max - euler_angle_bin_min + 1 ) << std::endl;

		if ( verbose_ ){ //Don't really need this......May 1, 2010...
			std::string const foldername = "test/";
			if ( system( std::string( "rm -r " + foldername ).c_str() ) == -1 ) {
				TR.Error << "Shell process failed!" << std::endl;
			}
			if ( system( std::string( "mkdir -p " + foldername ).c_str() ) == -1 ) {
				TR.Error << "Shell process failed!" << std::endl;
			}
			Analyze_base_bin_map( base_bin_map, foldername );
		}

		pose_data_list_ = pose_data_list;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_ResidueSampler::standard_sampling_WRAPPER( core::pose::Pose & pose,
														FloatingBaseChainClosureJobParameter const & prev_sugar_FB_JP,
														FloatingBaseChainClosureJobParameter const & curr_sugar_FB_JP,
														FloatingBaseChainClosureJobParameter const & five_prime_CB_sugar_FB_JP,
														FloatingBaseChainClosureJobParameter const & three_prime_CB_sugar_FB_JP ){

		using namespace ObjexxFCL;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace core::scoring;


		pose::Pose const pose_copy = pose;

		utility::vector1< pose_data_struct2 > pose_data_list;
		if ( pose_data_list_.size() > 0 ) {
			TR.Debug << "Previous poses exist in sampler: " << pose_data_list_.size() << std::endl;
			pose_data_list = pose_data_list_;
		}


		if (	( prev_sugar_FB_JP.sample_sugar || curr_sugar_FB_JP.sample_sugar || five_prime_CB_sugar_FB_JP.sample_sugar || three_prime_CB_sugar_FB_JP.sample_sugar ) == false ){

			standard_sampling( pose, pose_data_list, "" );

		} else{ //Case where have to sample virtual sugar...

			if ( prev_sugar_FB_JP.PDL.size() == 0 && curr_sugar_FB_JP.PDL.size() == 0 && five_prime_CB_sugar_FB_JP.PDL.size() == 0 && three_prime_CB_sugar_FB_JP.PDL.size() == 0 ){
				utility_exit_with_message( "pose_data_list is empty for all 4 possible virtual sugar!" );
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			utility::vector1< pose_data_struct2 > starting_pose_data_list;

			Size count = 0;

			for ( Size prev_sugar_ID = 1; prev_sugar_ID <= prev_sugar_FB_JP.PDL.size() || prev_sugar_ID == 1; prev_sugar_ID++ ){
				for ( Size curr_sugar_ID = 1; curr_sugar_ID <= curr_sugar_FB_JP.PDL.size() || curr_sugar_ID == 1; curr_sugar_ID++ ){
					for ( Size five_prime_CB_sugar_ID = 1; five_prime_CB_sugar_ID <= five_prime_CB_sugar_FB_JP.PDL.size() || five_prime_CB_sugar_ID == 1; five_prime_CB_sugar_ID++ ){
						for ( Size three_prime_CB_sugar_ID = 1; three_prime_CB_sugar_ID <= three_prime_CB_sugar_FB_JP.PDL.size() || three_prime_CB_sugar_ID == 1; three_prime_CB_sugar_ID++ ){

							count++;

							pose_data_struct2 start_pose_data;

							start_pose_data.pose_OP = new pose::Pose;
							( *start_pose_data.pose_OP ) = pose_copy;
							pose::Pose & start_pose = ( *start_pose_data.pose_OP );

							start_pose_data.score = 0;
							start_pose_data.tag = "";

							if ( prev_sugar_FB_JP.PDL.size() > 0 ) {
								start_pose_data.tag += prev_sugar_FB_JP.PDL[prev_sugar_ID].tag;
								copy_bulge_res_and_ribose_torsion( prev_sugar_FB_JP, start_pose, ( *prev_sugar_FB_JP.PDL[prev_sugar_ID].pose_OP ) );
							} else{
								start_pose_data.tag += "_null";
							}

							if ( curr_sugar_FB_JP.PDL.size() > 0 ) {
								start_pose_data.tag += curr_sugar_FB_JP.PDL[curr_sugar_ID].tag;
								copy_bulge_res_and_ribose_torsion( curr_sugar_FB_JP, start_pose, ( *curr_sugar_FB_JP.PDL[curr_sugar_ID].pose_OP ) );
							} else{
								start_pose_data.tag += "_null";
							}

							if ( five_prime_CB_sugar_FB_JP.PDL.size() > 0 ){
								start_pose_data.tag += five_prime_CB_sugar_FB_JP.PDL[five_prime_CB_sugar_ID].tag;
							 	copy_bulge_res_and_ribose_torsion( five_prime_CB_sugar_FB_JP, start_pose, ( *five_prime_CB_sugar_FB_JP.PDL[five_prime_CB_sugar_ID].pose_OP ) );
							} else{
								start_pose_data.tag += "_null";
							}

							if ( three_prime_CB_sugar_FB_JP.PDL.size() > 0 ){
								start_pose_data.tag += three_prime_CB_sugar_FB_JP.PDL[three_prime_CB_sugar_ID].tag;
							 	copy_bulge_res_and_ribose_torsion( three_prime_CB_sugar_FB_JP, start_pose, ( *three_prime_CB_sugar_FB_JP.PDL[three_prime_CB_sugar_ID].pose_OP ) );
							} else{
								start_pose_data.tag += "_null";
							}

							starting_pose_data_list.push_back( start_pose_data );
						}
					}
				}
			}

			utility::vector1< FloatingBaseChainClosureJobParameter > sampled_sugar_FB_JP_list;
			if ( prev_sugar_FB_JP.PDL.size() > 0 ) sampled_sugar_FB_JP_list.push_back( prev_sugar_FB_JP );
			if ( curr_sugar_FB_JP.PDL.size() > 0 ) sampled_sugar_FB_JP_list.push_back( curr_sugar_FB_JP );
			if ( five_prime_CB_sugar_FB_JP.PDL.size() > 0 ) sampled_sugar_FB_JP_list.push_back( five_prime_CB_sugar_FB_JP );
			if ( three_prime_CB_sugar_FB_JP.PDL.size() > 0 ) sampled_sugar_FB_JP_list.push_back( three_prime_CB_sugar_FB_JP );


			///////Ok, finally have to remove clashes that may arise due to the fact that the floating base sugar sampling and minimization were done individually of each other///
			minimize_all_sampled_floating_bases( pose, sampled_sugar_FB_JP_list, starting_pose_data_list, sampling_scorefxn_, job_parameters_, true /*virtual_ribose_is_from_prior_step*/ );

			TR.Debug << "starting_pose_data_list.size() = " << starting_pose_data_list.size() << std::endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			SilentFileData silent_file_data;
			for ( Size n = 1; n <= starting_pose_data_list.size(); n++ ){
				pose = ( *starting_pose_data_list[n].pose_OP ); //set viewer_pose;

				////Debug///////////////////////
				if ( verbose_ ){ //Umm..rna_sugar_close score adn geom_sol doesn't check for virtual_res?? lead to slight variation...fix this! May 13, 2010

					std::string starting_pose_tag = "starting_pose" + starting_pose_data_list[n].tag;

					utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();

					pose::Pose debug_pose = pose;

					for ( Size ii = 1; ii <= working_moving_partition_pos.size(); ii++ ){
						pose::add_variant_type_to_pose_residue( debug_pose, "VIRTUAL_RNA_RESIDUE", working_moving_partition_pos[ii] );
					}

					if ( job_parameters_->gap_size() == 0 ) pose::add_variant_type_to_pose_residue( debug_pose, "VIRTUAL_PHOSPHATE", job_parameters_->five_prime_chain_break_res() + 1 );


					( *sampling_scorefxn_ )( debug_pose );
					Output_data( silent_file_data, "ribose_sampling.out", starting_pose_tag, false, debug_pose, get_native_pose(), job_parameters_ );
					Output_data( silent_file_data, "SCORE_ribose_sampling.out", starting_pose_tag, true, debug_pose, get_native_pose(), job_parameters_ );
				}
				//////////////////////////////////////////////////

				standard_sampling( pose, pose_data_list, starting_pose_data_list[n].tag );
			}
		}

		pose_data_list_ = pose_data_list;

		pose = pose_copy;


	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_ResidueSampler::standard_sampling( core::pose::Pose & pose, utility::vector1< pose_data_struct2 > & pose_data_list, std::string const sugar_tag ) {

		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::id;

		Output_title_text( "Enter StepWiseRNA_ResidueSampler::standard_sampling", TR.Debug );

		clock_t const time_start( clock() );

		SilentFileData silent_file_data;

		Size const moving_res(  job_parameters_->working_moving_res() ); // Might not corresponds to user input.
		Size const moving_suite(  job_parameters_->working_moving_suite() ); // dofs betweeen this value and value+1 actually move.
		bool const Is_prepend(  job_parameters_->Is_prepend() );
		bool const Is_internal(  job_parameters_->Is_internal() ); // no cutpoints before or after moving_res.
		Size const actually_moving_res( job_parameters_->actually_moving_res() ); //Now same as moving_res
		Size const gap_size( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */
		utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();
		Size const num_nucleotides(  job_parameters_->working_moving_res_list().size() );
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();


		//Somewhat hacky...basically this check if the pose is build from scratch (for build loop outward), Apr 22, 2010
		build_pose_from_scratch_ = ( job_parameters_->working_sequence().length() == ( num_nucleotides + 1 ) ) ? true: false;


		if ( num_nucleotides != 1 && num_nucleotides != 2 ){
			utility_exit_with_message( "num_nucleotides != 1 and num_nucleotides != 2" );
 		}

		bool const Is_dinucleotide = ( num_nucleotides == 2 );


		if ( Is_dinucleotide == true && allow_base_pair_only_centroid_screen_ == true ){ //Feb 09, 2012: FIXED BUG. Used to be "and" instead of "&&"

			Size const user_input_num_pose_kept = num_pose_kept_;
			num_pose_kept_ = 4*num_pose_kept_;

			//TR.Debug << "Accessible conformational space is larger when sampling a dinucleotide " << std::endl;
			TR.Debug << "allow_base_pair_only_centroid_screen_ == true + dinucleotide sampling" << std::endl;
			TR.Debug << "Note that allow_base_pair_only_centroid_screen_ doesn't effect the screening in standard sampling mode." << std::endl;
			TR.Debug << "Just keeping more pose to be consistent with floating base mode + keep high score basepairing conformations." << std::endl;
			TR.Debug << "Increase num_pose_kept by 4 folds" << std::endl;

			TR.Debug << " user_input_num_pose_kept = " << user_input_num_pose_kept << " num_pose_kept_ " << num_pose_kept_ << std::endl;

		}

		if ( build_pose_from_scratch_ ){
			TR.Debug << "Since build_pose_from_scratch, choose to increase NUM_POSE_KEPT by 36 fold. ";
			TR.Debug << " Somewhat hacky..since sample both sugar..want to make sure that we keep are good energy score states " << std::endl;
			TR.Debug << "Old_num_pose_kept_ = " << num_pose_kept_  << std::endl;
			num_pose_kept_ = 36* num_pose_kept_;
			TR.Debug << "New_num_pose_kept_ = " << num_pose_kept_  << std::endl;
		}

		if ( sample_both_sugar_base_rotamer_ ){
			TR.Debug << "Since build_pose_from_scratch, choose to increase NUM_POSE_KEPT by 12 fold. ";
			TR.Debug << " Somewhat hacky..since sample both sugar..want to make sure that we keep are good energy score states " << std::endl;
			TR.Debug << "Old_num_pose_kept_ = " << num_pose_kept_  << std::endl;
			num_pose_kept_ = 12* num_pose_kept_;
			TR.Debug << "New_num_pose_kept_ = " << num_pose_kept_  << std::endl;
		}


		TR.Debug << " NUM_NUCLEOTIDES = " <<  num_nucleotides << std::endl;
		Output_boolean( " IS_DINUCLEOTIDE = ", Is_dinucleotide, TR.Debug ); TR.Debug << std::endl;
		TR.Debug << " GAP SIZE " << gap_size << std::endl;
		TR.Debug << " MOVING RES " << moving_res << std::endl;
		TR.Debug << " MOVING SUITE " << moving_suite << std::endl;
		Output_boolean( " PREPEND ", Is_prepend, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( " INTERNAL ", Is_internal, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( " allow_bulge_at_chainbreak_ ", allow_bulge_at_chainbreak_, TR.Debug ); TR.Debug << std::endl;
		Output_boolean( " build_pose_from_scratch_ = ", build_pose_from_scratch_, TR.Debug ); TR.Debug << std::endl;

		//for combine_long_loop_mode_
		Size const last_append_res = ( Is_prepend ) ? moving_res - 1: moving_res;
		Size const last_prepend_res = ( Is_prepend ) ? moving_res: moving_res + 1;
		Real const atom_atom_overlap_dist_cutoff =  - 1.0; //value taken from CombineLongLoopFilterer...two atoms in contact if there VDW edge are within 1 angstrom of each other.

		if ( combine_long_loop_mode_ && gap_size != 0 ){ //residue-residue contact screen;
			TR.Debug << "combine_long_loop_mode_ && gap_size == 0" << std::endl;
			TR.Debug << "Enforcing contact between LAST_APPEND_RES: " << last_append_res << " and LAST_PREPEND_RES: " << last_prepend_res  << std::endl;
			TR.Debug << "atom_atom_overlap_dist_cutoff " << atom_atom_overlap_dist_cutoff << std::endl;
		}




		/////////////////////////////// O2star sampling/virtualization //////////////////////////
		Pose pose_with_virtual_O2star_hydrogen = pose;
		Add_virtual_O2Star_hydrogen( pose_with_virtual_O2star_hydrogen );


		//	pose.set_torsion( TorsionID( moving_res, id::CHI, 4 ), 0 );  //This torsion is not sampled. Arbitary set to zero to prevent randomness
		//	Mod out on Apr 3, 2010 Important so that pose can remember the o2star torsion from previous minimization step.
		//This is true only in the INTERNAL CASE!...May 31, 2010

		//if perform_o2star_pack, pose 2'-OH torsion will be sampled!
		pose::Pose pose_with_original_HO2star_torsion;
		pose::Pose o2star_pack_pose;
		if ( perform_o2star_pack_ ) {
			pose_with_original_HO2star_torsion = pose;
			o2star_pack_pose = pose;
			if ( use_green_packer_ ){
				initialize_o2star_green_packer( o2star_pack_pose );
			} else {
				initialize_o2star_packer_task( o2star_pack_pose );
			}
		} else {
			// Otherwise, virtualize the 2-OH.
			pose = pose_with_virtual_O2star_hydrogen;
		}

		//Real base_rep_score(-9999999999), base_atr_score(-9999999999); //Feb 02, 2012 This might lead to server-test error at R47200
		Real base_rep_score(  - 999999 ), base_atr_score(  - 999999 ); //Feb 02, 2012


		get_base_atr_rep_score( pose_with_virtual_O2star_hydrogen, base_atr_score, base_rep_score );
	  //////////////////////////////////////////Setup Atr_rep_screening/////////////////////////////////////////////////

		pose::Pose screening_pose = pose_with_virtual_O2star_hydrogen; //Hard copy

		//Necessary for the case where gap_size == 0. In this case, Pose_setup does not automatically create a VIRTUAL_PHOSPHATE.///
		//However since screening_pose is scored before the CCD corrrectly position the chain_break phosphate atoms, ///////////////
		// the VIRTUAL_PHOSPHATE is needed to prevent artificial crashes. Parin Jan 28, 2010////////////////////////////////////
		if ( gap_size == 0 ) pose::add_variant_type_to_pose_residue( screening_pose, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res + 1 );

		///////////////////////////////////////Setup chainbreak_screening//////////////////////////////////////////////////////
		pose::Pose chain_break_screening_pose = pose_with_virtual_O2star_hydrogen; //Hard copy

		if ( gap_size == 0 )	{ //harmonic angle and distnace constraints are used ONLY by chainbreak_screening
			TR.Debug << "five_prime_chain_break_res = " << five_prime_chain_break_res << std::endl;
		 	Add_harmonic_chainbreak_constraint( chain_break_screening_pose, five_prime_chain_break_res );
		}


		/////Get the Rotamer Sampler/////
		rotamer_sampler::rna::RNA_SuiteRotamerOP
				sampler_op( setup_rotamer_sampler( pose ) );
		rotamer_sampler::rna::RNA_SuiteRotamer & sampler( *sampler_op );

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MAIN LOOP --> rotamer sampling.
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////

		Real /*current_score( 0.0 ), */ delta_rep_score( 0.0 ), delta_atr_score( 0.0 );
		Size ntries( 0 ), max_ntries( 10000 ), num_success( 0 ); // used in choose_random mode.

		for ( sampler.reset(); sampler.not_end(); ++sampler ) {
			++ntries;
			if ( choose_random_ && ntries > max_ntries ) break;

			sampler.apply( screening_pose );
			count_data_.tot_rotamer_count++;

			if ( fast_ && count_data_.both_count >= 100 ) break;
			if ( medium_fast_ && count_data_.both_count >= 1000 ) break;
			if ( integration_test_mode_ && count_data_.full_score_count >= 10 ) native_rmsd_screen_ = true;
			if ( integration_test_mode_ && count_data_.rmsd_count >= 10 ) break;

			std::string tag = create_tag( "U" + sugar_tag, ntries );


			if ( native_rmsd_screen_ && get_native_pose() ){
				//This assumes that screening_pose and native_pose are already superimposed.
				if ( suite_rmsd( *get_native_pose(), screening_pose, actually_moving_res, Is_prepend ) > ( native_screen_rmsd_cutoff_ ) ) continue;
				if ( rmsd_over_residue_list( *get_native_pose(), screening_pose, job_parameters_, false ) > ( native_screen_rmsd_cutoff_ ) ) continue; //Oct 14, 2010

				count_data_.rmsd_count++;
				if ( verbose_ ) TR.Debug << "rmsd_count = " << count_data_.rmsd_count << " total count = " << count_data_.tot_rotamer_count << std::endl;
			}

			if ( combine_long_loop_mode_ && gap_size != 0 ){ //residue-residue contact screen;
//Nov 18,2010
				if ( Is_residues_in_contact( last_append_res, screening_pose, last_prepend_res, screening_pose, atom_atom_overlap_dist_cutoff, 1 /*num_atom_contacts_cutoff*/ ) == false ){
					continue;
				}
				count_data_.residues_contact_screen++; //mistakenly put this inside the if loop, fix on Sept 22, 2010
			}

			bool is_possible_bulge = false;

			if ( centroid_screen_ ){

				//Reminder note of independency: Important that base_stub_list is updated even in the case where gap_size == 0 (and bulge is allowed) ////
				// since base_stub_list is used later below in the chain_break_screening section Jan 28, 2010 Parin S. ///////////////////////////////////

				// updates base_stub_list.
				bool found_a_centroid_interaction_partner( false );
				found_a_centroid_interaction_partner = base_centroid_screener_->Update_base_stub_list_and_Check_centroid_interaction( screening_pose, count_data_ );

				if ( gap_size == 0 ){ //special case allow for bulges
					if ( !found_a_centroid_interaction_partner ){ //does not stack or base_pair
						if ( working_moving_partition_pos.size() == 1 ) is_possible_bulge = true;
					}
				}

				if ( num_nucleotides > 1 && is_possible_bulge == true ) utility_exit_with_message( "num_nucleotides > 1 but is_possible_bulge == true!" );

				if ( ( gap_size > 0 || force_centroid_interaction_ ) && !found_a_centroid_interaction_partner ) continue;
				//Essential this doesn't screen for centroid interaction at chain_break.
				//The chain break can be at both a single strand and a double strand. The statement below is stricter and doesn't screen for centroid interaction only if
				//the chainbreak is single stranded.
				//	if(!found_a_centroid_interaction_partner && !is_possible_bulge) continue; //This is the new version

				// Note that is does not update base_stub_list. To do that, use Update_base_stub_list_and_Check_that_terminal_res_are_unstacked
				if ( !base_centroid_screener_->Check_that_terminal_res_are_unstacked() ) continue;

			}

			//////////////////////////////////////////////////////////////////////////////////////////
			///////////////Chain_break_screening -- distance cut                     /////////////////
			//////////////////////////////////////////////////////////////////////////////////////////
			if ( gap_size <= 1 ){

				if ( gap_size == 0 && finer_sampling_at_chain_closure_ == true ){ //hacky, use strict version of check_chain_closable when using finer_sampling.
					if ( !Check_chain_closable_floating_base( screening_pose, screening_pose, five_prime_chain_break_res, gap_size ) ) continue;
				} else{
					if ( !Check_chain_closable( screening_pose, five_prime_chain_break_res, gap_size ) ) continue;
				}
				count_data_.chain_closable_count++;
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			/////////////// Van_der_Waals_screening                        /////////////////
			//////////////////////////////////////////////////////////////////////////////////////////
			if ( !Full_atom_van_der_Waals_screening( screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, gap_size, Is_internal ) ) continue;


			if ( ( user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ) && ( gap_size != 0 ) && ( Is_internal == false ) ){
				//Does not work for chain_closure move and Is_internal move yet...
				//Residue at 3' of building region have a phosphate that is NOT VIRTUALIZED. This Residue should not be excluded in VDW_bin_screen_pose! Feb 21, 2011.
				//Residue at 5' of building region also have O3' atom that is covalently bond to the phosphate atom of loop res next to it. VDW_rep doesn't realize this and this lead to clash. Hence This Residue should not be excluded in the VDW_bin_screen_pose as well. Feb 21, 2011.

				if ( user_input_VDW_bin_screener_->VDW_rep_screen( screening_pose, moving_res ) == false ){
					//TR.Debug << tag << " pass Full_atom_VDW_screening but fail user_input_VDW_bin_screening! " << std::endl;
					continue;
				}
				count_data_.good_bin_rep_count++;
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			// Almost ready to actually score pose.
			//////////////////////////////////////////////////////////////////////////////////////////
			sampler.apply( pose );
			if ( perform_o2star_pack_ ) sampler.apply( o2star_pack_pose );

			//////////////////////////////////////////////////////////////////////////////////////////
			///////////////Chain_break_screening -- CCD closure /////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////
			bool bulge_added( false );
			if ( gap_size == 0 /*really need to close it!*/ ){

				sampler.apply( chain_break_screening_pose );
				if ( ! Chain_break_screening( chain_break_screening_pose, chainbreak_scorefxn_ ) ) continue;

				// Need to be very careful here -- do CCD torsions ever overlap with pose torsions?
				Copy_CCD_torsions( pose, chain_break_screening_pose );
				if ( perform_o2star_pack_ ) Copy_CCD_torsions( o2star_pack_pose, chain_break_screening_pose );

				if ( is_possible_bulge ){
					bulge_added = apply_bulge_variant( pose, delta_atr_score ); /*further cut on atr, inside*/
					if ( perform_o2star_pack_ ) apply_bulge_variant( o2star_pack_pose, delta_atr_score ); /*further cut on atr, inside*/
				}
			}
			////////////////////////////////////////////////////////////////////

			///////Add pose to pose_data_list if pose have good score///////////
			if ( perform_o2star_pack_ ){
				sample_o2star_hydrogen( o2star_pack_pose, pose_with_original_HO2star_torsion );
				copy_all_o2star_torsions( pose, o2star_pack_pose ); //Copy the o2star torsions from the o2star_pack_pose to the pose!
			}

			if ( include_torsion_value_in_tag_ ) tag += create_rotamer_string( pose );

			Pose_selection_by_full_score( pose_data_list, pose, tag );

			if ( verbose_ ){
				TR.Debug << tag <<  std::endl;
				Output_data( silent_file_data, silent_file_, tag, true, pose, get_native_pose(), job_parameters_ );
			}

			if ( bulge_added ){
				remove_virtual_rna_residue_variant_type( pose, job_parameters_->working_moving_res() );
				if ( perform_o2star_pack_ ) remove_virtual_rna_residue_variant_type( o2star_pack_pose, job_parameters_->working_moving_res() );
			}

			num_success++;
			if ( choose_random_ && num_success >= num_random_samples_ ) break;

	 	} //while( rotamer_generator->has_another_rotamer() )

		//		if ( choose_random_ ) TR << "Number of tries: " << ntries << std::endl;

		Output_title_text( "Final sort and clustering", TR.Debug );
		std::sort( pose_data_list.begin(), pose_data_list.end(), sort_criteria );
		cluster_pose_data_list( pose_data_list );
		if ( pose_data_list.size() > num_pose_kept_ ) pose_data_list.erase( pose_data_list.begin() + num_pose_kept_, pose_data_list.end() );
		TR.Debug << "after erasing.. pose_data_list = " << pose_data_list.size() << std::endl;

		//Hacky..temporary until we fix the reroot atom problem..This is just for calculating rmsd purposes... Apr 27, 2010 Parin/////////////////////////////////////////////
		if ( build_pose_from_scratch_ && get_native_pose() ){
			utility::vector1< core::Size > const & working_best_alignment( job_parameters_->working_best_alignment() );
			pose::Pose const & native_pose = *get_native_pose();

			for ( Size n = 1; n <= pose_data_list.size(); n++ ){ //align all other pose to first pose
				pose::Pose & current_pose = ( *pose_data_list[n].pose_OP );
				std::string const & tag = pose_data_list[n].tag;

				align_poses( current_pose, tag, native_pose, "native", working_best_alignment );
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		TR.Debug << "FINAL COUNTS" << std::endl;
		output_count_data();

		TR.Debug << "Total time in StepWiseRNA_ResidueSampler: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::output_count_data(){

		Size const gap_size( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */

		if ( gap_size <= 1 ) TR.Debug << " chain_closable_count = " << count_data_.chain_closable_count << std::endl;
		if ( gap_size == 0 ){
		 	TR.Debug << " angle_n = " << count_data_.good_angle_count << " dist_n = " << count_data_.good_distance_count;
			TR.Debug << " chain_break_screening = " << count_data_.chain_break_screening_count << std::endl;
		}
		if ( combine_long_loop_mode_ && gap_size != 0 ) TR.Debug << "res_contact = " << count_data_.residues_contact_screen << " ";

		TR.Debug << "stack = " << count_data_.base_stack_count << " pair = " << count_data_.base_pairing_count;
		TR.Debug << " strict_pair_n = " << count_data_.strict_base_pairing_count;
		TR.Debug << " atr = " << count_data_.good_atr_rotamer_count;
		TR.Debug << " rep = " << count_data_.good_rep_rotamer_count;
		TR.Debug << " both = " << count_data_.both_count;
		TR.Debug << " bulge = " << count_data_.bulge_at_chain_closure_count;
		TR.Debug << " rmsd = " << count_data_.rmsd_count << " tot = " << count_data_.tot_rotamer_count << std::endl;

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< pose_data_struct2 >
	StepWiseRNA_ResidueSampler::previous_floating_base_chain_closure( pose::Pose & viewer_pose, FloatingBaseChainClosureJobParameter const & FB_job_params, std::string const name ){

		return sample_virtual_ribose_and_bulge_and_close_chain( viewer_pose, FB_job_params, name,
																							 					scorefxn_, sampling_scorefxn_, atr_rep_screening_scorefxn_, chainbreak_scorefxn_,
																												job_parameters_, true /*virtual_ribose_is_from_prior_step*/ );


	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::kinematics::Stub
	StepWiseRNA_ResidueSampler::get_reference_stub( core::Size const reference_res, pose::Pose const & pose ) const{

		std::string const reference_stub_type = "base"; //"ribose"
		core::kinematics::Stub reference_stub;

		TR.Debug << "-----------------------get reference stub-----------------------" << std::endl;
		if ( reference_stub_type == "ribose" ){
			reference_stub = Get_ribose_stub( pose.residue( reference_res ), job_parameters_->Is_prepend(), true );
		} else{ //Use the base
			reference_stub.v = core::chemical::rna::get_rna_base_centroid(  pose.residue( reference_res ), true );
			reference_stub.M = core::chemical::rna::get_rna_base_coordinate_system( pose.residue( reference_res ), reference_stub.v );
		}

		TR.Debug << " reference_stub.v: x = " << reference_stub.v[0] << " y = " << reference_stub.v[1] << " z = " << reference_stub.v[2] << std::endl;
		TR.Debug << "---------------------------------------------------------------------" << std::endl;

		return reference_stub;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_ResidueSampler::Is_previous_sugar_virtual( core::pose::Pose const & pose ) const {

		//Check if previous sugar is virtual, if virtual then need to sample it.
		bool const Is_prepend(  job_parameters_->Is_prepend() ); // If true, moving_suite+1 is fixed. Otherwise, moving_suite is fixed.
		Size const moving_res(  job_parameters_->working_moving_res() ); // Corresponds to user input.
		Size const num_nucleotides(  job_parameters_->working_moving_res_list().size() );
		Size const previous_moving_res = ( Is_prepend ) ? ( moving_res + num_nucleotides ) : ( moving_res - num_nucleotides );
		Size const previous_bulge_res = ( Is_prepend ) ? ( moving_res + ( num_nucleotides + 1 ) ) : ( moving_res - ( num_nucleotides + 1 ) );

		return Is_ribose_virtual(  pose, previous_moving_res, previous_bulge_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////New June 12, 2011/////
	bool
	StepWiseRNA_ResidueSampler::Is_current_sugar_virtual( core::pose::Pose const & pose ) const {

		// Check if curr sugar is virtual. If virtual then need to sample it.
		// This occur when combining two chunk and the moving_res
		// in the moving_chunk was built with a dinucleotide move.
		bool const Is_prepend( job_parameters_->Is_prepend() );
		Size const moving_res( job_parameters_->working_moving_res() );
		Size const virtual_ribose_res = moving_res;
		Size const bulge_res = ( Is_prepend ) ? ( moving_res - 1 ) : ( moving_res + 1 );

		return Is_ribose_virtual(  pose, virtual_ribose_res, bulge_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_ResidueSampler::Is_five_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) const {

		Size const moving_res( job_parameters_->working_moving_res() );
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
		Size const gap_size( job_parameters_->gap_size() );

		if ( gap_size != 0 ) return false;

		// Make sure to not over count number of virtual_ribose-virtual_bulge
		// pairs to be build!
		if ( moving_res == five_prime_chain_break_res ) return false;

		Size const five_prime_CB_bulge_res = ( five_prime_chain_break_res - 1 );

		bool sugar_is_virtual = Is_ribose_virtual( pose, five_prime_chain_break_res, five_prime_CB_bulge_res );

		return sugar_is_virtual;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_ResidueSampler::Is_three_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) const {

		Size const moving_res(  job_parameters_->working_moving_res() );
		Size const three_prime_chain_break_res = job_parameters_->five_prime_chain_break_res() + 1;
		Size const gap_size( job_parameters_->gap_size() );

		if ( gap_size != 0 ) return false;

		// Make sure to not over count number of virtual_ribose-virtual_bulge
		// pairs to be build!
		if ( moving_res == three_prime_chain_break_res ) return false;

		Size const three_prime_CB_bulge_res = ( three_prime_chain_break_res + 1 );

		bool sugar_is_virtual = Is_ribose_virtual(  pose, three_prime_chain_break_res, three_prime_CB_bulge_res );

		return sugar_is_virtual;

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	void
	StepWiseRNA_ResidueSampler::get_base_atr_rep_score( core::pose::Pose const & pose, core::Real & base_atr_score, core::Real & base_rep_score ){

		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace ObjexxFCL;

		Size const working_moving_suite(  job_parameters_->working_moving_suite() );
		Size const working_moving_res(  job_parameters_->working_moving_res() );

		Size const nres = job_parameters_->working_sequence().size();

		bool const Is_prepend(  job_parameters_->Is_prepend() );

		///////////////////////////////Old_way////////////////////////////////////////////

		pose::Pose base_pose_screen = pose; //hard copy

		if ( output_pdb_ ) base_pose_screen.dump_pdb( "base_atr_rep_before.pdb" );

//		if(working_moving_suite>=nres) utility_exit_with_message( "working_moving_suite " + string_of(working_moving_suite) + " >= nres " + string_of(nres) );

		pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_PHOSPHATE", working_moving_res ); //May 7...

		if ( ( working_moving_res + 1 ) <= nres ){
			pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_PHOSPHATE", working_moving_res + 1 ); //May 7...
		}

		if ( sample_both_sugar_base_rotamer_ == true ){ //Nov 15, 2010
			Size const extra_sample_sugar_base_res = ( Is_prepend ) ? ( working_moving_res + 1 ) : ( working_moving_res - 1 );
			if ( verbose_ ) TR.Debug << "extra_sample_sugar_base_res = " << extra_sample_sugar_base_res << std::endl;
			pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_RIBOSE", extra_sample_sugar_base_res );
		}

		// I think this should work... push apart different parts of the structure so that whatever fa_atr, fa_rep is left is due to "intra-domain" interactions.
		// Crap this doesn't work when building 2 or more nucleotides.

		Size const jump_at_moving_suite = make_cut_at_moving_suite( base_pose_screen, working_moving_suite );
		kinematics::Jump j = base_pose_screen.jump( jump_at_moving_suite );
		j.set_translation( Vector( 1.0e4, 0.0, 0.0 ) );
		base_pose_screen.set_jump( jump_at_moving_suite, j );

		( *atr_rep_screening_scorefxn_ )( base_pose_screen );

		EnergyMap const & energy_map = base_pose_screen.energies().total_energies();
		base_atr_score = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[ scoring::fa_atr ]; //
		base_rep_score = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[ scoring::fa_rep ];
		TR.Debug << "base_rep = " << base_rep_score << " base_atr = " << base_atr_score << std::endl;

		if ( output_pdb_ ) base_pose_screen.dump_pdb( "base_atr_rep_after.pdb" );

		////////////////////////////////////////////////////////////////////////////////////

		/*

		Size const nres = job_parameters_->working_sequence().size();
		utility::vector1< Size > const working_moving_res_list = job_parameters_->working_moving_res_list();


		utility::vector1< Size > partition_0_seq_num_list;
		utility::vector1< Size > partition_1_seq_num_list;


		bool const root_partition = partition_definition( rerooted_fold_tree.root() );

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
//			if(working_moving_res_list.has_value(seq_num)) continue; //Exclude working_moving_residues (the one being sampled..)

			if ( partition_definition( seq_num ) == 0 ){
			 	partition_0_seq_num_list.push_back( seq_num );
			} else if ( partition_definition( seq_num ) == 1 ){
			 	partition_1_seq_num_list.push_back( seq_num );
			} else{
				utility_exit_with_message( "seq_num " + string_of( seq_num ) + " is not both in either partition!!" );
			}
		}

		sort_seq_num_list( partition_0_seq_num_list ); //Low seq_num on the top of the list [1,2,3,4,5]
		sort_seq_num_list( partition_1_seq_num_list ); //Low seq_num on the top of the list [1,2,3,4,5]

		*/

	}


	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::initialize_scorefunctions(){

		initialize_common_scorefxns( scorefxn_, sampling_scorefxn_, atr_rep_screening_scorefxn_, chainbreak_scorefxn_, o2star_pack_scorefxn_ );

	}

	////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< pose_data_struct2 > &
	StepWiseRNA_ResidueSampler::get_pose_data_list(){
		return pose_data_list_;
	}



	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::Copy_CCD_torsions( pose::Pose & pose, pose::Pose const & template_pose ) const {

 		using namespace core::chemical;
		using namespace core::conformation;
	  using namespace core::id;

		Size const five_prime_res = job_parameters_->five_prime_chain_break_res();
		Size const three_prime_res = five_prime_res + 1;


		//Even through there is the chain_break, alpha of 3' and epl and gamma of 5' should be defined due to the existence of the upper and lower variant type atoms.
		Copy_CCD_torsions_general( pose, template_pose, five_prime_res, three_prime_res );

	}


	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::Copy_CCD_torsions_general( pose::Pose & pose, pose::Pose const & template_pose, Size const five_prime_res, Size const three_prime_res ) const {

 		using namespace core::chemical;
		using namespace core::conformation;
	  using namespace core::id;

		if ( ( five_prime_res ) != ( three_prime_res - 1 ) ) utility_exit_with_message( "( five_prime_res ) != ( three_prime_res - 1 )" );

		conformation::Residue const & lower_res = template_pose.residue( five_prime_res );
		conformation::Residue const & upper_res = template_pose.residue( three_prime_res );

		for ( Size n = 1; n <= 3; n++ ){ //alpha, beta, gamma of 3' res
			pose.set_torsion( TorsionID( three_prime_res, id::BB,  n ), upper_res.mainchain_torsion( n ) );
		}

		for ( Size n = 5; n <= 6; n++ ){ //epsilon and zeta of 5' res
			pose.set_torsion( TorsionID( five_prime_res, id::BB,  n ), lower_res.mainchain_torsion( n ) );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_ResidueSampler::Chain_break_screening_general( pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn, Size const five_prime_res ){

		using namespace core::scoring;

 		static protocols::rna::RNA_LoopCloser rna_loop_closer;

		if ( chain_break_screening_pose.residue( five_prime_res ).has_variant_type( chemical::CUTPOINT_LOWER ) == false ) {
			utility_exit_with_message( "chain_break_screening_pose.residue( five_prime_chain_break_res ).has_variant_type(  chemical::CUTPOINT_LOWER ) == false" );
		}

		if ( chain_break_screening_pose.residue( five_prime_res + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) == false ) {
			utility_exit_with_message( "chain_break_screening_pose.residue( five_prime_chain_break_res + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) == false" );
		}

		if ( reinitialize_CCD_torsions_ ) set_CCD_torsions_to_zero( chain_break_screening_pose, five_prime_res );

		//		Real const mean_dist_err=rna_loop_closer.apply( chain_break_screening_pose, five_prime_res);
		rna_loop_closer.apply( chain_break_screening_pose, five_prime_res );

		( *chainbreak_scorefxn )( chain_break_screening_pose );

		scoring::EMapVector & energy_map = chain_break_screening_pose.energies().total_energies();
		Real const angle_score = energy_map[scoring::angle_constraint];
		Real const distance_score = energy_map[scoring::atom_pair_constraint];

		if ( angle_score < 5 ) count_data_.good_angle_count++;
		if ( distance_score < 5 ) count_data_.good_distance_count++;
		if ( ( angle_score < 5 ) && ( distance_score < 5 ) ){
			count_data_.chain_break_screening_count++;
			if ( verbose_ ){
				//				TR.Debug << " C5_O3= " << C5_O3_distance << " C5_O3_n= " << count_data_.C5_O3_distance_count;
				TR.Debug << "  chain_closable_count = " << count_data_.chain_closable_count;
				TR.Debug << " angle = " << angle_score << " dist = " << distance_score;
				TR.Debug << " angle_n = " << count_data_.good_angle_count;
				TR.Debug << " dist_n = " << count_data_.good_distance_count;
				TR.Debug << " chain_break_screening = " << count_data_.chain_break_screening_count;
				TR.Debug << " tot = " << count_data_.tot_rotamer_count << std::endl;
			}
			return true;
		} else {
			return false;
		}
	}


	bool
	StepWiseRNA_ResidueSampler::Chain_break_screening( pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn ){

		Size const five_prime_res = job_parameters_->five_prime_chain_break_res();

		return ( Chain_break_screening_general( chain_break_screening_pose, chainbreak_scorefxn, five_prime_res ) );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_ResidueSampler::apply_bulge_variant( core::pose::Pose & pose, Real const & delta_atr_score ){

		using namespace ObjexxFCL;

		if ( rebuild_bulge_mode_ ) return false; //Hacky want to output sample diverse bulge conformation

//		static Real const atr_cutoff_for_bulge( -1.0 );
//		static Real const atr_cutoff_for_bulge( -1.0 );
//		static Real const atr_cutoff_for_bulge( 0.0 );

		//static Real const atr_cutoff_for_bulge( -9999999999999999999.0 ); //Feb 02, 2012 This might lead to server-test error at R47200
		static Real const atr_cutoff_for_bulge( -999999.0 ); //Feb 02, 2012



//		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
		Size const working_moving_res(  job_parameters_->working_moving_res() );

		if ( delta_atr_score > (  + 0.01 ) ){
			utility_exit_with_message( "delta_atr_score > (  + 0.01 ). delta_atr_score = " + string_of( delta_atr_score ) );
		}

		if ( Is_virtual_base( pose.residue( working_moving_res ) ) ) { //Check that the residue is not be already virtualized...
			utility_exit_with_message( "The base at " + string_of( working_moving_res ) + " is already virtualized!!" );
		}

		bool bulge_added = false;

		if ( allow_bulge_at_chainbreak_ ) {
			if ( delta_atr_score >= atr_cutoff_for_bulge ) {


				//Note that there is problem in that even after applying virtual_rna_residue, the chain break torsion potential is still scored for the chain_break torsions.
				//The should_score_torsion function in RNA_torsional_potential returns true (indicating that the score should be scored) if it finds a chain_break torsion,
				//even if this torsion contain virtual atoms.. May 4, 2010
				apply_virtual_rna_residue_variant_type( pose, working_moving_res, true );


				count_data_.bulge_at_chain_closure_count++;
				bulge_added = true;

				if ( verbose_ ){
					TR.Debug << "delta_atr " << delta_atr_score << " passes cutoff for bulge. " << atr_cutoff_for_bulge;
					TR.Debug << "  bulge = " << count_data_.bulge_at_chain_closure_count << "  both = " << count_data_.both_count << " tot = " << count_data_.tot_rotamer_count << std::endl;
				}

			} else {

				bulge_added = false;
				if ( verbose_ ) TR.Debug << "delta_atr " << delta_atr_score << " DOES NOT PASS cutoff for bulge " << atr_cutoff_for_bulge << std::endl;

			}
		}

		return bulge_added;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::Update_pose_data_list( std::string const & tag, utility::vector1< pose_data_struct2 > & pose_data_list, pose::Pose const & current_pose, Real const & current_score ) const{

		bool add_pose_to_list = false;

		add_pose_to_list = ( current_score < current_score_cutoff_ ) ? true: false; //May 12, 2010: updated code to be more robust (Parin S).


		//The order of evaluation of the two expression in the if statement is important!
		if ( add_pose_to_list ){

			if ( verbose_ ){
				TR.Debug << "tag = " << tag << " current_score_cutoff_ " << current_score_cutoff_ << " score = " << current_score;
			}

			pose_data_struct2 current_pose_data;
			current_pose_data.pose_OP = new pose::Pose;
			( *current_pose_data.pose_OP ) = current_pose;
			current_pose_data.score = current_score;
			current_pose_data.tag = tag;

			//if ( get_native_pose())  { //ACTUALLY DEPRECATED SINCE EARLY 2010!! Comment out on March 16, 2012
			//	setPoseExtraScores( *current_pose_data.pose_OP, "all_rms",
			//											core::scoring::rms_at_corresponding_heavy_atoms( *current_pose_data.pose_OP, *get_native_pose() ) );
			//}

			pose_data_list.push_back( current_pose_data );
			if ( verbose_ ) TR.Debug << " pose_data_list.size = " << pose_data_list.size() << std::endl;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	StepWiseRNA_ResidueSampler::Pose_selection_by_full_score( utility::vector1< pose_data_struct2 > & pose_data_list, pose::Pose & current_pose, std::string const & tag ){

		using namespace core::scoring;

		count_data_.full_score_count++;

		Real const current_score = ( *sampling_scorefxn_ )( current_pose );

		Update_pose_data_list( tag, pose_data_list, current_pose, current_score );

		if ( ( pose_data_list.size() == num_pose_kept_*multiplier_ ) ){
			std::sort( pose_data_list.begin(), pose_data_list.end(), sort_criteria );
			cluster_pose_data_list( pose_data_list );
			if ( pose_data_list.size() > num_pose_kept_ ){
				pose_data_list.erase( pose_data_list.begin() + num_pose_kept_, pose_data_list.end() );
				TR.Debug << "after erasing.. pose_data_list.size() = " << pose_data_list.size() << std::endl;
			} else{
				TR.Debug << "pose_data_list.size() = " << pose_data_list.size() << std::endl;
			}
		}

		return current_score;

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Dec 18, 2009...took off alot of optimization from this code since it is very fast (not rate limiting) anyways.
	void
	StepWiseRNA_ResidueSampler::cluster_pose_data_list( utility::vector1< pose_data_struct2 > & pose_data_list ){


		//Super hacky...take this out once solve the root atom problem
		if ( build_pose_from_scratch_ ) return; //no clustering at all;

		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const actually_moving_res = job_parameters_->actually_moving_res();

		utility::vector1< bool > pose_state_list( pose_data_list.size(), true );

		Size num_clustered_pose = 0;

		for ( Size i = 1; i <= pose_data_list.size(); i++ ){

			if ( pose_state_list[i] == true ){
				num_clustered_pose++;
				for ( Size j = i + 1; j <= pose_data_list.size(); j++ ){


					Real rmsd;
					if ( PBP_clustering_at_chain_closure_ && job_parameters_->gap_size() == 0 ){ //new option Aug 15, 2010..include both phosphates in rmsd calculation at chain_break
						rmsd =	 phosphate_base_phosphate_rmsd( ( *pose_data_list[i].pose_OP ), ( *pose_data_list[j].pose_OP ), actually_moving_res,  false /*ignore_virtual_atom*/ );
					} else{
						rmsd = suite_rmsd( ( *pose_data_list[i].pose_OP ), ( *pose_data_list[j].pose_OP ), actually_moving_res, Is_prepend, false /*ignore_virtual_atom*/ );
					}

					bool const same_pucker = Is_same_ribose_pucker( ( *pose_data_list[i].pose_OP ), ( *pose_data_list[j].pose_OP ), actually_moving_res );

					if ( rmsd < cluster_rmsd_ && ( same_pucker || !distinguish_pucker_ ) ){
						pose_state_list[j] = false;
						if ( verbose_ ) {
							TR.Debug << "rmsd = " << rmsd << "  pose " << pose_data_list[j].tag << " is a neighbor of pose " << pose_data_list[i].tag;
							TR.Debug << " same_pucker = "; Output_boolean( same_pucker, TR.Debug );
							print_ribose_pucker_state( " center_pucker = ", Get_residue_pucker_state( ( *pose_data_list[i].pose_OP ), actually_moving_res ), TR.Debug );
							print_ribose_pucker_state( " curr_pucker = ", Get_residue_pucker_state( ( *pose_data_list[j].pose_OP ), actually_moving_res ), TR.Debug );
							TR.Debug << std::endl;
						}
					}
				}
			}
		}


		utility::vector1< pose_data_struct2 > clustered_pose_data_list;

		for ( Size i = 1; i <= pose_data_list.size(); i++ ) {
			if ( pose_state_list[i] == true ){
				clustered_pose_data_list.push_back( pose_data_list[i] );
			}
		}

		pose_data_list = clustered_pose_data_list;

		//check if pose_data_list size is equal to or exceed num_pose_kept_. Important to get score_cutoff here right after clustering.
		if ( pose_data_list.size() >= num_pose_kept_ ){
			current_score_cutoff_ = pose_data_list[num_pose_kept_].score;
		} else{
			//keep on adding pose to list if there are still not enough clusters

			//current_score_cutoff_=99999999999.9999; //Feb 02, 2012 This might lead to server-test error at R47200
			current_score_cutoff_ = 999999.9; //Feb 02, 2012
		}
		////////////////////////////////////////////


	}
	///////////////////////////////////////////////////////////////////////////////
	//TEMPORARY
	bool
	StepWiseRNA_ResidueSampler::Full_atom_van_der_Waals_screening_REPLICATE( pose::Pose & current_pose_screen,
																																          Real const & base_rep_score,
																																          Real const & base_atr_score,
																																          Real & delta_atr_score,
																																          Real & delta_rep_score,
																																          Size const & gap_size,
																																          bool const & Is_internal ){

		using namespace core::scoring;

		if ( VDW_atr_rep_screen_ == false ) return true;

		bool close_chain = ( gap_size == 0 ) ? true: false;

		if ( close_chain && Is_internal ) return true; //Don't screen at all Mar 1, 2010

		( *atr_rep_screening_scorefxn_ )( current_pose_screen );

		EnergyMap const & energy_map = current_pose_screen.energies().total_energies();

		Real rep_score = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[scoring::fa_rep];
		Real atr_score = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[scoring::fa_atr];

		delta_rep_score = rep_score - base_rep_score;
		delta_atr_score = atr_score - base_atr_score;

		Real actual_rep_cutoff = rep_cutoff_; //defualt
		if ( close_chain ) actual_rep_cutoff = 10; //Parin's old parameter
		if ( close_chain && Is_internal ) actual_rep_cutoff = 200; //Bigger chunk..easier to crash...not using this right now.

		bool pass_rep_screen = false;

		if ( delta_rep_score < actual_rep_cutoff ){
			pass_rep_screen = true;
		}


		bool pass_atr_rep_screen = false;

		if ( close_chain ){
			pass_atr_rep_screen = pass_rep_screen;
		} else{
			if ( delta_atr_score < (  - 1 ) && ( delta_rep_score + delta_atr_score ) < 0 ) pass_atr_rep_screen = true;
		}


		if ( pass_atr_rep_screen ) {
			if ( verbose_ ) {
				TR.Debug << " rep = " << delta_rep_score << " atr = " << delta_atr_score;
				TR.Debug << "  stack_n = " << count_data_.base_stack_count << " pair_n = " << count_data_.base_pairing_count;
				TR.Debug << "  strict_pair_n = " << count_data_.strict_base_pairing_count;
				TR.Debug << "  centroid_n = " << count_data_.pass_base_centroid_screen;
				TR.Debug << "  bin_rep_n = " << count_data_.good_bin_rep_count;
				TR.Debug << "  atr_n = " << count_data_.good_atr_rotamer_count;
				TR.Debug << "  rep_n = " << count_data_.good_rep_rotamer_count;
				TR.Debug << "  both = " << count_data_.both_count << " tot = " << count_data_.tot_rotamer_count;
				TR.Debug << "  closable = " << count_data_.chain_closable_count;
				TR.Debug << "  non_clash_ribose = " << count_data_.non_clash_ribose;
				TR.Debug << std::endl;
			}
			return true;
		} else {
			return false;
		}

	}
	///////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_ResidueSampler::Full_atom_van_der_Waals_screening( pose::Pose & current_pose_screen,
																																Real const & base_rep_score,
																																Real const & base_atr_score,
																																Real & delta_rep_score,
																																Real & delta_atr_score,
																																Size const & gap_size,
																																bool const & Is_internal ){

		using namespace core::scoring;
		using namespace ObjexxFCL;

		if ( VDW_atr_rep_screen_ == false ) return true;

		bool close_chain = ( gap_size == 0 ) ? true: false;

		if ( close_chain && Is_internal ) return true; //Don't screen at all Mar 1, 2010

		( *atr_rep_screening_scorefxn_ )( current_pose_screen );

		EnergyMap const & energy_map = current_pose_screen.energies().total_energies();

		Real rep_score = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[scoring::fa_rep];
		Real atr_score = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[scoring::fa_atr];

		delta_rep_score = rep_score - base_rep_score;
		delta_atr_score = atr_score - base_atr_score;

		if ( delta_rep_score < (  - 0.01 ) ){
			std::string const message = "delta_rep_score = " + string_of( delta_rep_score ) + " rep_score = " + string_of( rep_score ) + " base_rep_score = " + string_of( base_rep_score );
			utility_exit_with_message( "delta_rep_score < (  - 0.01 ), " + message );
		}

		if ( delta_atr_score > (  + 0.01 ) ){
			std::string const message = "delta_atr_score = " + string_of( delta_atr_score ) + " atr_score = " + string_of( atr_score ) + " base_atr_score = " + string_of( base_atr_score );
			utility_exit_with_message( "delta_atr_score > (  + 0.01 ), " + message );
		}

		Real actual_rep_cutoff = rep_cutoff_; //defualt
		if ( close_chain ) actual_rep_cutoff = 10; //Parin's old parameter
		if ( Is_internal ) actual_rep_cutoff = 200; //Bigger chunk..easier to crash (before May 4 used to be (close_chain && Is_internal) actual_rep_cutoff=200

		bool pass_rep_screen = false;

		if ( delta_rep_score < actual_rep_cutoff ){
			pass_rep_screen = true;
			count_data_.good_rep_rotamer_count++;
		}

		if ( delta_atr_score < (  - 1 ) || close_chain ) count_data_.good_atr_rotamer_count++;

		bool pass_atr_rep_screen = false;

		if ( close_chain ){
			pass_atr_rep_screen = pass_rep_screen;
		} else if ( Is_internal ){
			if ( delta_atr_score < (  - 1 ) && ( delta_rep_score + delta_atr_score ) < ( actual_rep_cutoff - rep_cutoff_ ) ) pass_atr_rep_screen = true;
		} else{
			if ( delta_atr_score < (  - 1 ) && ( delta_rep_score + delta_atr_score ) < 0 ) pass_atr_rep_screen = true;
		}


		if ( pass_atr_rep_screen ) {
			//	if((delta_atr_score<(-1)) && ((delta_rep_score+delta_atr_score) < 200) ) { //This causes about 5times more pose to pass the screen (50,000 poses vs 10,000 poses)
			count_data_.both_count++;
			if ( verbose_ ) {
				TR.Debug << " rep = " << delta_rep_score << " atr = " << delta_atr_score;
				if ( combine_long_loop_mode_ && ( job_parameters_->gap_size() != 0 ) ) TR.Debug << " res_contact = " << count_data_.residues_contact_screen;
				TR.Debug << "  stack_n = " << count_data_.base_stack_count << " pair_n = " << count_data_.base_pairing_count;
				TR.Debug << "  strict_pair_n = " << count_data_.strict_base_pairing_count;
				TR.Debug << "  centroid_n = " << count_data_.pass_base_centroid_screen;
				TR.Debug << "  bin_rep_n = " << count_data_.good_bin_rep_count;
				TR.Debug << "  atr_n = " << count_data_.good_atr_rotamer_count;
				TR.Debug << "  rep_n = " << count_data_.good_rep_rotamer_count;
				TR.Debug << "  both = " << count_data_.both_count << " tot = " << count_data_.tot_rotamer_count << std::endl;
			}
			return true;
		} else {
			return false;
		}

	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::initialize_o2star_packer_task( core::pose::Pose const & pose ){

		utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();

		utility::vector1< core::Size > const O2star_pack_seq_num = get_surrounding_O2star_hydrogen( pose, working_moving_partition_pos, false /*verbose*/ );

		o2star_pack_task_ = create_standard_o2star_pack_task( pose, O2star_pack_seq_num );

	}

	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::initialize_o2star_green_packer( core::pose::Pose & pose )
	{
		using namespace protocols::simple_moves;
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		o2star_green_packer_ = new protocols::simple_moves::GreenPacker;

		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
		bool const root_partition = partition_definition( pose.fold_tree().root() );

		Size const nres = pose.total_residue();
		UserDefinedGroupDiscriminatorOP user_defined_group_discriminator( new UserDefinedGroupDiscriminator );
		utility::vector1< Size > group_ids;

		Size current_group = 0;
		Size spectator_group = 1;

		for ( Size i = 1; i <= nres; i++ ) {

			if ( partition_definition( i ) != root_partition ) {
				current_group = 0;
				TR.Debug << "GREENPACKER SAMPLER " << i << std::endl;
			} else {
				TR.Debug << "GREENPACKER SPECTATOR   " << i <<  " --> group " << spectator_group << std::endl;
			}
			group_ids.push_back( current_group );
		}

		user_defined_group_discriminator->set_group_ids( group_ids );
		o2star_green_packer_->set_scorefunction( *o2star_pack_scorefxn_ );
		o2star_green_packer_->set_group_discriminator( user_defined_group_discriminator );

		TaskFactoryOP task_factory( new TaskFactory );
		task_factory->push_back( new InitializeFromCommandline );
		task_factory->push_back( new RestrictToRepacking );
		task_factory->push_back( new IncludeCurrent );
		for ( Size i = 1; i <= nres; i++ ) {
			if ( !pose.residue( i ).is_RNA() ) continue;
			task_factory->push_back( new ExtraChiCutoff( i, 0 ) );
			task_factory->push_back( new ExtraRotamers( i, 4 /*ex4*/ ) );
		}

		o2star_green_packer_->set_task_factory( task_factory );
		o2star_green_packer_->set_reference_round_task_factory( task_factory );

		// This should also initialize rotamers, etc...
		o2star_green_packer_->apply( pose );
	}



	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::sample_o2star_hydrogen( core::pose::Pose & pose, core::pose::Pose & pose_with_original_HO2star_torsion ){

		using namespace core::id;
		using namespace core::conformation;


		//reset the HO2star torsion to its starting value as to prevent randomness due to conformations sampling order...
		copy_all_o2star_torsions( pose, pose_with_original_HO2star_torsion );

		//TR.Debug << "Packing 2'-OH ... ";
		if ( use_green_packer_ ) {
			o2star_green_packer_->apply( pose );
		} else {

			//problem with bulge variant -- need to initialize_o2star_packer_task each time.
			initialize_o2star_packer_task( pose );

			pack::rotamer_trials( pose, *o2star_pack_scorefxn_, o2star_pack_task_ );

		}

	}

	////////////////////////////////////////////////////////////////////////
	std::string //silly function to convert to real to string
	StepWiseRNA_ResidueSampler::create_torsion_value_string( core::Real const & torsion_value ) const{

		using namespace ObjexxFCL;

		std::string torsion_string = "";

		core::Real	const principal_torsion = numeric::principal_angle_degrees( torsion_value );

		Size const principal_torsion_SIZE = Size( std::abs( principal_torsion + 0.00001 ) ); //0.00001 is to prevent random ambiguity if the torsion decimal value is exactly .0000 Oct 12, 2010


		if ( principal_torsion > 0 ){
			torsion_string = "p" + lead_zero_string_of( principal_torsion_SIZE, 3 );
		} else{
			torsion_string = "n" + lead_zero_string_of( principal_torsion_SIZE, 3 );
		}

		return torsion_string;
	}

	////////////////////////////////////////////////////////////////////////
	std::string //silly function used for appending the rotamer value to the tag
	StepWiseRNA_ResidueSampler::create_rotamer_string( core::pose::Pose const & pose ) const{

		std::string rotamer_tag = "";

		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const moving_res(  job_parameters_->working_moving_res() );

		conformation::Residue const & five_prime_rsd = ( Is_prepend ) ? pose.residue( moving_res ): pose.residue( moving_res - 1 );
		conformation::Residue const & three_prime_rsd = ( Is_prepend ) ?  pose.residue( moving_res + 1 ) : pose.residue( moving_res );


		rotamer_tag.append( "_E" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 5  ) ) );
		rotamer_tag.append( "_Z" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 6  ) ) );
		rotamer_tag.append( "_A" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 1 ) ) );
		rotamer_tag.append( "_B" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 2 ) ) );
		rotamer_tag.append( "_G" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 3 ) ) );


		if ( Is_prepend ){
			rotamer_tag.append( "_D" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 4 ) ) );
			rotamer_tag.append( "_C" + create_torsion_value_string( five_prime_rsd.chi(  1 ) ) );

		} else{
			rotamer_tag.append( "_D" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 4 ) ) );
			rotamer_tag.append( "_C" + create_torsion_value_string( three_prime_rsd.chi( 1 ) ) );
		}

		return rotamer_tag;

	}
	////////////////////////////////////////////////////////////////////////

	std::string
	StepWiseRNA_ResidueSampler::create_tag( std::string const & prestring, Size const i ) const {

		using namespace ObjexxFCL;

		std::string tag = prestring;

		tag.append( "_" + lead_zero_string_of( i, 9 ) );

		return tag;
	}

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
  StepWiseRNA_ResidueSampler::set_fast( bool const & setting ){
    fast_ = setting;
		if ( fast_ ) num_pose_kept_ = 2;
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_medium_fast( bool const & setting ){
    medium_fast_ = setting;
		if ( medium_fast_ ) num_pose_kept_ = 20;
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
			set_fast( true );
		}

  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_verbose( bool const & setting ){
    verbose_ = setting;
  }
  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_perform_o2star_pack( bool const & setting ){
    perform_o2star_pack_ = setting;
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
			Output_data( silent_file_data, final_sampler_output_silent_file, pose_data_list_[n].tag, false, *( pose_data_list_[n].pose_OP ), get_native_pose(), job_parameters_ );
		}

	}

	//////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::set_base_centroid_screener( StepWiseRNA_BaseCentroidScreenerOP & screener ){
		base_centroid_screener_ = screener;
	}

	//////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::set_cluster_rmsd( Real const & setting ){
		cluster_rmsd_ = setting;
		TR.Debug << "Set cluster_rmsd to " << cluster_rmsd_ << std::endl;
	}

	//////////////////////////////////////////////////////////////////
	rotamer_sampler::rna::RNA_SuiteRotamerOP
	StepWiseRNA_ResidueSampler::setup_rotamer_sampler( Pose const & pose ) const {
		using namespace rotamer_sampler::rna;
		using namespace chemical::rna;
		using namespace pose::rna;
		/////Load in constants being used/////
		bool const Is_prepend( job_parameters_->Is_prepend() );
		bool const Is_internal( job_parameters_->Is_internal() );
		Size const gap_size = job_parameters_->gap_size();
		utility::vector1<Size> const & working_moving_suite_list(
				job_parameters_->working_moving_suite_list() );
		utility::vector1<Size> const & syn_chi_res(
			job_parameters_->working_force_syn_chi_res_list() );
		utility::vector1<Size> const & north_puckers(
			job_parameters_->working_force_north_ribose_list() );
		utility::vector1<Size> const & south_puckers(
			job_parameters_->working_force_south_ribose_list() );

		Size const n_rsd( working_moving_suite_list.size() );
		runtime_assert( n_rsd == 1 );

		Size const moving_suite( working_moving_suite_list[1] );

		/////Get the base and pucker state/////
		utility::vector1<bool> sample_sugar( 2, false );
		utility::vector1<Size> base_state( 2, WHATEVER );
		utility::vector1<Size> pucker_state( 2, WHATEVER );

		if ( build_pose_from_scratch_ || sample_both_sugar_base_rotamer_) {
			sample_sugar[1] = true;
			sample_sugar[2] = true;
		} else if ( !Is_internal  ) {
			if ( Is_prepend ) {
				sample_sugar[1] = true;
			} else {
				sample_sugar[2] = true;
			}
		}

		for ( Size i = 1; i <= 2; ++i ) {
			Size const curr_rsd( moving_suite + i - 1 );
			if ( sample_sugar[i] ) {
				bool is_north ( north_puckers.has_value( curr_rsd ) );
				bool is_south ( south_puckers.has_value( curr_rsd ) );
				bool is_syn ( syn_chi_res.has_value( curr_rsd ) );
				runtime_assert( !is_north || !is_south );
				if ( is_north ) pucker_state[i] = NORTH;
				if ( is_south ) pucker_state[i] = SOUTH;
				if ( !allow_syn_pyrimidine_ && !is_purine( pose.residue( curr_rsd ) ) ) {
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
		RNA_SuiteRotamerOP sampler = new RNA_SuiteRotamer( moving_suite,
				pucker_state[1], pucker_state[2], base_state[1], base_state[2] );
		sampler->set_skip_same_pucker( false ); //TODO Should enable this later
		sampler->set_idealize_coord( false ); //TODO Add and use the ideal-geo flag
		sampler->set_sample_nucleoside_lower( sample_sugar[1] );
		sampler->set_sample_nucleoside_upper( sample_sugar[2] );
		sampler->set_fast( fast_ );
		sampler->set_extra_epsilon( extra_epsilon_rotamer_ );
		sampler->set_extra_beta( extra_beta_rotamer_ );
		sampler->set_extra_chi(	extra_chi_ );
		sampler->set_random( choose_random_ );
		if ( gap_size == 0 && finer_sampling_at_chain_closure_ == true )
				sampler->set_bin_size( 10 );
		sampler->init();

		return sampler;
	}

}
}
}
