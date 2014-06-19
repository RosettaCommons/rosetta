// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Clusterer
/// @detailed
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


//////////////////////////////////
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Clusterer.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/stepwise/legacy/sampling/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/stepwise/legacy/sampling/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/stepwise/sampling/output_util.hh>
//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/rna/util.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
//#include <basic/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/pdb/pose_io.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <list>
#include <time.h>

using namespace core;
using namespace core::chemical::rna;
using core::Real;
using basic::T;

static basic::Tracer TR( "protocols.stepwise.sampling.rna.StepWiseRNA_Clusterer" );


namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {

// @brief Auto-generated virtual destructor
SlicedPoseWorkingParameters::~SlicedPoseWorkingParameters() {}


  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseRNA_Clusterer::StepWiseRNA_Clusterer( utility::vector1< std::string > const & silent_files_in )
  {
		initialize_parameters_and_input();
		input_->set_record_source( true );
		input_->filenames( silent_files_in ); //triggers read in of files, too.
  }

  StepWiseRNA_Clusterer::StepWiseRNA_Clusterer( std::string const & silent_file_in )
	{
		initialize_parameters_and_input();
		input_->set_record_source( true );

		utility::vector1< std::string > silent_files_;
		silent_files_.push_back( silent_file_in );
		input_->filenames( silent_files_ ); //triggers read in of files, too.
	}

  StepWiseRNA_Clusterer::StepWiseRNA_Clusterer( core::io::silent::SilentFileDataOP & sfd )
	{
		initialize_parameters_and_input();
		input_->set_silent_file_data( sfd ); // triggers reordering by energy and all that.
	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseRNA_Clusterer::~StepWiseRNA_Clusterer()
  {}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseRNA_Clusterer::initialize_parameters_and_input(){
		input_  = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		input_->set_order_by_energy( true );

		//max_decoys_ = 9999999999; Feb 02, 2012; This lead to server-test error at R47198
		//score_diff_cut_ = 1000000000.0; Feb 02, 2012; This might lead to server-test error at R47198
		max_decoys_ = 999999; //Feb 02, 2012;
		score_diff_cut_ = 100000.0; //Feb 02, 2012;
		perform_score_diff_cut_ = false; //Jan 23, 2012: I rarely use score_diff_cut_ in SWA RNA, so PLEASE leave this false as DEFAULT to be safe!

		whole_struct_cluster_radius_ = 0.5;
		suite_cluster_radius_ = 999.99;
		loop_cluster_radius_ = 999.99;

		rename_tags_ = false;
		working_parameters_exist_ = false;
		distinguish_pucker_ = true;
		add_lead_zero_to_tag_ = false; //For easier pdb selection in pymol
		quick_alignment_ = false; //new option May 29, 2010...speed up clustering code...however only work if alignment residues are fixed res
		align_only_over_base_atoms_ = true; //Set to true for backward compatibility. Add option in Aug 20, 2011
		optimize_memory_usage_ = false;
		verbose_ = true;
		keep_pose_in_memory_ = true; //Can save memory at the expense of speed by not keeping the pose
		keep_pose_in_memory_hydrid_ = true; //basically keep pose in memory until memory runs out!
		max_memory_pose_num_ = 0;
		two_stage_clustering_ = false; //Cluster is two stage using triangle inequaility to speed up clustering. Need keep_pose_in_memory mode==false or else code will be too slow.
		use_triangle_inequality_ = false; //This is turned on during the second stage of the two_stage_clustering mode;
		PBP_clustering_at_chain_closure_ = false;
		quick_alignment_pose_is_intialized_ = false;
		skip_clustering_ = false; //I know this is weird...basically, this is for using clusterer to recalculate rmsd.
		perform_VDW_rep_screen_ = false; //March 20, 2011
		perform_filters_ = false; //June 20, 2011.
		VDW_rep_screen_info_.clear(); //make sure that this is empty, Marc 21, 2011
		full_length_loop_rmsd_clustering_ = false;
		ignore_FARFAR_no_auto_bulge_tag_ = false; //Sept 06, 2011..for post-processing.
		ignore_FARFAR_no_auto_bulge_parent_tag_ = false; //Sept 06, 2011..for post-processing.
		ignore_unmatched_virtual_res_ = false; //Sept 07, 2011...for post-processing
		output_pdb_ = false; //Sept 24, 2011
		min_num_south_sugar_filter_ = 0; //Oct 02, 2011
		rsd_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );
	}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseRNA_Clusterer::cluster()
	{
		using namespace core::scoring;
		using namespace core::import_pose::pose_stream;
		using namespace core::chemical;
		using namespace core::pose;

		clock_t const time_start( clock() );
		output_title_text( "StepWiseRNA_Clusterer::cluster()", TR );

		output_boolean( "verbose_ = ", verbose_, TR ); TR << std::endl;
		output_boolean( "skip_clustering_ = ", skip_clustering_, TR ); TR << std::endl;

		if ( skip_clustering_ ){

			//Commented out Dec 11, 2011. CHANGE TO SET THIS FROM COMMAND_LINE (See for example SWA_cluster.py!
			//TR << "skip_clustering==true --> set keep_pose_in_memory to false" << std::endl;
			//keep_pose_in_memory_=false;

			TR << "skip_clustering == true --> set keep_pose_in_memory_hydrid_ to false" << std::endl; //HACKY!
			keep_pose_in_memory_hydrid_ = false; //HACKY!
		}


		///March 20, 2011/////////////////////
		if ( ( perform_VDW_rep_screen_ == true ) && ( VDW_rep_screen_info_.size() == 0 ) ){
			TR << "User pass in perform_VDW_rep_screen_ == true but VDW_rep_screen_info_.size() == 0" << std::endl;
			TR << "Override and set perform_VDW_rep_screen_ to false" << std::endl;
			perform_VDW_rep_screen_ = false;
		}

		output_boolean( "perform_VDW_rep_screen_ = ", perform_VDW_rep_screen_, TR ); TR << std::endl;
		output_boolean( "perform_filters_ = ", perform_filters_, TR ); TR << std::endl; ///June 14, 2011 Perform other filters aside from VDW_rep_screen.

		if ( perform_filters_ && ( skip_clustering_ == false ) ) utility_exit_with_message( "perform_filters_ but skip_clustering_ == false" );
		if ( perform_VDW_rep_screen_ && ( skip_clustering_ == false ) )  utility_exit_with_message( "perform_VDW_rep_screen_ but skip_clustering_ == false" );
		/////////////////////////////////////

		TR << "suite_cluster_radius_ = " << suite_cluster_radius_ << std::endl;
		TR << "loop_cluster_radius_ = " << loop_cluster_radius_ << std::endl;
		output_boolean( "working_parameters_exist_ = ", working_parameters_exist_, TR ); TR << std::endl;
		output_boolean( "quick_alignment_ = ", quick_alignment_, TR ); TR << std::endl;
		output_boolean( "align_only_over_base_atoms_ = ", align_only_over_base_atoms_, TR ); TR << std::endl;
		output_boolean( "two_stage_clustering_ = ", two_stage_clustering_, TR ); TR << std::endl;
		output_boolean( "keep_pose_in_memory_ = ", keep_pose_in_memory_, TR ); TR << std::endl;
		output_boolean( "keep_pose_in_memory_hydrid_ = ", keep_pose_in_memory_hydrid_, TR ); TR << std::endl;
		output_boolean( "optimize_memory_usage_( by slicing out fixed region of the pose ) = ", optimize_memory_usage_, TR ); TR << std::endl;
		output_boolean( "distinguish_pucker_ = ", distinguish_pucker_, TR ); TR << std::endl;
		output_boolean( "add_lead_zero_to_tag_ = ", add_lead_zero_to_tag_, TR ); TR << std::endl;
		output_boolean( "PBP_clustering_at_chain_closure_ = ", PBP_clustering_at_chain_closure_, TR ); TR << std::endl;
		output_boolean( "full_length_loop_rmsd_clustering_ = ", full_length_loop_rmsd_clustering_, TR ); TR << std::endl;
		output_boolean( "ignore_FARFAR_no_auto_bulge_tag_ = ", ignore_FARFAR_no_auto_bulge_tag_, TR ); TR << std::endl;
		output_boolean( "ignore_FARFAR_no_auto_bulge_parent_tag_ = ", ignore_FARFAR_no_auto_bulge_parent_tag_, TR ); TR << std::endl;
		output_boolean( "ignore_unmatched_virtual_res_ = ", ignore_unmatched_virtual_res_, TR ); TR << std::endl;

		TR << "max_decoys_ = " << max_decoys_ << std::endl;
		TR << "score_diff_cut_ = " << score_diff_cut_ << std::endl;
		output_boolean( "perform_score_diff_cut_ = ", perform_score_diff_cut_, TR ); TR << std::endl;

		//////////basic initialization///////////////
		pose_output_list_.clear();
		tag_output_list_.clear();
		silent_struct_output_list_.clear();
		/////////////////////////////////////////////

		if ( optimize_memory_usage_ ){
			if ( !working_parameters_exist_ ) utility_exit_with_message( "optimize_memory_usage = True but working_parameters_exist_ = False!" );
			sliced_pose_job_params_.setup( working_parameters_ );
		}

		if ( ignore_FARFAR_no_auto_bulge_tag_ || ignore_FARFAR_no_auto_bulge_parent_tag_ ){
			create_tags_map();
		}

		if ( quick_alignment_ ) initialize_quick_alignment_pose();

		if ( perform_VDW_rep_screen_ ) initialize_VDW_rep_checker();

		initialize_max_memory_pose_num();

		if ( skip_clustering_ ){
			create_silent_file_and_tag_list();
		} else if ( two_stage_clustering_ ){
			two_stage_clustering();
		} else {
			do_some_clustering();
		}


		if ( tag_output_list_.size() != silent_struct_output_list_.size() ) utility_exit_with_message( "tag_output_list_.size() != silent_struct_output_list_.size()" );

		if ( ( keep_pose_in_memory_ == true ) && ( keep_pose_in_memory_hydrid_ == false ) ){
			if ( pose_output_list_.size() != tag_output_list_.size() ) utility_exit_with_message( "pose_output_list_.size() != tag_output_list_.size()" );
		}


		TR << "Final cluster_pose_list size = " << silent_struct_output_list_.size() << std::endl;
		TR << "Total clustering time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::initialize_max_memory_pose_num(){

		using namespace core::pose;
		using namespace ObjexxFCL;

		clock_t const time_start( clock() );


		pose::Pose first_pose;
		pose::Pose first_pose_before_slicing;

		Size num_silent_struct = 0;
		bool found_valid_struct = false;

		input_->reset(); //reset the silentfile stream to the beginning..


		//get the first pose in the silent_file_stream.
		while ( input_->has_another_pose() ) {
			num_silent_struct++;

			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );

			if ( found_valid_struct == false ){
				PoseOP pose_op( new Pose );
				silent_struct->fill_pose( *pose_op, *rsd_set_ );

				std::string const & tag( silent_struct->decoy_tag() );

				if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( (*pose_op), tag ) == true ) continue;

				first_pose_before_slicing = (*pose_op);

				if ( optimize_memory_usage_ ) (*pose_op) = sliced_pose_job_params_.create_sliced_pose(*pose_op);

				//OK found a valid (non-messed) pose. Will use this pose as the "global" quick alignment pose_
				first_pose = (*pose_op);
				found_valid_struct = true;
			}
		}

		input_->reset(); //reset the silentfile stream to the beginning..

		Size const total_res_before_slicing = first_pose_before_slicing.total_residue();

		Size const total_res = first_pose.total_residue();

		//OK, one example of crash due to insufficient memory:
		//Building region 11_4 of J5/J5a hinge 2r8s
		//12 nucleotides pose, 4G, 2C, 3U and 3A
		//209,155 silent struct  280 large clusters pose, and 9889 normal pose.
		//the size of the silent_file is: 1.4G REGION_11_4/start_from_region_11_2_sample_filtered.out
		//So 6.69 KB. 0.55 KB per nucleotide
		//Memory limit on Biox is
 		//MEMLIMIT
 		//4000000 KB, 4G.
		//So mememory used to store pose is 2.4G for (9889+280=10169 pose)
		//236 KB per pose. 19 KB per nucleotide.

		//Consistency check:
		//Finished reading 164975 structures from REGION_4_8/start_from_region_5_8_sample_filtered.out
		//553M    REGION_4_8/start_from_region_5_8_sample_filtered.out
		//This is 1ZIH 4-8 is res 5-9, (5,6,7,8,9)-> 5 res. -> 0.6 KB per nucleotide!!
		//Finished reading 247508 structures from REGION_4_8/start_from_region_4_7_sample_filtered.out
		//277M    REGION_4_8/start_from_region_4_7_sample_filtered.out
		//->0.2 KB per nucleotide... WHY DOESN THE VALUE FLUCTUATE SO MUCH??

		//4,000,000=(max_memory_pose_num_)*19*(total_res) + (num_silent_struct)*(0.078)*(total_res)
		Real const total_memory = 4000000;

		Real const memory_taken_by_silent_struct = ( num_silent_struct*0.55*total_res_before_slicing );

		if ( memory_taken_by_silent_struct > total_memory ){
			max_memory_pose_num_ = 0;
			TR << "memory_taken_by_silent_struct ( " << memory_taken_by_silent_struct << " ) > specified_total_memory( " << total_memory << " )" << std::endl;
		} else{
			max_memory_pose_num_ = int( 0.7*( ( total_memory - memory_taken_by_silent_struct )/( 19*total_res ) ) );  //0.7 is to be on the safe side
		}

		TR << "--------------StepWiseRNA_Clusterer::initialize_max_memory_pose_num----------" << std::endl;
		output_boolean( "optimize_memory_usage_ ( by slicing ) = ", optimize_memory_usage_, TR ); TR << std::endl;
		TR << "first_pose total_res ( before_slicing ) = " << total_res_before_slicing << std::endl;
		TR << "first_pose total_res ( already account for slicing ) = " << total_res << std::endl;
		TR << "num_silent_struct = " << num_silent_struct << std::endl;
		TR << "memory_taken_by_silent_struct = " << memory_taken_by_silent_struct << std::endl;
		TR << "max_memory_pose_num_ = " << max_memory_pose_num_ << std::endl;
		TR << "time in function = " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
		TR << "--------------StepWiseRNA_Clusterer::initialize_max_memory_pose_num----------" << std::endl;

	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::initialize_VDW_rep_checker(){

		using namespace core::pose;
		using namespace ObjexxFCL;

		if ( !working_parameters_exist_ ) utility_exit_with_message( "perform_VDW_rep_screen_ = True but working_parameters_exist_ = False!" );

		if ( optimize_memory_usage_ ) utility_exit_with_message( "perform_VDW_rep_screen_ = True and optimize_memory_usage_ = True!" );

		if ( working_parameters_->is_simple_full_length_job_params() == true ) utility_exit_with_message( "working_parameters_->is_simple_full_length_job_params() == true!" );

		input_->reset(); //reset the silentfile stream to the beginning..

		//get the first pose in the silent_file_stream.
		while ( input_->has_another_pose() ) {

			PoseOP pose_op( new Pose );
			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );
			silent_struct->fill_pose( *pose_op, *rsd_set_ );

			std::string const & tag( silent_struct->decoy_tag() );

			if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( (*pose_op), tag ) == true ) continue;

			if ( optimize_memory_usage_ ) (*pose_op) = sliced_pose_job_params_.create_sliced_pose(*pose_op);

			//OK found a valid (non-messed) pose. Will use this pose as the "global" quick alignment pose_

			user_input_VDW_bin_checker_->setup_using_user_input_VDW_pose( VDW_rep_screen_info_, (*pose_op), working_parameters::StepWiseWorkingParametersCOP( working_parameters_ ) );

			break;

		}

		input_->reset(); //reset the silentfile stream to the beginning..

	}

	/////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_Clusterer::initialize_quick_alignment_pose(){

		using namespace core::pose;
		using namespace ObjexxFCL;

		if ( !working_parameters_exist_ ) utility_exit_with_message( "quick_alignment_ = True but working_parameters_exist_ = False!" );


		//OK first check that it valid to use the quick_alignment_pose mode... in this mode, all alignment must be fixed res..However, this check itself is not enough to gaurantee that quicj_alignemnt_mode will work. Another requirement is that all the residue in working_best_alignment must be fixed in space with respect to each other. I try to ensure that this is always the case by making sure that every residues in the working_best_alignment is in the root_partition. (See StepWiseWorkingParametersSetup.cc)

		utility::vector1< core::Size > const working_best_alignment = working_parameters_->working_best_alignment();
		utility::vector1< core::Size > const working_fixed_res = working_parameters_->working_fixed_res();

		for ( Size n = 1; n <= working_best_alignment.size(); n++ ){
			Size const seq_num = working_best_alignment[n];

			if ( working_fixed_res.has_value( seq_num ) == false ) {

				output_seq_num_list( "working_best_alignment = ", working_best_alignment, TR, 30 );
				output_seq_num_list( "working_fixed_res = ", working_fixed_res, TR, 30 );

				utility_exit_with_message( "quick_alignment_mode is true. However: seq_num " + string_of( seq_num ) + " is a element of working_best_alignment BUT not a element of working_fixed_res " );

			}
		}


		input_->reset(); //reset the silentfile stream to the beginning..

		quick_alignment_pose_is_intialized_ = true;

		//get the first pose in the silent_file_stream.
		while ( input_->has_another_pose() ) {


			PoseOP pose_op( new Pose );
			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );
			silent_struct->fill_pose( *pose_op, *rsd_set_ );

			std::string const & tag( silent_struct->decoy_tag() );

			if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( (*pose_op), tag ) == true ) continue;

			if ( optimize_memory_usage_ ) (*pose_op) = sliced_pose_job_params_.create_sliced_pose(*pose_op);

			//OK found a valid (non-messed) pose. Will use this pose as the "global" quick alignment pose_

			quick_alignment_pose_ = (*pose_op);
			quick_alignment_tag_ = tag;
			TR << "found quick alignment_pose, tag = " << tag << std::endl;

			break;

		}

		if ( output_pdb_ ) quick_alignment_pose_.dump_pdb( "quick_alignment_pose_" + quick_alignment_tag_ + ".pdb" );

		input_->reset(); //reset the silentfile stream to the beginning..


	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::align_to_quick_alignment_pose( core::pose::Pose & pose, std::string const & tag ) const {

		using namespace core::pose;

		if ( quick_alignment_pose_is_intialized_ == false ) utility_exit_with_message( "quick_alignment_pose_is_intialized_ == false" );

		utility::vector1< core::Size > const & alignment_res =  get_act_alignment_res();

		align_poses( pose, tag, quick_alignment_pose_, "quick_alignment_tag_" + quick_alignment_tag_, alignment_res, align_only_over_base_atoms_ );

	}

  //////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_Clusterer::two_stage_clustering(){

		output_title_text( "Enter two_stage_clustering function", TR );

		Real const whole_struct_cluster_radius_actual = whole_struct_cluster_radius_;
		Real const loop_cluster_radius_actual  = loop_cluster_radius_;
		Real const suite_cluster_radius_actual = suite_cluster_radius_;
		bool const keep_pose_in_memory_actual = keep_pose_in_memory_;
		bool const keep_pose_in_memory_hydrid_actual = keep_pose_in_memory_hydrid_;

		whole_struct_cluster_radius_ = 2.0; //hard code
		loop_cluster_radius_ = 2.0; //hard_code
		suite_cluster_radius_ = 999; //hard_code ...no suite_cluster...
		keep_pose_in_memory_ = true;
		keep_pose_in_memory_hydrid_ = false;
		use_triangle_inequality_ = false;

		output_title_text( "First stage: large RMSD clustering", TR );

		input_->reset(); //reset the silentfile stream to the beginning..
		do_some_clustering();

		large_cluster_pose_list_ = pose_output_list_;
		pose_output_list_.clear();
		tag_output_list_.clear();
		silent_struct_output_list_.clear();

		//////////////////////////////////////////////////////


		input_->reset(); //reset the silentfile stream to the beginning..
		create_large_cluster_centers_member_list();

		//////////////////////////////////////////////////////
		output_title_text( "Second stage: Actual clustering", TR );


		//Reset to actual (user specified) value
		whole_struct_cluster_radius_ = whole_struct_cluster_radius_actual;
		loop_cluster_radius_  = loop_cluster_radius_actual;
		suite_cluster_radius_ = suite_cluster_radius_actual;
		keep_pose_in_memory_ = keep_pose_in_memory_actual;
		keep_pose_in_memory_hydrid_ = keep_pose_in_memory_hydrid_actual;

		use_triangle_inequality_ = true;

		input_->reset(); //reset the silentfile stream to the beginning..
		do_some_clustering();

	}

	/////////////////////////////////////////////////////////////////////
	//The member is both dimension of the vector are sorted so that lowest energy appear first.

	void
	StepWiseRNA_Clusterer::create_large_cluster_centers_member_list(){

		using namespace core::pose;


		output_title_text( "create_large_cluster_centers_member_list", TR );
		clock_t const time_start( clock() );

		cluster_centers_neighbor_list_.clear();

		utility::vector1< Cluster_Member > empty_vector;
		cluster_centers_neighbor_list_.assign( large_cluster_pose_list_.size(), empty_vector );


		utility::vector1< core::Size > const & alignment_res =  get_act_alignment_res();
		utility::vector1 < core::Size > const & calc_rms_res = get_act_calc_rms_res();
		std::map< core::Size, core::Size > const & full_to_sub = get_act_full_to_sub();
		std::map< core::Size, bool > const & is_prepend_map = get_act_is_prepend_map();

		Size input_ID = 0;

		Real last_cluster_center_score( 0.0 ); //does slicing the pose change the score?
		getPoseExtraScores( *( large_cluster_pose_list_[large_cluster_pose_list_.size()] ), "score", last_cluster_center_score );

		while ( input_->has_another_pose() ) {

			input_ID++; //count messed up poses as well

			if ( ( input_ID % 1000 ) == 0 ){
				TR << "input_ID = " << input_ID << " time taken so far = " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
			}


			PoseOP pose_op( new Pose );
			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );

			if ( pass_FARFAR_no_auto_bulge_filter( silent_struct ) == false ) continue;

			silent_struct->fill_pose( *pose_op, *rsd_set_ ); //umm is the pose still connected to the silent_struct? Meaning that if the pose change, does the silent struct get changed as well? Apr 23, 2010 Parin.

			std::string const & tag( silent_struct->decoy_tag() );

			//Hacky thing. Ignore messed structure until we find a fix
			if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( (*pose_op), tag ) == true ) continue;

			if ( optimize_memory_usage_ ) (*pose_op) = sliced_pose_job_params_.create_sliced_pose(*pose_op);

			//need to align pose_op to the global alignment pose.
			if ( quick_alignment_ ) align_to_quick_alignment_pose( *pose_op, tag );

			Real score( 0.0 ); //does slicing the pose change the score?
			getPoseExtraScores( *pose_op, "score", score );

			if ( score > ( last_cluster_center_score + 0.001 ) ) break; //Exclude bad score poses that will never to be part of the final output_pose_list.

			for ( Size n = 1; n <= large_cluster_pose_list_.size(); n++ ){

				pose::Pose const & cluster_center_pose = *( large_cluster_pose_list_[n] );

				if ( quick_alignment_ == false ) align_poses( *pose_op, "current_pose", cluster_center_pose, "large_cluster_center", alignment_res, align_only_over_base_atoms_ );

				Real const RMSD = rmsd_over_residue_list( *pose_op, cluster_center_pose, calc_rms_res, full_to_sub, is_prepend_map, false );

				if ( RMSD < ( loop_cluster_radius_*1.5 ) ){ //A neigbor/member of this cluster_center
					Cluster_Member member;
					member.ID = input_ID;
					member.RMSD = RMSD;
					member.score = score;
					cluster_centers_neighbor_list_[n].push_back( member );
				}
			}
		}

		TR << "check large_cluster_pose_list_ member size " << std::endl;
		for ( Size n = 1; n <= large_cluster_pose_list_.size(); n++ ){
			TR << "cluster center " << n << " has " << cluster_centers_neighbor_list_[n].size() << " members " << std::endl;
		}


	}


	/////////////////////////////////////////////////////////////////////


	void
	StepWiseRNA_Clusterer::create_silent_file_and_tag_list(){

		using namespace core::pose;
		using namespace ObjexxFCL;


		output_title_text( "StepWiseRNA_Clusterer::create_silent_file_and_tag_list()", TR );

		input_->reset(); //reset the silentfile stream to the beginning..

		tag_output_list_.clear();
		silent_struct_output_list_.clear();
		pose_output_list_.clear();

		utility::vector1 < core::Size > working_global_sample_res_list;
		utility::vector1 < core::Size > working_filter_virtual_res_list;


		if ( perform_VDW_rep_screen_ || perform_filters_ ){

			if ( working_parameters_exist_ == false ) utility_exit_with_message( "( perform_VDW_rep_screen_ || perform_filters_ ) but working_parameters_exist_ == false!" );

			working_global_sample_res_list = working_parameters_->working_global_sample_res_list();
			working_filter_virtual_res_list = apply_full_to_sub_mapping( filter_virtual_res_list_, working_parameters_ );

			output_seq_num_list( "filter_virtual_res_list_ = ", filter_virtual_res_list_, TR, 50 );
			output_seq_num_list( "working_filter_virtual_res_list = ", working_filter_virtual_res_list, TR, 50 );
			output_seq_num_list( "working_global_sample_res_list = ", working_global_sample_res_list, TR, 50 );
			TR << "min_num_south_sugar_filter_ = " << min_num_south_sugar_filter_ << std::endl;
		}

		Size input_ID = 0;

		bool filter_verbose = true;

		while ( input_->has_another_pose() ) {
			input_ID++;

			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );
			std::string const & tag( silent_struct->decoy_tag() );

			if ( pass_FARFAR_no_auto_bulge_filter( silent_struct ) == false ) continue;

			if ( perform_VDW_rep_screen_ || perform_filters_ ){

				PoseOP pose_op( new Pose );
				silent_struct->fill_pose( *pose_op, *rsd_set_ ); //umm is the pose still connected to the silent_struct?

				if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( (*pose_op), tag ) == true ) continue;

				///Jan 12, 2012:Consistency check://///
				if ( working_parameters_exist_ ){
					if ( (*pose_op).total_residue() != working_parameters_->working_sequence().size() ){
						utility_exit_with_message( "(*pose_op).total_residue() = ( " + string_of( (*pose_op).total_residue() ) + " ) != ( " + string_of( working_parameters_->working_sequence().size() ) + " ) = working_parameters_working_sequence().size()" );
					}
				}
				////////////////////////////////////////

				if ( perform_VDW_rep_screen_ ){

					if ( user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() != true ){
						utility_exit_with_message( "user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() != true" );
					}

					bool const pass_VDW_rep_screen = user_input_VDW_bin_checker_->VDW_rep_screen_with_act_pose( (*pose_op), working_global_sample_res_list, false /*local verbose*/ );
					if ( pass_VDW_rep_screen == false ){
						if ( filter_verbose ) TR << "tag = " << tag << " fail VDW_rep_screen! " << std::endl;
						continue;
					}
				}

				if ( perform_filters_ ){
					bool pass_filter = true;

					utility::vector1< core::Size > const & force_north_sugar_list = working_parameters_->working_force_north_sugar_list();
					utility::vector1< core::Size > const & force_south_sugar_list = working_parameters_->working_force_south_sugar_list();
					utility::vector1< core::Size > const & force_syn_chi_res_list = working_parameters_->working_force_syn_chi_res_list();

					for ( Size n = 1; n <= force_north_sugar_list.size(); n++ ){
						Size const seq_num = force_north_sugar_list[n];
						if ( (*pose_op).residue( seq_num ).has_variant_type( "BULGE" ) ) continue;
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RIBOSE" ) ) continue;
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
						if ( get_residue_pucker_state( (*pose_op), seq_num ) != NORTH ){
							pass_filter = false;
							if ( filter_verbose ) TR << "pose = " << tag << " doesn't have north_sugar at seq_num = " << seq_num << std::endl;
						}
					}

					for ( Size n = 1; n <= force_south_sugar_list.size(); n++ ){
						Size const seq_num = force_south_sugar_list[n];
						if ( (*pose_op).residue( seq_num ).has_variant_type( "BULGE" ) ) continue;
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RIBOSE" ) ) continue;
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
						if ( get_residue_pucker_state( (*pose_op), seq_num ) != SOUTH ){
							 pass_filter = false;
							if ( filter_verbose ) TR << "pose = " << tag << " doesn't have south_sugar at seq_num = " << seq_num << std::endl;
						}
					}

					for ( Size n = 1; n <= force_syn_chi_res_list.size(); n++ ){
						Size const seq_num = force_syn_chi_res_list[n];
						if ( (*pose_op).residue( seq_num ).has_variant_type( "BULGE" ) ) continue;
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RIBOSE" ) ) continue;
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
						if ( get_residue_base_state( (*pose_op), seq_num ) != SYN ){
						 	pass_filter = false;
							if ( filter_verbose ) TR << "pose = " << tag << " doesn't have syn_chi at seq_num = " << seq_num << std::endl;
						}
					}

					for ( Size n = 1; n <= working_filter_virtual_res_list.size(); n++ ){
						Size const seq_num = working_filter_virtual_res_list[n];
						if ( (*pose_op).residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == false ){
							pass_filter = false;
							if ( filter_verbose ) TR << "pose = " << tag << " doesn't have virtual_rna_residue variant_type at seq_num = " << seq_num << std::endl;
						}
					}

					if ( min_num_south_sugar_filter_ != 0 ){
						Size num_south_sugar = 0;
						for ( Size n = 1; n <= working_global_sample_res_list.size(); n++ ){
							Size const seq_num = working_global_sample_res_list[n];
							if ( get_residue_pucker_state( (*pose_op), seq_num ) == SOUTH ){
								num_south_sugar += 1;
							}
						}
						//if(filter_verbose) TR << "pose= " << tag << " have " << num_south_sugar << " south_pucker_sugar." << std::endl;
						if ( num_south_sugar < min_num_south_sugar_filter_ ) pass_filter = false;
					}

					if ( pass_filter == false ) continue;

				}

			}

			if ( verbose_ ) TR << "Adding " << tag << " ID = " << input_ID << std::endl;

			tag_output_list_.push_back(  tag );
			silent_struct_output_list_.push_back(  silent_struct  );


			if ( keep_pose_in_memory_ ){
				if ( ( keep_pose_in_memory_hydrid_ == false ) || ( pose_output_list_.size() < max_memory_pose_num_ )  ){
					PoseOP localized_pose_op( new Pose );
					silent_struct->fill_pose( *localized_pose_op, *rsd_set_ );
					pose_output_list_.push_back(  localized_pose_op );
				}
			}

		}

		if ( perform_VDW_rep_screen_ || perform_filters_ ){
			TR << tag_output_list_.size() << " out of " << input_ID << " poses pass the filters." << std::endl;
		}

		input_->reset(); //reset the silentfile stream to the beginning..

		output_title_text( "", TR );


	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::do_some_clustering() {

		using namespace core::pose;
		using namespace ObjexxFCL;

		clock_t const time_start( clock() );

		input_->reset(); 				 					//Dec 11, 2011.
		tag_output_list_.clear();	 					//Dec 11, 2011.
		silent_struct_output_list_.clear(); //Dec 11, 2011.
		pose_output_list_.clear(); 					//Dec 11, 2011.

		if ( use_triangle_inequality_ ) all_pose_to_output_pose_ID_map_.clear();

		bool is_first_pose = true;
		bool score_min_defined = false;
		Real score_min = 0.0;

		Size num_pose_clustered = 0;
		Size input_ID = 0; //this count messed up pose where as num_pose_clustered doesn't
		while ( input_->has_another_pose() ){

			input_ID++; //count messed up poses as well
			if ( use_triangle_inequality_ ) all_pose_to_output_pose_ID_map_.push_back( 0 ); //If is a output_pose, ID value will be updated below

			PoseOP pose_op( new Pose );
			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );

			if ( pass_FARFAR_no_auto_bulge_filter( silent_struct ) == false ) continue;

			silent_struct->fill_pose( *pose_op, *rsd_set_ );

			Real score( 0.0 );
			getPoseExtraScores( *pose_op, "score", score );

			if ( score_min_defined == false ){
				score_min = score;
				score_min_defined = true;
			}

			if ( perform_score_diff_cut_ && ( score > ( score_min + score_diff_cut_ ) ) ) break;

			std::string const & tag( silent_struct->decoy_tag() );
			TR << "CHECKING " << tag << " with score " << score << " ( score_min = " << score_min << " ) against list of size " << silent_struct_output_list_.size();
			TR << " Num_pose_clustered so far " << num_pose_clustered << std::endl;

			//Hacky thing. Ignore messed structure until we find a fix...ideally should just remove messed up pose from input_ at beginning of the Class.
			if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( (*pose_op), tag ) == true ) continue;

			///Jan 12, 2012:Consistency check://///
			if ( working_parameters_exist_ ){
				if ( working_parameters_ -> add_virt_res_as_root() ) {
					if ( (*pose_op).total_residue() - 1 != working_parameters_->working_sequence().size() ){
						utility_exit_with_message( "(*pose_op).total_residue() = ( " + string_of( (*pose_op).total_residue() ) + " ) != ( " + string_of( working_parameters_->working_sequence().size() ) + " ) = working_parameters_working_sequence().size()" );
					}
				} else {
					if ( (*pose_op).total_residue() != working_parameters_->working_sequence().size() ){
						utility_exit_with_message( "(*pose_op).total_residue() = ( " + string_of( (*pose_op).total_residue() ) + " ) != ( " + string_of( working_parameters_->working_sequence().size() ) + " ) = working_parameters_working_sequence().size()" );
					}
				}
			}
			////////////////////////////////////////

			if ( optimize_memory_usage_ ) (*pose_op) = sliced_pose_job_params_.create_sliced_pose(*pose_op);

			if ( is_first_pose ){
				first_pose_ = (*pose_op);
				is_first_pose = false;
			}


			/////////////////////////////////////////////////////////


			bool const OK = check_for_closeness( pose_op, tag );

			if ( OK )  {
				TR << "ADDING " << tag << std::endl;

				tag_output_list_.push_back(  tag );

				if ( keep_pose_in_memory_ == true ){
					if ( ( keep_pose_in_memory_hydrid_ == false ) || ( pose_output_list_.size() < max_memory_pose_num_ )  ){

					 	pose_output_list_.push_back(  pose_op );

					}
				}

				silent_struct_output_list_.push_back(  silent_struct  );

				if ( use_triangle_inequality_ ) {

					if ( input_ID > all_pose_to_output_pose_ID_map_.size() ){
						utility_exit_with_message( "input_ID > all_pose_to_output_pose_ID_map_.size(), input_ID = " + string_of( input_ID ) + ", all_pose_to_output_pose_ID_map_.size() = " + string_of( all_pose_to_output_pose_ID_map_.size() ) );
					}

					all_pose_to_output_pose_ID_map_[input_ID] = silent_struct_output_list_.size();
				}

				if ( silent_struct_output_list_.size() >= max_decoys_ ) break;

			}

			num_pose_clustered++;

			if ( ( num_pose_clustered % 100 ) == 0 ){
				TR << "num_pose_clustered = " << num_pose_clustered << " num_cluster_centers = " << silent_struct_output_list_.size() << " time_taken_so_far = " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
			}
		}

	  TR << "After clustering, number of decoys: " << silent_struct_output_list_.size() << " from " << num_pose_clustered << " input poses " << std::endl;
		return;

	}

	/////////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_Clusterer::is_old_individual_suite_cluster( pose::Pose const & current_pose,
                                                     pose::Pose const & cluster_center_pose,
                                                     utility::vector1 < core::Size > const & calc_rms_res,
																									std::map< core::Size, core::Size > const & full_to_sub,
																									std::map< core::Size, bool > const & is_prepend_map,
																									core::Real const & cluster_radius ) const{


		utility::vector1< Real > rmsd_list( calc_rms_res.size(), 9999.99 );
		utility::vector1< bool > same_sugar_pucker_list( calc_rms_res.size(), false );


		for ( Size i = 1; i <= calc_rms_res.size(); i++ ){

 			Size const full_seq_num = calc_rms_res[i];

			if ( full_to_sub.find( full_seq_num ) == full_to_sub.end() ) utility_exit_with_message( "full_to_sub.find( full_seq_num ) == full_to_sub.end()!" );
			if ( is_prepend_map.find( full_seq_num ) == is_prepend_map.end() ) utility_exit_with_message( "is_prepend_map.find( full_seq_num ) == is_prepend_map.end()!" );

 			Size const seq_num = full_to_sub.find( full_seq_num )->second;
			bool is_prepend = is_prepend_map.find( full_seq_num )->second;

			//Important only if both pose are real
			same_sugar_pucker_list[i] = is_same_sugar_pucker( current_pose, cluster_center_pose, seq_num );

			bool const current_is_virtual_res = current_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" );
			bool const center_is_virtual_res = cluster_center_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" );

			bool const current_is_virtual_sugar = current_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RIBOSE" );
			bool const center_is_virtual_sugar = cluster_center_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RIBOSE" );


			if ( ignore_unmatched_virtual_res_ == false ){ //Sep 07. 2011

				if ( current_is_virtual_res != center_is_virtual_res ){
					return false; //current_pose is not part of this cluster center
				}

			}

			if ( current_is_virtual_res && center_is_virtual_res ){
				rmsd_list[i] = 8888.88;
				continue;
			}


			if ( PBP_clustering_at_chain_closure_ && working_parameters_->gap_size() == 0 ){ //new option Aug 15, 2010..include both phosphates in rmsd calculation at chain_break
				rmsd_list[i] =	 phosphate_base_phosphate_rmsd( current_pose, cluster_center_pose, seq_num, false /*ignore_virtual_atom*/ );
			} else{
				rmsd_list[i] =  suite_rmsd( current_pose, cluster_center_pose, seq_num, is_prepend, false /*ignore_virtaul_atom*/ );
			}

			if ( rmsd_list[i] > cluster_radius ) return false; //current_pose is not part of this cluster center


			if ( distinguish_pucker_ ){

				if ( current_is_virtual_sugar != center_is_virtual_sugar ){
					//New on Oct 09, 2011. This should NOT lead to any new changes, since virtual_sugar is usually accompanied by virtual_res at the neighoring nucleotide.
					return false;
				}

				bool check_pucker = true;

				if ( current_is_virtual_sugar && center_is_virtual_sugar ){
					//New on Oct 09, 2011. This should remove "false" new cluster where current pose and cluster center pose differ only by the pucker of a virtual_sugar.
					check_pucker = false;
				}

				if ( check_pucker && ( same_sugar_pucker_list[i] == false ) ){
					return false;
				}

			}


			/*
			if ( distinguish_pucker_ ){
				if ( rmsd_list[i] > cluster_radius || ( same_sugar_pucker_list[i] == false ) ) return false; //current_pose is not part of this cluster center
			} else{
				if ( rmsd_list[i] > cluster_radius ) return false; //current_pose is not part of this cluster center
			}
			*/

		}

		if ( verbose_ ){
			for ( Size i = 1; i <= calc_rms_res.size(); i++ ){

 				Size const full_seq_num = calc_rms_res[i];
 				Size const seq_num = full_to_sub.find( full_seq_num )->second;
				bool is_prepend = is_prepend_map.find( full_seq_num )->second;
				bool both_pose_res_is_virtual = false;
				if ( current_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) && cluster_center_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
					both_pose_res_is_virtual = true;
				}
				TR << "full_seq_num = " << full_seq_num << " seq_num = " << seq_num;
				output_boolean( " is_prepend = ", is_prepend, TR );
				output_boolean( " both_pose_res_is_virtual = ", both_pose_res_is_virtual, TR );
				TR << " same_pucker[" << i << "] = "; output_boolean( same_sugar_pucker_list[i], TR );
				print_sugar_pucker_state( " curr_pucker = ", get_residue_pucker_state( current_pose, seq_num ), TR );
				print_sugar_pucker_state( " center_pucker = ", get_residue_pucker_state( cluster_center_pose, seq_num ), TR );
				TR << " rmsd_list[" << i << "] = " << rmsd_list[i];

				if ( ignore_unmatched_virtual_res_ ){
					if ( current_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) != cluster_center_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
						TR << " Ignoring unmatched_virtual_res = ";
						output_boolean( " curr_virt = ", current_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ), TR );
						output_boolean( " center_virt = ", cluster_center_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ), TR );
					}
				}

				TR << std::endl;
			}

		}

		return true; //current_pose is not part of this cluster center
	}

	//////////////////////////////////////////////////////////////////
	core::pose::PoseOP
	StepWiseRNA_Clusterer::get_poseOP( Size const n ){

		using namespace core::pose;
		using namespace ObjexxFCL;

//		TR << "enter get_poseOP(" << n << ")" << std::endl;

		if ( keep_pose_in_memory_ ){
			if ( keep_pose_in_memory_hydrid_ == false || n <= max_memory_pose_num_ ){

				if ( pose_output_list_.size() < n ) utility_exit_with_message( "pose_output_list_.size() ( " + string_of( pose_output_list_.size() ) + " ) < n ( " + string_of( n ) + " )" );

				return pose_output_list_[ n ]; //if quick_alignment is true then this pose is already aligned to the quick_alignment_pose?
			}
		}

		//OK if reach this point means that pose is not stored in the pose_output_list_, need to extract it from the silent_file.
		core::pose::PoseOP pose_op( new Pose );
		silent_struct_output_list_[n]->fill_pose( *pose_op, *rsd_set_ );

		if ( optimize_memory_usage_ ) (*pose_op) = sliced_pose_job_params_.create_sliced_pose(*pose_op);

		if ( quick_alignment_ ) align_to_quick_alignment_pose( (*pose_op), tag_output_list_[n] );

		return pose_op;

	}


	//////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::setup_fail_triangle_inequailty_list( pose::Pose & current_pose, std::string const & tag, utility::vector1< bool > & fail_triangle_inequality_list ){

		using namespace core::scoring;
		using namespace ObjexxFCL;

		utility::vector1< core::Size > const & alignment_res =  get_act_alignment_res();
		utility::vector1 < core::Size > const & calc_rms_res = get_act_calc_rms_res();
		std::map< core::Size, core::Size > const & full_to_sub = get_act_full_to_sub();
		std::map< core::Size, bool > const & is_prepend_map = get_act_is_prepend_map();

		Size num_fail_triangle_inequality = 0;
		Size num_cluster_center_used = 0;

		Real current_score( 0.0 );

		getPoseExtraScores( current_pose, "score", current_score ); //Is this slow?

		fail_triangle_inequality_list.assign( silent_struct_output_list_.size(), false );

		for ( Size n = 1; n <= large_cluster_pose_list_.size(); n++ ){ //lowest score cluster center at the beginning of the list

			Real cluster_center_score( 0.0 );

			pose::Pose & cluster_center_pose = *( large_cluster_pose_list_[n] );

			if ( quick_alignment_ == false ) align_poses( current_pose, tag, cluster_center_pose, "large_cluster_center", alignment_res, align_only_over_base_atoms_ );

			getPoseExtraScores( cluster_center_pose, "score", cluster_center_score ); //Is this slow?

			if ( ( cluster_center_score + 0.001 ) > current_score ) break; //0.001 to make account for round off error.  Umm maybe faster without this break statement?

			// TR << "cluster_center_score= " << cluster_center_score << " current_score= " << current_score << std::endl;

			num_cluster_center_used++;

			Real const RMSD = rmsd_over_residue_list( current_pose, cluster_center_pose, calc_rms_res, full_to_sub, is_prepend_map, false );
			//problem is that bulge residues are excluded?? The weight of the RMSD and member.RMSD might not be the same... Aug 9, 2010

			utility::vector1< Cluster_Member > const & member_list = cluster_centers_neighbor_list_[n];

			for ( Size ii = 1; ii <= member_list.size(); ii++ ){ //lowest socre cluster center member at the beginning of the list.

				Cluster_Member const & member = member_list[ii];

				if ( ( member.score + 0.001 ) > current_score ) break; //0.001 to account for round off errors.

				if ( ( RMSD - member.RMSD ) > ( loop_cluster_radius_ + 0.02 ) ){ //satisfies triangle inequality

					if ( ( member.ID > all_pose_to_output_pose_ID_map_.size() ) || ( member.ID < 1 ) ){
						utility_exit_with_message( "member.ID ( " + string_of( member.ID )  + " ) > all_pose_to_output_pose_ID_map_.size() ( " +  string_of( all_pose_to_output_pose_ID_map_.size() ) + " ) " );
					}

					Size const output_pose_ID = all_pose_to_output_pose_ID_map_[member.ID];

					if ( output_pose_ID == 0 ) continue; //member is not a output_pose..

					if ( ( output_pose_ID > silent_struct_output_list_.size() ) || ( output_pose_ID < 1 ) ){
						utility_exit_with_message( "output_pose_ID ( " + string_of( output_pose_ID )  + " ) > silent_struct_output_list_.size() ( " +  string_of( silent_struct_output_list_.size() ) + " ) " );
					}

					if ( fail_triangle_inequality_list[output_pose_ID] == false ) num_fail_triangle_inequality++;

					fail_triangle_inequality_list[output_pose_ID] = true;

				}
			}
		}
		//TR << "num_cluster_center_used= " << num_cluster_center_used;
		TR << "num_fail_triangle_inequality = " << num_fail_triangle_inequality << " out_of = " << silent_struct_output_list_.size() << std::endl;

	}

	//////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_Clusterer::is_new_cluster_center_with_working_parameters( core::pose::PoseOP const & pose_op, std::string const & tag ){

		using namespace core::scoring;
		using namespace ObjexxFCL;

		//////////////////////////////////////////////////////////////////

		utility::vector1< core::Size > const & alignment_res =  get_act_alignment_res();
		utility::vector1 < core::Size > const & calc_rms_res = get_act_calc_rms_res();
		std::map< core::Size, core::Size > const & full_to_sub = get_act_full_to_sub();
		std::map< core::Size, bool > const & is_prepend_map = get_act_is_prepend_map();

		pose::Pose & current_pose = *( pose_op );

		if ( quick_alignment_ ) align_to_quick_alignment_pose( current_pose, tag );

		utility::vector1< bool > fail_triangle_inequality_list;

		if ( use_triangle_inequality_ ) setup_fail_triangle_inequailty_list( current_pose, tag, fail_triangle_inequality_list );

		for ( Size n = silent_struct_output_list_.size(); n >= 1; n-- ) {

			if ( use_triangle_inequality_ && ( fail_triangle_inequality_list[n] == true ) ) continue;

			pose::PoseOP const cluster_center_poseOP = get_poseOP( n );
			pose::Pose const & cluster_center_pose = ( *cluster_center_poseOP );
			std::string const & cluster_center_tag = tag_output_list_[n];


			if ( quick_alignment_ == false ) align_poses( current_pose, tag, cluster_center_pose, cluster_center_tag, alignment_res, align_only_over_base_atoms_ );

			//////////////////////////////////////////////////////////

			bool old_suite_cluster = is_old_individual_suite_cluster( current_pose, cluster_center_pose, calc_rms_res, full_to_sub, is_prepend_map, suite_cluster_radius_ );


			Real loop_rmsd = 99.99;

			if ( full_length_loop_rmsd_clustering_ ){
				if ( optimize_memory_usage_ ) utility_exit_with_message( "Both full_length_loop_rmsd_clustering_ and optimize_memory_usage_ equal true" );
				std::string const & full_sequence = working_parameters_->full_sequence();
				loop_rmsd = full_length_rmsd_over_residue_list( current_pose, cluster_center_pose, calc_rms_res, full_sequence, false /*verbose*/, false /*ignore_virtual_atom*/ );
			} else{
				loop_rmsd = rmsd_over_residue_list( current_pose, cluster_center_pose, calc_rms_res, full_to_sub, is_prepend_map, false /*verbose*/, false /*ignore_virtual_atom*/ );
			}

			bool old_loop_cluster = ( loop_rmsd < loop_cluster_radius_ );

			if ( verbose_ ){
				TR << "Between " << tag << " AND " << cluster_center_tag << ": loop_rmsd = " << loop_rmsd << " ";
				output_boolean( "is_old_suite_cluster = ", old_suite_cluster, TR ); TR << std::endl;
			}

			if ( old_suite_cluster == true && old_loop_cluster == true ){
				TR << tag << " is a neighbor of " << cluster_center_tag << std::endl;
				return false;
			}

		}

		return true; //new cluster center!
	}

	//////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_Clusterer::check_for_closeness_without_working_parameters( core::pose::PoseOP const & pose_op )
	{
		using namespace core::scoring;

		// go through the list backwards, because poses may be grouped by similarity --
		// the newest pose is probably closer to poses at the end of the list.
		for ( Size n = silent_struct_output_list_.size(); n >= 1; n-- ) {

			Real rmsd = all_atom_rmsd( *( get_poseOP( n ) ), *pose_op );

			if ( rmsd < whole_struct_cluster_radius_ )	return false;
		}
		return true;
	}

	//////////////////////////////////////////////////////////////////

	bool
	StepWiseRNA_Clusterer::check_for_closeness( core::pose::PoseOP const & pose_op, std::string const & tag ){

		if ( skip_clustering_ == true ){
			utility_exit_with_message( "skip_clustering == true but StepWiseRNA_Clusterer::check_for_closeness() is called! " );
		}

		if ( working_parameters_exist_ ){
			return is_new_cluster_center_with_working_parameters( pose_op, tag );
		} else{
			return check_for_closeness_without_working_parameters( pose_op );
		}

	}



	/////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::output_silent_file( std::string const & silent_file ){

		using namespace core::io::silent;

		SilentFileData silent_file_data;

		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {

			SilentStructOP & s( silent_struct_output_list_[ n ] );

			if ( rename_tags_ ){
				s->add_comment( "PARENT_TAG", s->decoy_tag() );

				std::string tag;
				if ( add_lead_zero_to_tag_ ){
					tag = "S_" + ObjexxFCL::lead_zero_string_of( n - 1 /* start with zero */, 6 );
				} else{
					tag = "S_" + ObjexxFCL::string_of( n - 1 /* start with zero */ );
				}

				s->set_decoy_tag( tag );
			}

			silent_file_data.write_silent_struct( *s, silent_file, false /*write score only*/ );

		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::recalculate_rmsd_and_output_silent_file( std::string const & silent_file,
				                                                     protocols::stepwise::legacy::sampling::rna::StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup,
																													bool const write_score_only ){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		clock_t const time_start( clock() );

		output_title_text( "ENTER StepWiseRNA_Clusterer::recalculate_rmsd_and_output_silent_file()", TR );

		if ( working_parameters_exist_ == false ) utility_exit_with_message( "working_parameters_exist_ == false!" );

		///This could actually work...but it is just yet tested!
		if ( working_parameters_->is_simple_full_length_job_params() == true ) utility_exit_with_message( "working_parameters_->is_simple_full_length_job_params() == true!" );

		output_boolean( "write_score_only = ", write_score_only, TR ); TR << std::endl;

		utility::vector1< core::Size > const & working_best_alignment = working_parameters_->working_best_alignment();
		utility::vector1< core::Size > const & working_native_alignment = working_parameters_->working_native_alignment();
		std::string const & full_sequence = working_parameters_->full_sequence();


		stepwise_rna_pose_setup->set_verbose( true ); //New OPTION, Mar 22

		if ( tag_output_list_.size() != silent_struct_output_list_.size() ) utility_exit_with_message( "pose_output_list_.size() != silent_struct_output_list_!" );

		SilentFileData silent_file_data;

		utility::vector1 < core::Size > const & calc_rms_res = working_parameters_->calc_rms_res();
		std::map< core::Size, core::Size > const & full_to_sub = working_parameters_->const_full_to_sub();

		//bool const ignore_min_decoys=true; //Over the keep min_decoy mode...Comment out on Dec 11, 2011.

		//float best_score=9999999999999; //lead to server-test build error at R47198; Feb 02, 2012
		//Real best_score=999999999; //Should Fix server-test build error BUT yet not tested; Feb 02, 2012

		std::map< core::Size, bool > is_prepend_map;
		is_prepend_map.clear();

		is_prepend_map = working_parameters_->is_prepend_map();

		bool is_full_length_pose = true; //Will be init first time the loop is tranversed.

		bool is_valid_first_struct = true;

		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {

			std::string tag = tag_output_list_[n];
			SilentStructOP s( silent_struct_output_list_[ n ] );

			if ( ( n % 100 ) == 0 ){
				TR << "recalculate rmsd for " << tag << " n = " << n << " taken time " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC  << std::endl;
			}
			core::pose::PoseOP const pose_op = get_poseOP( n );
			core::pose::Pose pose = (*pose_op);

			if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( pose, tag ) == true ) continue;

			Real score( 0.0 );
			getPoseExtraScores( pose, "score", score );

			//This kinda weird in that setup_native_pose actually set the working_native_pose in the working_parameters...this interdependency is not good!
			if ( is_valid_first_struct ){
				//best_score = score;
				stepwise_rna_pose_setup->setup_native_pose( pose ); //Setup native_pose;
				is_full_length_pose = ( pose.total_residue() == full_sequence.size() ) ? true : false;
				output_boolean( "is_full_length_pose = ", is_full_length_pose, TR ); TR << std::endl;

				is_valid_first_struct = false;
			}

			//if(score > best_score + score_diff_cut_) break; //Comment out on Dec 11, 2011.

			PoseOP native = new Pose;
			( *native ) = ( *working_parameters_->working_native_pose() ); //Hard copy...

			align_poses( ( *native ), "native", pose, tag, working_best_alignment, align_only_over_base_atoms_ );


			s->add_energy( "NEW_all_rms", rms_at_corresponding_heavy_atoms( pose, *native ) );
			s->add_energy( "NEW_loop_rmsd", rmsd_over_residue_list( pose, *native, calc_rms_res, full_to_sub, is_prepend_map, false, false ) );

			///////////////////////////////////////////////////////////////////////////////////////////////

			if ( working_native_alignment.size() != 0 ){ //user specify which residue to align with native.
				align_poses( ( *native ), "native", pose, tag, working_native_alignment, align_only_over_base_atoms_ );
			} else{ //default
				align_poses( ( *native ), "native", pose, tag, working_best_alignment, align_only_over_base_atoms_ ); //REDUNDANT
			}
			s->add_energy( "NEW_O_loop_rmsd", rmsd_over_residue_list( pose, *native, calc_rms_res, full_to_sub, is_prepend_map, false, false ) );

			if ( is_full_length_pose ){
				s->add_energy( "NEW_Full_L_rmsd", full_length_rmsd_over_residue_list( pose, *native, calc_rms_res, full_sequence, false, false ) );
			}

			////////Simple loop RMSD exclude only virtual atoms in native_pdb (mostly just the native virtual_res)//////////////
			core::pose::Pose curr_pose_no_variants = pose;
			remove_all_variant_types( curr_pose_no_variants ); //This remove all virtual_atoms!

			if ( working_native_alignment.size() != 0 ){ //user specify which residue to align with native.
				align_poses( ( *native ), "native", curr_pose_no_variants, tag + "_no_variants", working_native_alignment, align_only_over_base_atoms_ );
			} else{ //default
				align_poses( ( *native ), "native", curr_pose_no_variants, tag + "_no_variants", working_best_alignment, align_only_over_base_atoms_ );
			}

			s->add_energy( "NEW_NAT_rmsd", rmsd_over_residue_list( curr_pose_no_variants, *native, calc_rms_res, full_to_sub, is_prepend_map, false /*verbose*/, true /*ignore_virtual_atom*/ ) );

			////March 7, 2011....Output BASE-PAIRS STATISTIC///////////////////////////////
			//utility::vector1< core::Size > const working_calc_rms_res=apply_full_to_sub_mapping(calc_rms_res, working_parameters_);

			//Nov 01, 2011 WARNING THIS currently does not work if there is protonated Adenosine!
			//add_base_pair_stats( s, pose, *native, working_calc_rms_res);


			if ( rename_tags_ ){
				s->add_comment( "PARENT_TAG", s->decoy_tag() );
				if ( add_lead_zero_to_tag_ ) tag = "S_" + ObjexxFCL::lead_zero_string_of( n - 1 /* start with zero */, 6 );
				s->set_decoy_tag( tag );
			}

			///////////////////////////////////////////////////////////////////////////////////////////////

			silent_file_data.write_silent_struct( *s, silent_file, write_score_only );

		}

		TR << "Total # pose alignment and rmsd recalculation= " << silent_struct_output_list_.size() << std::endl;
		TR << "Total recalculate rmsd time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		output_title_text( "EXIT StepWiseRNA_Clusterer::recalculate_rmsd_and_output_silent_file()", TR );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::get_best_neighboring_shift_RMSD_and_output_silent_file( std::string const & silent_file ){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		clock_t const time_start( clock() );

		output_title_text( "ENTER StepWiseRNA_Clusterer::get_best_neighboring_shift_RMSD_and_output_silent_file()", TR );

		if ( working_parameters_exist_ == false ) utility_exit_with_message( "working_parameters_exist_ == false!" );

		///This could actually work...but it is just yet tested!
		if ( working_parameters_->is_simple_full_length_job_params() == true ) utility_exit_with_message( "working_parameters_->is_simple_full_length_job_params() == true!" );

		utility::vector1< core::Size > const & working_best_alignment = working_parameters_->working_best_alignment();

		utility::vector1 < core::Size > const & calc_rms_res = working_parameters_->calc_rms_res();
		std::map< core::Size, core::Size > const & full_to_sub = working_parameters_->const_full_to_sub();

		std::map< core::Size, bool > is_prepend_map;
		is_prepend_map.clear();

		is_prepend_map = working_parameters_->is_prepend_map();

		std::string const & full_sequence = working_parameters_->full_sequence();

		SilentFileData silent_file_data;

		TR << "loop_cluster_radius_ = "  << loop_cluster_radius_  << std::endl;
		TR << "suite_cluster_radius_ = "  << suite_cluster_radius_ << std::endl;

		bool is_full_length_pose = true; //Will be init first time the loop is tranversed.

		bool is_valid_first_struct = true;

		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {

			std::string tag = tag_output_list_[n];
			SilentStructOP s( silent_struct_output_list_[ n ] );

			if ( ( n % 100 ) == 0 ){
				TR << "find_best_neighboring_shift_rmsd for " << tag << " n = " << n << " taken time " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC  << std::endl;
			}

			core::pose::PoseOP const current_pose_op = get_poseOP( n );
			core::pose::Pose current_pose = ( *current_pose_op );

			if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( current_pose, tag ) == true ) continue;

			if ( is_valid_first_struct ){
				 is_full_length_pose = ( current_pose.total_residue() == full_sequence.size() ) ? true : false;
				 output_boolean( "is_full_length_pose = ", is_full_length_pose, TR ); TR << std::endl;

				is_valid_first_struct = false;

			}

			Real start_score( 0.0 );
			bool has_total_score = getPoseExtraScores( current_pose, "score", start_score );
			if ( has_total_score == false ) utility_exit_with_message( "current_pose ( " + tag + " ) missing total score!" );


			Real start_shift_score( 0.0 );
			bool has_shift_score = getPoseExtraScores( current_pose, "shift_score", start_shift_score );
			if ( has_shift_score == false ) utility_exit_with_message( "current_pose ( " + tag + " ) missing shift_score!" );

			if ( quick_alignment_ ) align_to_quick_alignment_pose( current_pose, tag );

			Real best_shift_score = start_shift_score;
			std::string best_shift_tag = tag;

			for ( Size other_pose_ID = 1 ; other_pose_ID  <= silent_struct_output_list_.size(); other_pose_ID ++ ) {

				std::string other_tag = tag_output_list_[ other_pose_ID ];

				core::pose::PoseOP const other_pose_op = get_poseOP( other_pose_ID );
				core::pose::Pose other_pose = ( *other_pose_op );

				if ( protocols::stepwise::sampling::rna::check_for_messed_up_structure( other_pose, other_tag ) == true ) continue;

				if ( quick_alignment_ == false ) align_poses( other_pose, other_tag, current_pose, tag, working_best_alignment, align_only_over_base_atoms_ );


				bool old_suite_cluster = is_old_individual_suite_cluster( current_pose, other_pose, calc_rms_res, full_to_sub, is_prepend_map, suite_cluster_radius_ );

				Real loop_rmsd = 99.99;

				if ( is_full_length_pose ){
					if ( optimize_memory_usage_ ) utility_exit_with_message( "Both full_length_loop_rmsd_clustering_ and optimize_memory_usage_ equal true" );
					loop_rmsd = full_length_rmsd_over_residue_list( current_pose, other_pose, calc_rms_res, full_sequence, false /*verbose*/, false /*ignore_virtual_atom*/ );
				} else{
					loop_rmsd = rmsd_over_residue_list( current_pose, other_pose, calc_rms_res, full_to_sub, is_prepend_map, false /*verbose*/, false /*ignore_virtual_atom*/ );
				}

				bool old_loop_cluster = ( loop_rmsd < loop_cluster_radius_ );

				if ( old_suite_cluster == true && old_loop_cluster == true ){
					if ( verbose_ ) TR << tag << " is a neighbor of " << other_tag << std::endl;
				} else{
					continue; //Not a neighor!
				}

				Real other_shift_score( 0.0 );
				bool has_shift_score = getPoseExtraScores( other_pose, "shift_score", other_shift_score );
				if ( has_shift_score == false ) utility_exit_with_message( "other_pose ( " + other_tag + " ) missing shift_score" );

				if ( other_shift_score < best_shift_score ){
					best_shift_score = other_shift_score;
					best_shift_tag = other_tag;
				}
			}

			float const new_score = start_score - start_shift_score + best_shift_score;

			//setPoseExtraScores(pose, "score", new_score);
			//setPoseExtraScores(pose, "shift_score", best_shift_score);
			//setPoseExtraScores(pose, "self_shift_score", start_shift_score);
			//add_score_line_string(pose, "src_shift_tag", best_shift_tag);

			SilentStructOP new_silent_struct = s->clone(); //Important to create new one since shift_score is being changed!!!

			new_silent_struct->add_energy( "score", new_score );
			new_silent_struct->add_energy( "shift_score", best_shift_score );
			new_silent_struct->add_energy( "self_shift_score", start_shift_score );
			new_silent_struct->add_string_value( "src_shift_tag", best_shift_tag );

			silent_file_data.write_silent_struct( *new_silent_struct, silent_file, false /*write_score_only*/ );

		}

		TR << "silent_struct_output_list_.size() = " << silent_struct_output_list_.size() << std::endl;
		TR << "Total get_best_neighboring_shift_RMSD_and_output_silent_file time : " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		output_title_text( "EXIT StepWiseRNA_Clusterer::get_best_neighboring_shift_RMSD_and_output_silent_file()", TR );


	}


	////////////////////////////////////Sept 06, 2011 (For post_processing)////////////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::create_tags_map(){


		output_title_text( "ENTER StepWiseRNA_Clusterer::create_tag_map()", TR );

		input_->reset(); //reset the silentfile stream to the beginning..

		current_tags_map_.clear();
		parent_tags_map_.clear();

		//mymap.count(c)

		while ( input_->has_another_pose() ) {

			core::io::silent::SilentStructOP const silent_struct( input_->next_struct() );

			std::string const tag = silent_struct->decoy_tag();

			if ( ignore_FARFAR_no_auto_bulge_tag_ ){

				if ( current_tags_map_.count( tag ) != 0 ) utility_exit_with_message( tag  + " already exist in current_tags_map_!" );

				current_tags_map_[tag] = true;

			}


			if ( ignore_FARFAR_no_auto_bulge_parent_tag_ ){

				if ( silent_struct->has_parent_remark( "PARENT_TAG" ) == false ){
					TR << "silent_struct ( " << tag  << " ) missing PARENT_TAG!" << std::endl;
					silent_struct->print_parent_remarks( TR );
					utility_exit_with_message( "silent_struct ( " + tag + " ) missing PARENT_TAG!" );
				}

				std::string const parent_tag = silent_struct->get_parent_remark( "PARENT_TAG" );

				if ( parent_tags_map_.count( parent_tag ) != 0 ) utility_exit_with_message( parent_tag  + " already exist in parent_tags_map_!" );

				parent_tags_map_[parent_tag] = true;

			}

		}

		input_->reset(); //reset the silentfile stream to the beginning..

		output_title_text( "EXIT StepWiseRNA_Clusterer::create_tag_map()", TR );

	}

	////////////////////////////////////Sept 06, 2011 (For post_processing)////////////////////////////////////////////////////
	bool
	StepWiseRNA_Clusterer::pass_FARFAR_no_auto_bulge_filter( core::io::silent::SilentStructOP const & silent_struct ) const{

		//This only effects to FARFAR models!
		//For the purpose of clustering, assume that the NO_AUTO_BULGE belong to the same same cluster as the WITH_AUTO_BULGE pose.
		//So if a instance of the pose with WITH_AUTO_BULGE exist in the silent_file then ignore the NO_AUTO_BULGE version of the pose!
		//The WITH_AUTO_BULGE pose always have better energy!

		std::string const NO_AUTO_BULGE_STR = "_NO_AUTO_BULGE";
		std::string const WITH_AUTO_BULGE_STR = "_WITH_AUTO_BULGE";

		std::string const tag = silent_struct->decoy_tag();

		if ( ignore_FARFAR_no_auto_bulge_tag_ ){

			if ( current_tags_map_.size() == 0 ) utility_exit_with_message( "current_tags_map_ is empty!" );

			size_t found_curr_tag;
			found_curr_tag = tag.find( NO_AUTO_BULGE_STR );

			if ( found_curr_tag != std::string::npos ){

				std::string WITH_AUTO_BULGE_curr_tag = tag;

				WITH_AUTO_BULGE_curr_tag.replace( found_curr_tag, 14, WITH_AUTO_BULGE_STR );

				if ( current_tags_map_.count( WITH_AUTO_BULGE_curr_tag ) > 0 ){

					TR << "Ignoring NO_AUTO_BULGE pose: " << tag << " since WITH_BULGE_curr_tag: " << WITH_AUTO_BULGE_curr_tag << " exist!" << std::endl;
					return false;

				}

			}

			//if(tag.substr(tag.size()-14, 14)==NO_AUTO_BULGE_STR){ //Problem with this is that it doesn't account for the renamed tag S_0 to S_0_1 possibility!
			//if(found==std::string::npos) utility_exit_with_message("CANNOT FIND NO_AUTO_BULGE_STR in current_tag: " + tag);

		}


		if ( ignore_FARFAR_no_auto_bulge_parent_tag_ ){

			if ( parent_tags_map_.size() == 0 ) utility_exit_with_message( "parent_tags_map_ is empty!" );

			if ( silent_struct->has_parent_remark( "PARENT_TAG" ) == false ) utility_exit_with_message( "silent_struct ( " + tag + " missing PARENT_TAG!" );

			std::string const parent_tag = silent_struct->get_parent_remark( "PARENT_TAG" );

			size_t found_parent_tag;
			found_parent_tag = parent_tag.find( NO_AUTO_BULGE_STR );

			if ( found_parent_tag != std::string::npos ){

				std::string WITH_AUTO_BULGE_parent_tag = parent_tag;

				WITH_AUTO_BULGE_parent_tag.replace( found_parent_tag, 14, WITH_AUTO_BULGE_STR );

				if ( parent_tags_map_.count( WITH_AUTO_BULGE_parent_tag ) > 0 ){

					TR << "Ignoring NO_AUTO_BULGE pose: " << tag << " with parent_tag: " << parent_tag;
					TR << ", since WITH_AUTO_BULGE_parent_tag " << WITH_AUTO_BULGE_parent_tag << " exist!" <<  std::endl;
					return false;

				}
			}
		}

		return true;

	}


	//////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::set_working_parameters( protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP & working_parameters ){

		working_parameters_ = working_parameters;

	}
	//////////////////////////////////////////////
	void
	StepWiseRNA_Clusterer::set_working_parameters_exist( bool const working_parameters_exist ){

		working_parameters_exist_ = working_parameters_exist;

	}

	//////////////////////////////////////////////
	utility::vector1< core::Size > const &
	StepWiseRNA_Clusterer::get_act_alignment_res()  const {
		utility::vector1< core::Size > const & alignment_res =   ( optimize_memory_usage_ ) ? sliced_pose_job_params_.sliced_pose_best_alignment: working_parameters_->working_best_alignment() ;
		return alignment_res;
	}

	utility::vector1 < core::Size > const &
	StepWiseRNA_Clusterer::get_act_calc_rms_res()	 const {
		utility::vector1 < core::Size > const & calc_rms_res = ( optimize_memory_usage_ ) ? sliced_pose_job_params_.sliced_pose_calc_rms_res : working_parameters_->calc_rms_res();
		return calc_rms_res;
	}

	std::map< core::Size, core::Size > const &
	StepWiseRNA_Clusterer::get_act_full_to_sub()	 const {
		std::map< core::Size, core::Size > const & full_to_sub = ( optimize_memory_usage_ ) ? sliced_pose_job_params_.sliced_pose_full_to_sub   : working_parameters_->const_full_to_sub();
		return full_to_sub;
	}

	std::map< core::Size, bool > const &
	StepWiseRNA_Clusterer::get_act_is_prepend_map() const 	{
		std::map< core::Size, bool > const & is_prepend_map =   ( optimize_memory_usage_ ) ? sliced_pose_job_params_.sliced_pose_is_prepend_map: working_parameters_->is_prepend_map();
		return is_prepend_map ;
	}



	//////////////////////////////////////////////
	void
	SlicedPoseWorkingParameters::setup( protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP & working_parameters ){

		output_title_text( "Enter SlicedPoseWorkingParameters::setup()", TR );

		is_setup_ = true;

		Size const nres = ( working_parameters->working_sequence() ).size();
		utility::vector1< core::Size > const & working_best_alignment( working_parameters->working_best_alignment() );
		utility::vector1 < core::Size > const & calc_rms_res = working_parameters->calc_rms_res();
		std::map< core::Size, bool > const & is_prepend_map = working_parameters->is_prepend_map();
		std::map< core::Size, core::Size > const & sub_to_full( working_parameters->const_sub_to_full() );


		utility::vector1< core::Size > working_calc_rms_res = apply_full_to_sub_mapping( calc_rms_res, working_parameters );

		Size sliced_seq_num = 1;
		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
			bool keep_res = false;

			if ( working_best_alignment.has_value( seq_num ) ) {
				TR << "seq_num " << seq_num << " is in working_best_alignment res " << std::endl;
				keep_res = true;
			}


			if ( working_calc_rms_res.has_value( seq_num ) ){
				TR << "seq_num " << seq_num << " is in working_calc_rms_res " << std::endl;
				keep_res = true;
			}

			if ( keep_res == false && ( seq_num + 1 ) <= nres && working_calc_rms_res.has_value( seq_num + 1 ) ){
				TR << "seq_num " << seq_num << " is in working_calc_rms_res - 1 " << std::endl;
				keep_res = true;
			}

			if ( keep_res == false && ( seq_num - 1 ) >= 1 && working_calc_rms_res.has_value( seq_num - 1 ) ){
				TR << "seq_num " << seq_num << " is in working_calc_rms_res + 1 " << std::endl;
				keep_res = true;
			}


			is_sliced_res_.push_back( keep_res );

			if ( keep_res == true ){
				working_to_sliced_res_map_.push_back( sliced_seq_num );
				sliced_to_working_res_map_.push_back( seq_num );
				sliced_seq_num++;
			} else{
				working_to_sliced_res_map_.push_back( 0 );
			}

		}

		TR << "------------Before slice to After slice seq_num------------" << std::endl;
		for ( Size seq_num = 1; seq_num <= working_to_sliced_res_map_.size(); seq_num++ ){
			TR << seq_num << "----> " << working_to_sliced_res_map_[seq_num] << std::endl;

			if ( working_best_alignment.has_value( seq_num ) ) sliced_pose_best_alignment.push_back( working_to_sliced_res_map_[seq_num] ) ;
			if ( working_calc_rms_res.has_value( seq_num ) ) sliced_pose_calc_rms_res.push_back( working_to_sliced_res_map_[seq_num] ) ;
		}
		TR << "-----------------------------------------------------------" << std::endl;

		TR << "------------After slice to Before slice seq_num------------" << std::endl;
		for ( Size seq_num = 1; seq_num <= sliced_to_working_res_map_.size(); seq_num++ ){
			TR << seq_num << "----> " << sliced_to_working_res_map_[seq_num] << std::endl;
			sliced_pose_full_to_sub[seq_num] = seq_num; //identity

			Size const working_seq_num = sliced_to_working_res_map_[seq_num];
			Size const full_seq_num = sub_to_full.find( working_seq_num )->second;
			bool const is_prepend = is_prepend_map.find( full_seq_num )->second;
			sliced_pose_is_prepend_map[seq_num] = ( is_prepend );

		}
		TR << "-----------------------------------------------------------" << std::endl;

		//////////////////////////////////////////////
		bool in_delete_range = false;

		Size range_end = 0;
		Size range_begin = 0;

		for ( Size seq_num = 1; seq_num <= is_sliced_res_.size() + 1; seq_num++ ){ //optimization for using delete_residue_range_slow instead of delete_residue_slow


			if ( in_delete_range == false ){

				if ( seq_num == ( is_sliced_res_.size() + 1 ) ) continue;

				if ( is_sliced_res_[seq_num] == false ){
					range_begin = seq_num;
					in_delete_range = true;
				}

			} else{
				if ( seq_num == ( is_sliced_res_.size() + 1 ) || is_sliced_res_[seq_num] == true ){
					range_end = seq_num - 1;	//This obviously fail if seq_num=0...but this cannot occur since in_delete_range is false at first cycle.
					in_delete_range = false;

					delete_res_range_list_.push_back( std::make_pair( range_begin, range_end ) );
					range_end = 0;
					range_begin = 0;
				}
			}
		}
		//////////////////////////////////////////////


		//output debug
		output_seq_num_list( "sliced_pose_best_alignment = ", sliced_pose_best_alignment, TR, 50 );
		output_seq_num_list( "sliced_pose_calc_rms_res = ", sliced_pose_calc_rms_res, TR, 50 );
		output_is_prepend_map( "sliced_pose_is_prepend_map = ", sliced_pose_is_prepend_map, working_to_sliced_res_map_.size(), TR, 50 );
		output_pair_size( delete_res_range_list_, "delete_res_range_list = ", TR, 50 );

		output_title_text( "Exit SlicedPoseWorkingParameters::setup()", TR );

	}

	//////////////////////////////////////////////
	core::pose::Pose
	SlicedPoseWorkingParameters::create_sliced_pose( core::pose::Pose const & working_pose ){

		using namespace core::conformation;
		using namespace core::pose;
		using namespace ObjexxFCL;

		if ( is_setup_ == false ){
			utility_exit_with_message( "is_setup_ == false" );
		}

		core::pose::Pose sliced_pose = working_pose;

		if ( is_sliced_res_.size() != working_pose.total_residue() ){
			utility_exit_with_message( "is_sliced_res.size() ( " + string_of( is_sliced_res_.size() ) + " ) != working_pose.total_residue() ( " + string_of( working_pose.total_residue() ) + " )" );
		}

//		for(Size seq_num=is_sliced_res_.size(); seq_num>=1; seq_num--){
//			if(is_sliced_res_[seq_num]==false){
//				sliced_pose.conformation().delete_residue_slow(seq_num);
//				sliced_pose.conformation().delete_polymer_residue(seq_num); //doesn't work at jump_point...
//			}
//		}

		for ( Size n = delete_res_range_list_.size(); n >= 1; n-- ){
			sliced_pose.conformation().delete_residue_range_slow( delete_res_range_list_[n].first, delete_res_range_list_[n].second );
		}

		if ( sliced_pose.total_residue() != sliced_to_working_res_map_.size() ){
			utility_exit_with_message( "working_pose.total_res() ( " + string_of( working_pose.total_residue() ) + " ) != sliced_to_working.size() ( " + string_of( sliced_to_working_res_map_.size() ) + " )" );
		}


//		working_pose.dump_pdb( "clusterer_working_pose.pdb");
//		sliced_pose.dump_pdb( "clusterer_sliced_pose.pdb");

//		exit(1);

		return sliced_pose;

	}

	//////////////////////////////////////////////

} //rna
} //sampling
} //stepwise
} //protocols
