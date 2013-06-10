// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_AnalyticalLoopCloseSampler
/// @brief Alternative SWA sampling using Analytical Loop Closure
/// @detailed
/// @author Fang-Chieh Chou

//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_AnalyticalLoopCloseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBase_Sampler_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_AnalyticalLoopCloseSampler.hh>
#include <protocols/swa/rna/RNA_LoopCloseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.hh>
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_Bin_Screener.hh>

#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <core/scoring/rna/RNA_Util.hh>
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
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
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
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/rna/RNA_IdealCoord.hh>
#include <core/io/pdb/pose_io.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>
#include <numeric/random/random.hh>
#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>


using namespace core;
using core::Real;
using io::pdb::dump_pdb;
static numeric::random::RandomGenerator RG(19912388);  // <- Magic number, do not change it!

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of RNA
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_AnalyticalLoopCloseSampler::StepWiseRNA_AnalyticalLoopCloseSampler ( StepWiseRNA_JobParametersCOP & job_parameters ) :
	job_parameters_ ( job_parameters ),
	sfd_ ( new core::io::silent::SilentFileData ),
	scorefxn_ ( core::scoring::ScoreFunctionFactory::create_score_function ( "rna_hires.wts" ) ), // can be replaced from the outside
	silent_file_ ( "silent_file.txt" ),
	output_filename_ ( "data.txt" ),
	bin_size_ ( 20 ),
	rep_cutoff_ ( 4.0 ),
	num_pose_kept_ ( 108 ),
	multiplier_ ( 2 ), //Sort and cluster poses when the number of pose is pose_data_list exceed multiplier*num_pose_kept,
	cluster_rmsd_ ( 0.5001 ),
	verbose_ ( false ),
	native_rmsd_screen_ ( false ),
	native_screen_rmsd_cutoff_ ( 2.0 ),
	o2star_screen_ ( true ),
	fast_ ( false ),
	medium_fast_ ( false ),
	centroid_screen_ ( true ),
	VDW_atr_rep_screen_ ( true ),
	distinguish_pucker_ ( true ),
	current_score_cutoff_ ( 99999 ), //New option May 12, 2010
	finer_sampling_at_chain_closure_ ( false ), //New option Jun 10 2010
	PBP_clustering_at_chain_closure_ ( false ), //New option Aug 15 2010
	extra_anti_chi_rotamer_(false), //Split to syn and anti on June 16, 2011
	extra_syn_chi_rotamer_(false), //Split to syn and anti on June 16, 2011
	use_phenix_geo_(false), //The standard geo from PHENIX
	choose_random_( false ),
	force_centroid_interaction_( false )
{
	set_native_pose ( job_parameters_->working_native_pose() );
	////////////////Parin Feb 28, 2010////////////////////////////////////////////////
	utility::vector1 < core::Size > const & rmsd_res_list = job_parameters_->rmsd_res_list();
	working_rmsd_res_ = apply_full_to_sub_mapping ( rmsd_res_list, job_parameters );
	std::map< core::Size, bool > const & Is_prepend_map = job_parameters_->Is_prepend_map();
	Output_is_prepend_map ( "Is_prepend_map= ", Is_prepend_map , job_parameters_->full_sequence().size(), 30 );
	Output_seq_num_list ( "rmsd_res= ", rmsd_res_list, 30 );
	Output_seq_num_list ( "working_rmsd_res= ", working_rmsd_res_, 30 );
	////////////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseRNA_AnalyticalLoopCloseSampler::~StepWiseRNA_AnalyticalLoopCloseSampler()
{}

//////////////////////////////////////////////////////////////////////////
bool
sort_criteria2 ( pose_data_struct2  pose_data_1, pose_data_struct2 pose_data_2 ) { //This function used to be call sort_criteria2
	return ( pose_data_1.score < pose_data_2.score );
}
////////////////////////////////////////////////////////////////////////////


void
StepWiseRNA_AnalyticalLoopCloseSampler::apply ( core::pose::Pose & pose ) {
	using namespace ObjexxFCL;
	Output_title_text ( "Enter StepWiseRNA_AnalyticalLoopCloseSampler::apply" );
	clock_t const time_start ( clock() );
	//output screen options
	std::cout << "--------SCREEN OPTIONS---------- " << std::endl;
	Output_boolean ( "native_rmsd_screen = ", native_rmsd_screen_ );
	std::cout << std::endl;
	std::cout << "native_screen_rmsd_cutoff = " << native_screen_rmsd_cutoff_ << std::endl;
	Output_boolean ( "o2star_screen = ", o2star_screen_ );
	std::cout << std::endl;
	Output_seq_num_list ( "working_moving_partition_pos= ", job_parameters_->working_moving_partition_pos(), 40 );
	Output_boolean ( "centroid_screen = ", centroid_screen_ );
	std::cout << std::endl;
	Output_boolean ( "VDW_atr_rep_screen = ", VDW_atr_rep_screen_ );
	std::cout << std::endl;
	std::cout << "--------------------------------" << std::endl;

	Pose const pose_save = pose;
	pose = pose_save; //this recopy is useful for triggering graphics.

	initialize_scorefunctions(); //////////////// Sets up scorefunctions for a bunch of different screens /////////

	utility::vector1< pose_data_struct2 > pose_data_list;
	if ( pose_data_list_.size() > 0 ) {
		std::cout << "Previous poses exist in sampler: " << pose_data_list_.size() << std::endl;
		pose_data_list = pose_data_list_;  // give the class a little memory ... later need to include reset() function.
	}

	standard_sampling ( pose, pose_data_list );

	std::cout << "Total time in StepWiseRNA_AnalyticalLoopCloseSampler::apply " << static_cast<Real> ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
	Output_title_text ( "Exit StepWiseRNA_AnalyticalLoopCloseSampler::apply" );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
StepWiseRNA_AnalyticalLoopCloseSampler::standard_sampling ( core::pose::Pose & pose , utility::vector1< pose_data_struct2 > & pose_data_list ) {
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace protocols::rna;
	using namespace core::id;
	Output_title_text ( "Enter StepWiseRNA_AnalyticalLoopCloseSampler::standard_sampling" );
	clock_t const time_start ( clock() );
	SilentFileData silent_file_data;
	Size const moving_res ( job_parameters_->working_moving_res() ); // Might not corresponds to user input.
	Size const moving_suite ( job_parameters_->working_moving_suite() ); // dofs betweeen this value and value+1 actually move.
	bool const Is_prepend ( job_parameters_->Is_prepend() );
	bool const Is_internal ( job_parameters_->Is_internal() ); // no cutpoints before or after moving_res.
	Size const actually_moving_res ( job_parameters_->actually_moving_res() ); //Now same as moving_res
	Size const gap_size ( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */
	utility::vector1 < core::Size > const & cutpoint_closed_list = job_parameters_->cutpoint_closed_list();
	Size const num_nucleotides ( job_parameters_->working_moving_res_list().size() );
	Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
	Size const cutpoint_closed = cutpoint_closed_list[1];
	std::cout << " NUM_NUCLEOTIDES= " <<  num_nucleotides << std::endl;
	std::cout << " GAP SIZE " << gap_size << std::endl;
	std::cout << " MOVING RES " << moving_res << std::endl;
	std::cout << " MOVING SUITE " << moving_suite << std::endl;
	Output_boolean ( " PREPEND ", Is_prepend );
	std::cout << std::endl;
	Output_boolean ( " INTERNAL ", Is_internal );
	std::cout << std::endl;
	Output_boolean ( " distinguish_pucker_ ", distinguish_pucker_ );
	std::cout << std::endl;
	/////////////////////////////// O2star sampling/virtualization //////////////////////////
	Pose pose_with_virtual_O2star_hydrogen = pose;
	Add_virtual_O2Star_hydrogen ( pose_with_virtual_O2star_hydrogen );
//		Mod out on Apr 3, 2010 Important so that pose can remember the o2star torsion from previous minimization step. Check if this cause problems!!
//    This is true only in the INTERNAL CASE!...May 31, 2010
	// if o2star_screen, pose 2'-OH will be sampled!
	pose::Pose pose_with_original_HO2star_torsion;

	if ( o2star_screen_ ) {
		pose_with_original_HO2star_torsion = pose;
		initialize_o2star_packer_task ( pose );
	} else {
		// Otherwise, virtualize the 2-OH.
		pose = pose_with_virtual_O2star_hydrogen;
	}

	Real base_rep_score ( -99999 ), base_atr_score ( -99999 );
	get_base_atr_rep_score ( pose_with_virtual_O2star_hydrogen, base_atr_score, base_rep_score );
	//////////////////////////////////////////Setup Atr_rep_screening/////////////////////////////////////////////////
	pose::Pose screening_pose = pose_with_virtual_O2star_hydrogen; //Hard copy
	//Necessary for the case where gap_size == 0. In this case, Pose_setup does not automatically create a VIRTUAL_PHOSPHATE.///
	//However since screening_pose is scored before the CCD corrrectly position the chain_break phosphate atoms, ///////////////
	// the VIRTUAL_PHOSPHATE is needed to prevent artificial crashes. Parin Jan 28, 2010////////////////////////////////////
	if ( gap_size == 0 ) pose::add_variant_type_to_pose_residue ( screening_pose, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res + 1 );

	//Start the Rotamer Sampling
	Real /*current_score ( 0.0 ),*/ delta_rep_score ( 0.0 ), delta_atr_score ( 0.0 );
	core::Real delta_pucker;
	core::Real nu2_pucker;
	core::Real nu1_pucker;
	core::scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;


	RNA_LoopCloseSampler rna_loop_close_sampler ( moving_suite, cutpoint_closed );

	if ( finer_sampling_at_chain_closure_ ) {
		rna_loop_close_sampler.set_bin_size ( 10 );
	} else {
		rna_loop_close_sampler.set_bin_size ( 20 );
		rna_loop_close_sampler.set_epsilon_range ( 40 );
	}

	rna_loop_close_sampler.set_sample_only ( true );
	rna_loop_close_sampler.set_include_current ( true );
	rna_loop_close_sampler.set_scorefxn ( scorefxn_ );
	rna_loop_close_sampler.set_choose_random( choose_random_ );

	Size pucker_id;
	Size total_count = 0;
	std::cout << "Start Generating Rotamer ..." << std::endl;
	protocols::rna::RNA_IdealCoord const ideal_coord;

	// the possibilities for pucker_id. Probably should refactor slightly to be pucker_states (NORTH and SOUTH)... will soon create pucker_sampler.
	utility::vector1< Size > pucker_ids;
	for ( pucker_id = 0; pucker_id < 2; pucker_id ++ ) pucker_ids.push_back( pucker_id );

	for ( Size k = 0; k < pucker_ids.size(); k++ ){

		Size pucker_id = k;
		if ( choose_random_ ) pucker_id = numeric::random::random_element( pucker_ids );
		std::cout << "pucker_id = " << pucker_id << std::endl;

		if (use_phenix_geo_) {
			if ( pucker_id == 0 ) { //Sample North Pucker
				ideal_coord.apply( screening_pose, moving_res, true);
			} else { //Sample South Pucker
				ideal_coord.apply( screening_pose, moving_res, false);
			}
		} else {
			if ( pucker_id == 0 ) { //Sample North Pucker
				delta_pucker = rna_fitted_torsion_info.ideal_delta_north();
				nu2_pucker = rna_fitted_torsion_info.ideal_nu2_north();
				nu1_pucker = rna_fitted_torsion_info.ideal_nu1_north();
			} else { //Sample South Pucker
				delta_pucker = rna_fitted_torsion_info.ideal_delta_south();
				nu2_pucker = rna_fitted_torsion_info.ideal_nu2_south();
				nu1_pucker = rna_fitted_torsion_info.ideal_nu1_south();
			}

			screening_pose.set_torsion ( TorsionID ( moving_res , id::BB, 4 ) , delta_pucker );
			screening_pose.set_torsion ( TorsionID ( moving_res , id::CHI, 2 ) , nu2_pucker );
			screening_pose.set_torsion ( TorsionID ( moving_res , id::CHI, 3 ) , nu1_pucker );
		}

		//Loop Closure
		std::cout << "Entering Analytical Loop Closing" << std::endl;
		rna_loop_close_sampler.clear_all();
		rna_loop_close_sampler.apply ( screening_pose );
		std::cout << "Exiting Analytical Loop Closing" << std::endl;

		// note that in chose_random, loop_close_sampler should return 1 loop.
		for ( Size ii = 1; ii <= rna_loop_close_sampler.n_construct(); ++ii ) {
			rna_loop_close_sampler.fill_pose ( screening_pose, ii );

			// following initialization could be outside inner loop, right?
			//Sample CHI torsion angle
			BaseState base_state;
			if ( allow_syn_pyrimidine_ ) {
				base_state = BOTH;
			} else {
				base_state = ( core::scoring::rna::is_purine ( pose.residue ( moving_res ) ) ) ? BOTH : ANTI;
			}
			PuckerState pucker_state;
			if ( pucker_id == 0 ) {
				pucker_state = NORTH;
			} else {
				pucker_state = SOUTH;
			}
			StepWiseRNA_Base_Sugar_RotamerOP base_sugar_rotamer = new StepWiseRNA_Base_Sugar_Rotamer ( base_state, pucker_state, rna_fitted_torsion_info );

			base_sugar_rotamer->set_extra_syn_chi ( extra_syn_chi_rotamer_ );
			base_sugar_rotamer->set_extra_anti_chi ( extra_anti_chi_rotamer_ );
			base_sugar_rotamer->set_choose_random( choose_random_ );

			while ( base_sugar_rotamer->get_next_rotamer() ) {

				Real chi = base_sugar_rotamer->chi();
				screening_pose.set_torsion ( TorsionID ( moving_res , id::CHI, 1 ) ,chi );

				//Native RMSD Screening
				if ( native_rmsd_screen_ && get_native_pose() ) {
					if ( suite_rmsd ( *get_native_pose(), screening_pose, actually_moving_res, Is_prepend ) > ( native_screen_rmsd_cutoff_ ) ) continue;

					if ( rmsd_over_residue_list ( *get_native_pose(), screening_pose, job_parameters_, false ) > ( native_screen_rmsd_cutoff_ ) ) continue; //Oct 14, 2010

					count_data_.rmsd_count++;

					if ( verbose_ ) std::cout << "rmsd_count = " << count_data_.rmsd_count << " total count= " << count_data_.tot_rotamer_count << std::endl;
				}

				//Centroid Screening
				if ( centroid_screen_ ) {
					bool found_a_centroid_interaction_partner ( false );
					found_a_centroid_interaction_partner = base_centroid_screener_->Update_base_stub_list_and_Check_centroid_interaction ( screening_pose, count_data_ );

					if ( (gap_size > 0 || force_centroid_interaction_ )  && !found_a_centroid_interaction_partner ) continue;

					if ( !base_centroid_screener_->Check_that_terminal_res_are_unstacked() ) continue;
				}

				//VDW Screening//
				if ( !Full_atom_van_der_Waals_screening ( screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, gap_size, Is_internal ) ) continue;

				if ( ( user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ) && ( gap_size != 0 ) && ( Is_internal == false ) ) {
					//not well calibrated for chain_closure move and Is_internal move yet...
					//Also Right now the code doesn't correctly account for the phosphate of the prepend case...
					//Ok in that it ignore the virtual phosphate but doesn't take the 3' prime phosphate into accoutn June 13, 2010..
					if ( user_input_VDW_bin_screener_->VDW_rep_screen ( screening_pose, moving_res ) == false ) continue;

					count_data_.good_bin_rep_count++;
				}


				//Copy the Torsion angles of screen_pose to pose
				copy_suite_torsion ( pose, screening_pose, moving_res );
				if (use_phenix_geo_) {
					if ( pucker_id == 0 ) { //Sample North Pucker
						ideal_coord.apply( pose, moving_res, true);
					} else { //Sample South Pucker
						ideal_coord.apply( pose, moving_res, false);
					}
				}

				////////////////Add pose to pose_data_list/////
				if ( o2star_screen_ ) sample_o2star_hydrogen ( pose , pose_with_original_HO2star_torsion );

				std::stringstream ss;
				ss << "U" << total_count;
				std::string tag = ss.str();
				total_count++;
				Real current_score = Pose_selection_by_full_score ( pose_data_list, pose, tag );
				std::cout << "Size of pose_data_list: " << pose_data_list.size() << std::endl;

				std::cout << "CURRENT SCORE " << current_score <<  std::endl;

				if ( verbose_ ) {
					std::cout << tag <<  std::endl;
					Output_data ( silent_file_data, silent_file_, tag, true, pose, get_native_pose(), job_parameters_ );
				}

				if( choose_random_ && total_count > 0 ) break;

			} // chi (side chain)

			if( choose_random_ && total_count > 0 ) break;
		} // closed main chain

		if( choose_random_ && total_count > 0 ) break;
	} // pucker


	Output_title_text ( "Final sort and clustering" );
	std::cout << "before erasing.. pose_data_list= " << pose_data_list.size() << std::endl;
	std::sort ( pose_data_list.begin(), pose_data_list.end(), sort_criteria2 );
	cluster_pose_data_list ( pose_data_list );

	if ( pose_data_list.size() > num_pose_kept_ ) {
		pose_data_list.erase ( pose_data_list.begin() + num_pose_kept_, pose_data_list.end() );
	}

	std::cout << "after erasing.. pose_data_list= " << pose_data_list.size() << std::endl;
	pose_data_list_ = pose_data_list;
	//Hacky..temporary until we fix the reroot atom problem..This is just for calculating rmsd purposes... Apr 27 , 2010 Parin
	//////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "FINAL COUNTS" << std::endl;

	if ( gap_size <= 1 ) std::cout << " chain_closable_count= " << count_data_.chain_closable_count << std::endl;

	if ( gap_size == 0 ) {
		std::cout << " angle_n= " << count_data_.good_angle_count << " dist_n= " << count_data_.good_distance_count;
		std::cout << " chain_break_screening= " << count_data_.chain_break_screening_count << std::endl;
	}

	std::cout << "stack= " << count_data_.base_stack_count << " pair= " << count_data_.base_pairing_count;
	std::cout << " strict_pair_n= " << count_data_.strict_base_pairing_count;
	std::cout << " atr= " << count_data_.good_atr_rotamer_count;
	std::cout << " rep= " << count_data_.good_rep_rotamer_count;
	std::cout << " both= " << count_data_.both_count;
	std::cout << " bulge= " << count_data_.bulge_at_chain_closure_count;
	std::cout << " rmsd= " << count_data_.rmsd_count << " tot= " << count_data_.tot_rotamer_count << std::endl;
	std::cout << "Total time in StepWiseRNA_AnalyticalLoopCloseSampler: " << static_cast<Real> ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
}



//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::copy_suite_torsion ( core::pose::Pose & pose, core::pose::Pose const & ref_pose, core::Size const & suite_num ) {
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace protocols::rna;
	using namespace core::id;
	pose.set_torsion ( TorsionID ( suite_num - 1, id::BB, 5 ), ref_pose.torsion ( TorsionID ( suite_num - 1, id::BB, 5 ) ) ); //epsilon-1
	pose.set_torsion ( TorsionID ( suite_num - 1, id::BB, 6 ), ref_pose.torsion ( TorsionID ( suite_num - 1, id::BB, 6 ) ) ); //zeta-1
	pose.set_torsion ( TorsionID ( suite_num, id::BB, 1 ), ref_pose.torsion ( TorsionID ( suite_num, id::BB, 1 ) ) ); //alpha
	pose.set_torsion ( TorsionID ( suite_num, id::BB, 2 ), ref_pose.torsion ( TorsionID ( suite_num, id::BB, 2 ) ) ); //beta
	pose.set_torsion ( TorsionID ( suite_num, id::BB, 3 ), ref_pose.torsion ( TorsionID ( suite_num, id::BB, 3 ) ) ); //gamma
	pose.set_torsion ( TorsionID ( suite_num, id::BB, 4 ), ref_pose.torsion ( TorsionID ( suite_num, id::BB, 4 ) ) ); //delta
	pose.set_torsion ( TorsionID ( suite_num, id::BB, 5 ), ref_pose.torsion ( TorsionID ( suite_num, id::BB, 5 ) ) ); //epsilon
	pose.set_torsion ( TorsionID ( suite_num + 1, id::BB, 1 ), ref_pose.torsion ( TorsionID ( suite_num + 1, id::BB, 1 ) ) ); //alpha+1
	pose.set_torsion ( TorsionID ( suite_num + 1, id::BB, 2 ), ref_pose.torsion ( TorsionID ( suite_num + 1, id::BB, 2 ) ) ); //beta+1
	pose.set_torsion ( TorsionID ( suite_num + 1, id::BB, 3 ), ref_pose.torsion ( TorsionID ( suite_num + 1, id::BB, 3 ) ) ); //gamma+1
	pose.set_torsion ( TorsionID ( suite_num, id::CHI, 1 ), ref_pose.torsion ( TorsionID ( suite_num, id::CHI, 1 ) ) ); //chi
	pose.set_torsion ( TorsionID ( suite_num, id::CHI, 2 ), ref_pose.torsion ( TorsionID ( suite_num, id::CHI, 2 ) ) ); //nu2
	pose.set_torsion ( TorsionID ( suite_num, id::CHI, 3 ), ref_pose.torsion ( TorsionID ( suite_num, id::CHI, 3 ) ) ); //nu1
	return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::get_base_atr_rep_score ( core::pose::Pose const & pose, core::Real & base_atr_score, core::Real & base_rep_score ) {
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace ObjexxFCL;
	Size const working_moving_suite ( job_parameters_->working_moving_suite() );
	Size const working_moving_res ( job_parameters_->working_moving_res() );
	Size const nres = job_parameters_->working_sequence().size();
	///////////////////////////////Old_way////////////////////////////////////////////
	pose::Pose base_pose_screen = pose; //hard copy
//		std::cout << "BLAH_BLAH_BLAH May 5, 2010" << std::endl;
	if ( verbose_ ) base_pose_screen.dump_pdb ( "base_atr_rep_before.pdb" );

//		if(working_moving_suite>=nres) utility_exit_with_message( "working_moving_suite " + string_of(working_moving_suite) + " >= nres " + string_of(nres) );
	pose::add_variant_type_to_pose_residue ( base_pose_screen, "VIRTUAL_PHOSPHATE", working_moving_res ); //May 7...

	if ( ( working_moving_res + 1 ) <= nres ) {
		pose::add_variant_type_to_pose_residue ( base_pose_screen, "VIRTUAL_PHOSPHATE", working_moving_res + 1 ); //May 7...
	}

	// I think this should work... push apart different parts of the structure so that whatever fa_atr, fa_rep is left is due to "intra-domain" interactions.
	// Crap this doesn't work when building 2 or more nucleotides.
	Size const jump_at_moving_suite = make_cut_at_moving_suite ( base_pose_screen, working_moving_suite );
	kinematics::Jump j = base_pose_screen.jump ( jump_at_moving_suite );
	j.set_translation ( Vector ( 1.0e4, 0.0, 0.0 ) );
	base_pose_screen.set_jump ( jump_at_moving_suite, j );
	( *atr_rep_screening_scorefxn_ ) ( base_pose_screen );
	EnergyMap const & energy_map = base_pose_screen.energies().total_energies();
	base_atr_score = atr_rep_screening_scorefxn_->get_weight ( fa_atr ) * energy_map[ scoring::fa_atr ]; //
	base_rep_score = atr_rep_screening_scorefxn_->get_weight ( fa_rep ) * energy_map[ scoring::fa_rep ];
	std::cout << "base_rep= " << base_rep_score << " base_atr= " << base_atr_score << std::endl;

	if ( verbose_ ) base_pose_screen.dump_pdb ( "base_atr_rep_after.pdb" );

	////////////////////////////////////////////////////////////////////////////////////
	/*

	Size const nres = job_parameters_->working_sequence().size();
	utility::vector1< Size > const working_moving_res_list=job_parameters_->working_moving_res_list();


	utility::vector1< Size > partition_0_seq_num_list;
	utility::vector1< Size > partition_1_seq_num_list;


	bool const root_partition = partition_definition( rerooted_fold_tree.root() );

	for (Size seq_num=1; seq_num<=nres; seq_num++){
	//			if(Contain_seq_num(seq_num, working_moving_res_list)) continue; //Exclude working_moving_residues (the one being sampled..)

		if ( partition_definition( seq_num ) == 0){
		 	partition_0_seq_num_list.push_back( seq_num );
		}else if( partition_definition( seq_num ) == 1){
		 	partition_1_seq_num_list.push_back( seq_num );
		}else{
			utility_exit_with_message("seq_num " + string_of(seq_num) + " is not both in either partition!!" );
		}
	}

	sort_seq_num_list(partition_0_seq_num_list); //Low seq_num on the top of the list [1,2,3,4,5]
	sort_seq_num_list(partition_1_seq_num_list); //Low seq_num on the top of the list [1,2,3,4,5]

	*/
}


////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::initialize_scorefunctions() {
	initialize_common_scorefxns(scorefxn_, sampling_scorefxn_, atr_rep_screening_scorefxn_, chainbreak_scorefxn_, o2star_pack_scorefxn_);
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< pose_data_struct2 > &
StepWiseRNA_AnalyticalLoopCloseSampler::get_pose_data_list() {
	return pose_data_list_;
}



////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::Update_pose_data_list ( std::string const & tag, utility::vector1< pose_data_struct2 > & pose_data_list, pose::Pose const & current_pose, Real const & current_score ) const {
	bool add_pose_to_list = false;

	add_pose_to_list = ( current_score < current_score_cutoff_ );

	//The order of evaluation of the two expression in the if statement is important!
	if ( add_pose_to_list ) {
		if ( verbose_ ) {
			std::cout << "tag= " << tag << " current_score_cutoff_ " << current_score_cutoff_ << " score= " << current_score;
		}

		pose_data_struct2 current_pose_data;
		current_pose_data.pose_OP = new pose::Pose;
		( *current_pose_data.pose_OP ) = current_pose;
		current_pose_data.score = current_score;
		current_pose_data.tag = tag;

		if ( get_native_pose() )  {
			setPoseExtraScores ( *current_pose_data.pose_OP, "all_rms",
			                     core::scoring::rms_at_corresponding_heavy_atoms ( *current_pose_data.pose_OP, *get_native_pose() ) );
		}

		pose_data_list.push_back ( current_pose_data );

		if ( verbose_ ) std::cout << " pose_data_list.size= " << pose_data_list.size() << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseRNA_AnalyticalLoopCloseSampler::Pose_selection_by_full_score ( utility::vector1< pose_data_struct2 >& pose_data_list, pose::Pose & current_pose, std::string const & tag ) {
	using namespace core::scoring;
	Real const current_score = ( *sampling_scorefxn_ ) ( current_pose );
	Update_pose_data_list ( tag, pose_data_list, current_pose, current_score );

	if ( ( pose_data_list.size() == num_pose_kept_ * multiplier_ ) ) {
		std::sort ( pose_data_list.begin(), pose_data_list.end(), sort_criteria2 );
		cluster_pose_data_list ( pose_data_list );

		if ( pose_data_list.size() > num_pose_kept_ ) {
			pose_data_list.erase ( pose_data_list.begin() + num_pose_kept_, pose_data_list.end() );
			std::cout << "after erasing.. pose_data_list.size()= " << pose_data_list.size() << std::endl;
		} else {
			std::cout << "pose_data_list.size()= " << pose_data_list.size() << std::endl;
		}
	}

	return current_score;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dec 18, 2009...took off alot of optimization from this code since it is very fast (not rate limiting) anyways.
void
StepWiseRNA_AnalyticalLoopCloseSampler::cluster_pose_data_list ( utility::vector1< pose_data_struct2 >& pose_data_list ) {
	bool const Is_prepend ( job_parameters_->Is_prepend() );
	Size const actually_moving_res = job_parameters_->actually_moving_res();
	utility::vector1< bool > pose_state_list ( pose_data_list.size(), true );
	Size num_clustered_pose = 0;

	for ( Size i = 1; i <= pose_data_list.size(); i++ ) {
		if ( pose_state_list[i] == true ) {
			num_clustered_pose++;

			for ( Size j = i + 1; j <= pose_data_list.size(); j++ ) {
				Real rmsd;

				if ( PBP_clustering_at_chain_closure_ && job_parameters_->gap_size() == 0 ) { //new option Aug 15, 2010..include both phosphates in rmsd calculation at chain_break
					rmsd =	 phosphate_base_phosphate_rmsd ( ( *pose_data_list[i].pose_OP ), ( *pose_data_list[j].pose_OP ), actually_moving_res,  false /*ignore_virtual_atom*/ );
				} else {
					rmsd = suite_rmsd ( ( *pose_data_list[i].pose_OP ), ( *pose_data_list[j].pose_OP ), actually_moving_res, Is_prepend , false /*ignore_virtual_atom*/ );
				}

				bool const same_pucker = Is_same_ribose_pucker ( ( *pose_data_list[i].pose_OP ), ( *pose_data_list[j].pose_OP ), actually_moving_res );

				if ( rmsd < cluster_rmsd_ && ( same_pucker || !distinguish_pucker_ ) ) {
					pose_state_list[j] = false;

					if ( verbose_ ) {
						std::cout << "rmsd= " << rmsd << "  pose " << pose_data_list[j].tag << " is a neighbor of pose " << pose_data_list[i].tag;
						std::cout << " same_pucker= ";
						Output_boolean ( same_pucker );
						print_ribose_pucker_state ( " center_pucker= ", Get_residue_pucker_state ( ( *pose_data_list[i].pose_OP ), actually_moving_res ) );
						print_ribose_pucker_state ( " curr_pucker= ", Get_residue_pucker_state ( ( *pose_data_list[j].pose_OP ), actually_moving_res ) );
						std::cout << std::endl;
					}
				}
			}
		}
	}

	utility::vector1< pose_data_struct2> clustered_pose_data_list;

	for ( Size i = 1; i <= pose_data_list.size(); i++ ) {
		if ( pose_state_list[i] == true ) {
			clustered_pose_data_list.push_back ( pose_data_list[i] );
		}
	}

	pose_data_list = clustered_pose_data_list;

	////////May 12, get the score cutoff////////
	//check if pose_data_list size is equal to or exceed num_pose_kept_. Important to get score_cutoff here right after clustering.
	if ( pose_data_list.size() >= num_pose_kept_ ) {
		current_score_cutoff_ = pose_data_list[num_pose_kept_].score;
	} else {
		current_score_cutoff_ = 99999; //keep on adding pose to list if there are still not enough clusters
	}

	////////////////////////////////////////////
}
///////////////////////////////////////////////////////////////////////////////
//TEMPORARY
bool
StepWiseRNA_AnalyticalLoopCloseSampler::Full_atom_van_der_Waals_screening_REPLICATE ( pose::Pose & current_pose_screen,
    Real const & base_rep_score,
    Real const & base_atr_score,
    Real & delta_atr_score,
    Real & delta_rep_score,
    Size const & gap_size,
    bool const & Is_internal ) {
	using namespace core::scoring;

	if ( VDW_atr_rep_screen_ == false ) return true;

	bool close_chain = ( gap_size == 0 ) ? true : false;

	if ( close_chain && Is_internal ) return true; //Don't screen at all Mar 1, 2010

	( *atr_rep_screening_scorefxn_ ) ( current_pose_screen );
	EnergyMap const & energy_map = current_pose_screen.energies().total_energies();
	Real rep_score = atr_rep_screening_scorefxn_->get_weight ( fa_rep ) * energy_map[scoring::fa_rep];
	Real atr_score = atr_rep_screening_scorefxn_->get_weight ( fa_atr ) * energy_map[scoring::fa_atr];
	delta_rep_score = rep_score - base_rep_score;
	delta_atr_score = atr_score - base_atr_score;
	Real actual_rep_cutoff = rep_cutoff_; //defualt

	if ( close_chain ) actual_rep_cutoff = 200; //Parin's old parameter

	if ( close_chain && Is_internal ) actual_rep_cutoff = 200; //Bigger chunk..easier to crash...not using this right now.

	bool pass_rep_screen = false;

	if ( delta_rep_score < actual_rep_cutoff ) {
		pass_rep_screen = true;
	}

	bool pass_atr_rep_screen = false;

	if ( close_chain ) {
		pass_atr_rep_screen = pass_rep_screen;
	} else {
		if ( delta_atr_score < ( -1 ) && ( delta_rep_score + delta_atr_score ) < 0 ) pass_atr_rep_screen = true;
	}

	if ( pass_atr_rep_screen ) {
		if ( verbose_ ) {
			std::cout << " rep= " << delta_rep_score << " atr= " << delta_atr_score;
			std::cout << "  stack_n= " << count_data_.base_stack_count << " pair_n= " << count_data_.base_pairing_count;
			std::cout << "  strict_pair_n= " << count_data_.strict_base_pairing_count;
			std::cout << "  centroid_n= " << count_data_.pass_base_centroid_screen;
			std::cout << "  bin_rep_n= " << count_data_.good_bin_rep_count;
			std::cout << "  atr_n= " << count_data_.good_atr_rotamer_count;
			std::cout << "  rep_n= " << count_data_.good_rep_rotamer_count;
			std::cout << "  both= " << count_data_.both_count << " tot= " << count_data_.tot_rotamer_count;
			std::cout << "  closable= " << count_data_.chain_closable_count;
			std::cout << "  non_clash_ribose= " << count_data_.non_clash_ribose;
			std::cout << std::endl;
		}

		return true;
	} else {
		return false;
	}
}
///////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_AnalyticalLoopCloseSampler::Full_atom_van_der_Waals_screening ( pose::Pose & current_pose_screen,
    Real const & base_rep_score,
    Real const & base_atr_score,
    Real & delta_rep_score,
    Real & delta_atr_score,
    Size const & gap_size,
    bool const & Is_internal ) {
	using namespace core::scoring;
	using namespace ObjexxFCL;

	if ( VDW_atr_rep_screen_ == false ) return true;

	bool close_chain = ( gap_size == 0 ) ? true : false;

	if ( close_chain && Is_internal ) return true; //Don't screen at all Mar 1, 2010

	( *atr_rep_screening_scorefxn_ ) ( current_pose_screen );
	EnergyMap const & energy_map = current_pose_screen.energies().total_energies();
	Real rep_score = atr_rep_screening_scorefxn_->get_weight ( fa_rep ) * energy_map[scoring::fa_rep];
	Real atr_score = atr_rep_screening_scorefxn_->get_weight ( fa_atr ) * energy_map[scoring::fa_atr];
	delta_rep_score = rep_score - base_rep_score;
	delta_atr_score = atr_score - base_atr_score;

	if ( delta_rep_score < ( -0.01 ) ) {
		std::string const message = "delta_rep_score= " + string_of ( delta_rep_score ) + " rep_score= " + string_of ( rep_score ) + " base_rep_score= " + string_of ( base_rep_score );
		utility_exit_with_message ( "delta_rep_score<(-0.01), " + message );
	}

	if ( delta_atr_score > ( +0.01 ) ) {
		std::string const message = "delta_atr_score= " + string_of ( delta_atr_score ) + " atr_score= " + string_of ( atr_score ) + " base_atr_score= " + string_of ( base_atr_score );
		utility_exit_with_message ( "delta_atr_score>(+0.01), " + message );
	}

	Real actual_rep_cutoff = rep_cutoff_; //defualt

	if ( close_chain ) actual_rep_cutoff = 200; //Parin's old parameter

	if ( Is_internal ) actual_rep_cutoff = 200; //Bigger chunk..easier to crash (before May 4 used to be (close_chain && Is_internal) actual_rep_cutoff=200

	bool pass_rep_screen = false;

	if ( delta_rep_score < actual_rep_cutoff ) {
		pass_rep_screen = true;
		count_data_.good_rep_rotamer_count++;
	}

	if ( delta_atr_score < ( -1 ) || close_chain ) count_data_.good_atr_rotamer_count++;

	bool pass_atr_rep_screen = false;

	if ( close_chain ) {
		pass_atr_rep_screen = pass_rep_screen;
	} else if ( Is_internal ) {
		if ( delta_atr_score < ( -1 ) && ( delta_rep_score + delta_atr_score ) < ( actual_rep_cutoff - rep_cutoff_ ) ) pass_atr_rep_screen = true;
	} else {
		if ( delta_atr_score < ( -1 ) && ( delta_rep_score + delta_atr_score ) < 0 ) pass_atr_rep_screen = true;
	}

	if ( pass_atr_rep_screen ) {
		//	if((delta_atr_score<(-1)) && ((delta_rep_score+delta_atr_score) < 200) ) { //This causes about 5times more pose to pass the screen (50,000 poses vs 10,000 poses)
		count_data_.both_count++;

		if ( verbose_ ) {
			std::cout << " rep= " << delta_rep_score << " atr= " << delta_atr_score;
			std::cout << "  stack_n= " << count_data_.base_stack_count << " pair_n= " << count_data_.base_pairing_count;
			std::cout << "  strict_pair_n= " << count_data_.strict_base_pairing_count;
			std::cout << "  centroid_n= " << count_data_.pass_base_centroid_screen;
			std::cout << "  bin_rep_n= " << count_data_.good_bin_rep_count;
			std::cout << "  atr_n= " << count_data_.good_atr_rotamer_count;
			std::cout << "  rep_n= " << count_data_.good_rep_rotamer_count;
			std::cout << "  both= " << count_data_.both_count << " tot= " << count_data_.tot_rotamer_count << std::endl;
		}

		return true;
	} else {
		return false;
	}
}

////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::initialize_o2star_packer_task ( core::pose::Pose const & pose ) {
//		utility::vector1< core::Size > const working_moving_res_list= job_parameters_->working_moving_res_list();
	utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();
	utility::vector1< core::Size > const O2star_pack_seq_num = get_surrounding_O2star_hydrogen ( pose, working_moving_partition_pos, false /*verbose*/ );
//		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	o2star_pack_task_ =  pack::task::TaskFactory::create_packer_task ( pose );

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {
		if ( Contain_seq_num ( seq_num, O2star_pack_seq_num ) && pose.residue ( seq_num ).is_RNA() ) { //pack this residue!
			o2star_pack_task_->nonconst_residue_task ( seq_num ).and_extrachi_cutoff ( 0 );
			o2star_pack_task_->nonconst_residue_task ( seq_num ).or_ex4 ( true ); //extra O2star sampling
			o2star_pack_task_->nonconst_residue_task ( seq_num ).or_include_current ( true );
			// How about bump check?
		} else {
			o2star_pack_task_->nonconst_residue_task ( seq_num ).prevent_repacking();
		}
	}

	/*
			o2star_pack_task_ =  pack::task::TaskFactory::create_packer_task( pose );

			//Commented off becuase now this function is called by sample_o2star_hydrogen and this is computationally expensive? Parin S. Jan 28, 2010
			//		o2star_pack_task_->initialize_from_command_line();

			for (Size i = 1; i <= pose.total_residue(); i++) {
				if ( !pose.residue(i).is_RNA() ) continue;
				o2star_pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
				// Following could be useful...
				o2star_pack_task_->nonconst_residue_task(i).or_ex4( true ); //extra rotamers?? Parin S. Jan 28, 2010
				o2star_pack_task_->nonconst_residue_task(i).or_include_current( true );
			}
	*/
}

////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::sample_o2star_hydrogen ( core::pose::Pose & pose , core::pose::Pose & pose_with_original_HO2star_torsion ) {
	using namespace core::id;
	using namespace core::conformation;

	//reset the HO2star torsion to its starting value as to prevent randomness due to conformations sampling order...
	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {
		if ( pose.residue ( seq_num ).aa() == core::chemical::aa_vrt ) continue; //FCC

		//SML PHENIX conference
		if (basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value()){
			if ( !pose.residue ( seq_num ).is_RNA()) continue;
		}

		pose.set_torsion ( TorsionID ( seq_num, id::CHI, 4 ), pose_with_original_HO2star_torsion.torsion ( TorsionID ( seq_num, id::CHI, 4 ) ) );
	}

	//std::cout << "Packing 2'-OH ... ";
	//problem with bulge variant -- need to initialize_o2star_packer_task each time.
	initialize_o2star_packer_task ( pose );
	pack::rotamer_trials ( pose, *o2star_pack_scorefxn_, o2star_pack_task_ );
}

////////////////////////////////////////////////////////////////////////
std::string //silly function to convert to real to string
StepWiseRNA_AnalyticalLoopCloseSampler::create_torsion_value_string ( core::Real const & torsion_value ) const {
	using namespace ObjexxFCL;
	std::string torsion_string = "";
	core::Real	const principal_torsion = numeric::principal_angle_degrees ( torsion_value );
	Size const principal_torsion_SIZE = Size( std::abs ( principal_torsion + 0.00001 ) ); //0.00001 is to prevent random ambiguity if the torsion decimal value is exactly .0000 Oct 12, 2010

	if ( principal_torsion > 0 ) {
		torsion_string = "p" + lead_zero_string_of ( principal_torsion_SIZE, 3 );
	} else {
		torsion_string = "n" + lead_zero_string_of ( principal_torsion_SIZE, 3 );
	}

	return torsion_string;
}

////////////////////////////////////////////////////////////////////////
std::string //silly function used for appending the rotamer value to the tag
StepWiseRNA_AnalyticalLoopCloseSampler::create_rotamer_string ( core::pose::Pose const & pose ) const {
	std::string rotamer_tag = "";
	bool const Is_prepend ( job_parameters_->Is_prepend() );
	Size const moving_res ( job_parameters_->working_moving_res() );
	conformation::Residue const & five_prime_rsd = ( Is_prepend ) ? pose.residue ( moving_res ) : pose.residue ( moving_res - 1 );
	conformation::Residue const & three_prime_rsd = ( Is_prepend ) ?  pose.residue ( moving_res + 1 ) : pose.residue ( moving_res );
	rotamer_tag.append ( "_E" + create_torsion_value_string ( five_prime_rsd.mainchain_torsion ( 5 ) ) );
	rotamer_tag.append ( "_Z" + create_torsion_value_string ( five_prime_rsd.mainchain_torsion ( 6 ) ) );
	rotamer_tag.append ( "_A" + create_torsion_value_string ( three_prime_rsd.mainchain_torsion ( 1 ) ) );
	rotamer_tag.append ( "_B" + create_torsion_value_string ( three_prime_rsd.mainchain_torsion ( 2 ) ) );
	rotamer_tag.append ( "_G" + create_torsion_value_string ( three_prime_rsd.mainchain_torsion ( 3 ) ) );

	if ( Is_prepend ) {
		rotamer_tag.append ( "_D" + create_torsion_value_string ( five_prime_rsd.mainchain_torsion ( 4 ) ) );
		rotamer_tag.append ( "_C" + create_torsion_value_string ( five_prime_rsd.chi ( 1 ) ) );
	} else {
		rotamer_tag.append ( "_D" + create_torsion_value_string ( three_prime_rsd.mainchain_torsion ( 4 ) ) );
		rotamer_tag.append ( "_C" + create_torsion_value_string ( three_prime_rsd.chi ( 1 ) ) );
	}

	return rotamer_tag;
}
////////////////////////////////////////////////////////////////////////

std::string
StepWiseRNA_AnalyticalLoopCloseSampler::create_tag ( std::string const prestring, StepWiseRNA_RotamerGenerator_WrapperOP const & rotamer_generator ) const {
	using namespace ObjexxFCL;
	std::string tag = prestring;

	for ( Size list_position = rotamer_generator->rotamer_generator_list_size(); list_position >= 2; list_position-- ) { //For dinucleotide
		tag.append ( "_" + lead_zero_string_of ( rotamer_generator->group_rotamer ( list_position ), 4 ) );
	}

	tag.append ( "_" + lead_zero_string_of ( rotamer_generator->group_rotamer ( 1 ), 4 ) );
	tag.append ( "_" + lead_zero_string_of ( rotamer_generator->subgroup_rotamer ( 1 ), 5 ) );
	return tag;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_silent_file ( std::string const & silent_file ) {
	silent_file_ = silent_file;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_fast ( bool const & setting ) {
	fast_ = setting;

	if ( fast_ ) num_pose_kept_ = 2;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_medium_fast ( bool const & setting ) {
	medium_fast_ = setting;

	if ( medium_fast_ ) num_pose_kept_ = 20;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_native_rmsd_screen ( bool const & setting ) {
	native_rmsd_screen_ = setting;
	//if (native_rmsd_screen_) num_pose_kept_ = 20; Aug 16 2010..Parin S. For python_rebuild_suite.py
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_verbose ( bool const & setting ) {
	verbose_ = setting;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_o2star_screen ( bool const & setting ) {
	o2star_screen_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_output_filename ( std::string const & output_filename ) {
	output_filename_ = output_filename;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_scorefxn ( core::scoring::ScoreFunctionOP const & scorefxn ) {
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
core::io::silent::SilentFileDataOP &
StepWiseRNA_AnalyticalLoopCloseSampler::silent_file_data() {
	return sfd_;
}


//////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::output_pose_data_list ( std::string const final_sampler_output_silent_file ) const {
	using namespace core::io::silent;

	if ( verbose_ == false ) { //consistency check Apr 3, 2010
		utility_exit_with_message ( "verbose_==false, but StepWiseRNA_AnalyticalLoopCloseSampler::output_pose_data_list is still called?!" );
	}

	SilentFileData silent_file_data;

	for ( Size n = 1; n <= pose_data_list_.size(); n++ ) {
		Output_data ( silent_file_data, final_sampler_output_silent_file, pose_data_list_[n].tag, false, * ( pose_data_list_[n].pose_OP ), get_native_pose(), job_parameters_ );
	}
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_base_centroid_screener ( StepWiseRNA_BaseCentroidScreenerOP & screener ) {
	base_centroid_screener_ = screener;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_AnalyticalLoopCloseSampler::set_user_input_VDW_bin_screener ( StepWiseRNA_VDW_Bin_ScreenerOP const & user_input_VDW_bin_screener )
{
	user_input_VDW_bin_screener_ = user_input_VDW_bin_screener;
}


void
StepWiseRNA_AnalyticalLoopCloseSampler::set_cluster_rmsd ( Real const & setting ) {
	cluster_rmsd_ = setting;
	std::cout << "Set cluster_rmsd to " << cluster_rmsd_ << std::endl;
}

void
StepWiseRNA_AnalyticalLoopCloseSampler::set_num_pose_kept ( core::Size const & num_pose_kept ) {
	num_pose_kept_ = num_pose_kept ;
	// PARIN THIS DOES NOT MAKE ANY SENSE??
	//If fast_ or native_rmsd_scren the num_pose_kept_= 10 or 20 respectively Parin Feb 13, 2010
	set_native_rmsd_screen ( native_rmsd_screen_ );
	set_fast ( fast_ );
}

}
}
}
