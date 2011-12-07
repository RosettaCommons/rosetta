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
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseJobParameters.hh>
#include <protocols/swa/StepWiseUtil.hh>
//#include <protocols/swa/rna/StepWiseRNA_Dinucleotide_Sampler_Util.hh>

//////////////////////////////////
#include <core/types.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/BinaryRNASilentStruct.hh>

// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <protocols/rna/RNA_LoopCloser.hh>
#include <core/io/pdb/pose_io.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

// AUTO-REMOVED #include <numeric/xyz.functions.hh> // APL TEMP

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>
#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
// AUTO-REMOVED #include <math.h>
// AUTO-REMOVED #include <stdlib.h>

//Auto Headers
#include <core/scoring/rms_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using namespace core;
using core::Real;
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
  StepWiseRNA_ResidueSampler::StepWiseRNA_ResidueSampler( StepWiseJobParametersCOP & job_parameters ):
		job_parameters_( job_parameters ),
	  sfd_( new core::io::silent::SilentFileData ),
		scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ),// can be replaced from the outside
		silent_file_( "silent_file.txt" ),
		output_filename_( "data.txt"),
		bin_size_( 20 ),
		rep_cutoff_( 4.0 ),
		num_pose_kept_( 108 ),
		multiplier_(2), //Sort and cluster poses when the number of pose is pose_data_list exceed multiplier*num_pose_kept,
		cluster_rmsd_(0.5),
		verbose_( false ),
		native_rmsd_screen_(false),
		native_screen_rmsd_cutoff_(2.0),
		o2star_screen_( true ),
		use_green_packer_( false ),
		allow_bulge_at_chainbreak_( false ),
		fast_( false )
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


	//////////////////////////////////////////////////////////////////////////
	bool
	sort_criteria(pose_data_struct2  pose_data_1, pose_data_struct2 pose_data_2) {  //This function used to be call sort_criteria2
		return (pose_data_1.score < pose_data_2.score);
	}
	////////////////////////////////////////////////////////////////////////////

	/* Now can build both single nucleotide and dinucleotide Parin Feb 6, 2010*/
	void
	StepWiseRNA_ResidueSampler::apply( core::pose::Pose &  pose ) {

		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::id;



		std::ofstream outfile;
		SilentFileData silent_file_data;
		outfile.open(output_filename_.c_str());

		// FOLLOWING IS NOT PERFECT. In the general case multiple residues might move.
		Size const moving_res(  job_parameters_->working_moving_res_list()[1] ); // corresponds to user input.
		Size const moving_suite(  job_parameters_->working_moving_suite_list()[1] ); // dofs betweeen this value and value+1 actually move.

		bool const Is_prepend(  job_parameters_->Is_prepend() ); // if true, moving_suite+1 is fixed. Otherwise, moving_suite is fixed.
		bool const Is_internal(  job_parameters_->Is_internal() ); // no cutpoints before or after moving_res.
		Size const actually_moving_res( job_parameters_->actually_moving_res() );
		Size const gap_size( job_parameters_->gap_size()); /* If this is zero or one, need to screen or closable chain break */
		utility::vector1 < core::Size > const & moving_positions = job_parameters_->moving_pos();
		Size const num_nucleotides(  job_parameters_->working_moving_res_list().size() );

		if(num_nucleotides!=1 && num_nucleotides!=2){
			utility_exit_with_message( "num_nucleotides!=1 and num_nucleotides!=2" );
 		}

		bool const Is_dinucleotide=(num_nucleotides==2);

		std::cout << " NUM_NUCLEOTIDES= " <<  num_nucleotides << std::endl;
		Output_boolean(" IS_DINUCLEOTIDE= ", Is_dinucleotide); std::cout << std::endl;
		std::cout << " GAP SIZE " << gap_size << std::endl;
		std::cout << " MOVING RES " << moving_res << std::endl;
		std::cout << " MOVING SUITE " << moving_suite << std::endl;
		Output_boolean(" PREPEND ", Is_prepend ); std::cout << std::endl;
		Output_boolean(" INTERNAL ", Is_internal); std::cout << std::endl;



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Pose const pose_save = pose;
		pose = pose_save; //this recopy is useful for triggering graphics.

		//////////////// Sets up scorefunctions for a bunch of different screens ////////////////////////////////////////
		initialize_scorefunctions();

		//////////////////////////////Setup Bulge residue///////////////////////////////////////////////////////////////////////////

		Size const bulge_moving_suite = (Is_prepend) ? moving_suite+1: moving_suite-1;
		Size const bulge_moving_res= (Is_prepend) ? moving_res+1 : moving_res-1;

		if(Is_dinucleotide){
			if(bulge_moving_res!=job_parameters_->working_moving_res_list()[2]){
				utility_exit_with_message( "bulge_moving_res!=working_moving_res_list[2]" );
			}
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", bulge_moving_res);
		}

		/////////////////////////////// O2star sampling/virtualization //////////////////////////
		Pose working_pose = pose;
		Add_virtual_O2Star_hydrogen( working_pose );
		working_pose.set_torsion( TorsionID( moving_res, id::CHI, 4 ), 0 );  //This torsion is not sampled. Arbitary set to zero to prevent randomness

		// if o2star_screen, pose 2'-OH will be sampled!
		if ( o2star_screen_ ) {
			if ( use_green_packer_ ){
				initialize_o2star_green_packer( pose );
			} else {
				initialize_o2star_packer_task( pose );
			}
		} else {
			// Otherwise, virtualize the 2-OH.
			pose = working_pose;
		}

		/////////////////////////////////////////Setup Base-stack/Base-pairing screening/////////////////////////////////////////
		bool const base_centroid_screening = ( base_centroid_screener_ != 0 );

		///////////////////////////////////////Setup chainbreak_screening//////////////////////////////////////////////////////
		//This assumes that the both the harmonic constraint and chemical::CUTPOINT_LOWER/UPPER is setup for the CCD is already set up
 		RNA_LoopCloser rna_loop_closer;

		Size const five_prime_chain_break_res = job_parameters_->first_chain_break_res();

		if ( gap_size == 0 )	{
			// This part is needed even though we have a "linear_chainbreak" term in the rna_minimizer  /////////////////////////////
			// since the harmonic chain_break score is used for chain_break_screening.	/////////////////////////////////////////////
			// Since there are two score terms (distance and angle constraints) in the harmonic term, ///////////////////////////////
			// the screening is more "robust" than just using the single term in the linear_chain_break , Parin Jan 2, 2010 /////////
			std::cout << "five_prime_chain_break_res= " << five_prime_chain_break_res << std::endl;
		 	Add_harmonic_chainbreak_constraint(working_pose, five_prime_chain_break_res );
		}

		//Does CCD work if the Virtual phosphate variant type exist?. Rhiju mentioned that he will fix this.

		//In this latest version, chain_break_screening_pose DOES NOT have a VIRTUAL_PHOSPHATE if gap_size == 0. ///////////////////
		//This is CONSISTENT with Parin's old code /////////////////////////////////////////////////////////////////////////////////
		//It neither hurt nor help for chain_break_screening_pose to have a VIRTUAL_PHOSPHATE. One of the concern was that with ////
		//a VIRTUAL_PHOSPHATE, the CCD code might but work but Parin try running with the virtual phosphate and "seem" to work  ////
		//Parin Jan 28, 2010 ///////////////////////////////////////////////////////////////////////////////////////////////////////

		pose::Pose chain_break_screening_pose=working_pose; //Hard copy


	  //////////////////////////////////////////Setup Atr_rep_screening/////////////////////////////////////////////////

		pose::Pose screening_pose = working_pose; //Hard copy

		//Necessary for the case where gap_size == 0. In this case, Pose_setup does not automatically create a VIRTUAL_PHOSPHATE.///
		//However since screening_pose is scored before the CCD corrrectly position the chain_break phosphate atoms, ///////////////
		// the VIRTUAL_PHOSPHATE is needed to prevent artificial crashes. Parin Jan 28, 2010////////////////////////////////////
		if (gap_size == 0) pose::add_variant_type_to_pose_residue( screening_pose, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res+1 );

		Pose base_pose_screen = screening_pose;

		// I think this should work... push apart different parts of the structure so that
		// whatever fa_atr, fa_rep is left is due to "intra-domain" interactions.
		Size const jump_at_moving_suite = make_cut_at_moving_suite( base_pose_screen, moving_suite); //Might want to change so the cut is at the first bulge_moving_suite Jan 30, 2010
		kinematics::Jump j = base_pose_screen.jump( jump_at_moving_suite );
		j.set_translation( Vector( 1.0e4, 0.0, 0.0 ) );
		base_pose_screen.set_jump( jump_at_moving_suite, j );

		(*atr_rep_screening_scorefxn_)(base_pose_screen);

		EnergyMap const & energy_map=base_pose_screen.energies().total_energies();
		Real base_atr_score = atr_rep_screening_scorefxn_->get_weight(fa_atr) * energy_map[ scoring::fa_atr ];
		Real base_rep_score = atr_rep_screening_scorefxn_->get_weight(fa_rep) * energy_map[ scoring::fa_rep ];
		std::cout << "base_rep= " << base_rep_score << " base_atr= " << base_atr_score << std::endl;

		base_pose_screen.dump_pdb( "start_extend.pdb" ); //for debugging..mod out longer in Rhiju version

		// Note: Is_prepend should be generalized to include two more possibilities:
		//   we are creating a dinucleotide from scratch --> sample sugar/chi for both moving_res and
		bool sample_sugar_and_base1( false ), sample_sugar_and_base2( false );
		if ( !Is_internal  ) {
			if ( Is_prepend ) {
				sample_sugar_and_base1 = true;
			} else {
				sample_sugar_and_base2 = true;
			}
		}


		utility::vector1< core::Size > const working_moving_suite_list(  job_parameters_->working_moving_suite_list() );

		StepWiseRNA_RotamerGenerator_WrapperOP rotamer_generator = new StepWiseRNA_RotamerGenerator_Wrapper( pose,
																																																				 working_moving_suite_list,
																																																				 sample_sugar_and_base1,
																																																				 sample_sugar_and_base2);
		rotamer_generator->set_fast( fast_ );

		utility::vector1< pose_data_struct2 > pose_data_list;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MAIN LOOP --> rotamer sampling.
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		Real current_score( 0.0 ), delta_rep_score( 0.0), delta_atr_score( 0.0 ); //Before use to just create these real time. Parin S, Jan 28, 2010

		clock_t const time_start( clock() ); //Remember to move this back to the top of the function after done testing Parin S. Jan 28, 2010

		while( rotamer_generator->has_another_rotamer() ){


			utility::vector1< Torsion_Info > current_rotamer = rotamer_generator->get_next_rotamer();
			apply_rotamer( screening_pose, current_rotamer);
			count_data_.tot_rotamer_count++;

			//				if(fast_ && count_data_.both_count==100) break;

			std::string tag=create_tag("U", rotamer_generator);

			if (native_rmsd_screen_ && get_native_pose()){
				// Following does not look right -- what if pose and native_pose are not superimposed corectly? -- rhiju
				if( suite_rmsd(*get_native_pose(), screening_pose, actually_moving_res, Is_prepend)>(native_screen_rmsd_cutoff_)) continue;
				count_data_.rmsd_count++;
				if(verbose_) std::cout << "rmsd_count = " << count_data_.rmsd_count << " total count= " << count_data_.tot_rotamer_count << std::endl;
			}

			bool found_a_centroid_interaction_partner( false );

			if( base_centroid_screening ){

				//Reminder note of independency: Important that base_stub_list is updated even in the case where gap_size == 0 (and bulge is allowed) ////
				// since base_stub_list is used later below in the chain_break_screening section Jan 28, 2010 Parin S. ///////////////////////////////////

				// updates base_stub_list.
				found_a_centroid_interaction_partner = base_centroid_screener_->Update_base_stub_list_and_Check_centroid_interaction( screening_pose);

				if ( gap_size > 0 && !found_a_centroid_interaction_partner ) continue;
				/*if its the "last residue", allow for bulge*/

				// Note that is does not update base_stub_list. To do that, use Update_base_stub_list_and_Check_that_terminal_res_are_unstacked
				if ( !base_centroid_screener_->Check_that_terminal_res_are_unstacked() ) continue;

			}

			//////////////////////////////////////////////////////////////////////////////////////////
			///////////////Chain_break_screening -- distance cut                     /////////////////
			//////////////////////////////////////////////////////////////////////////////////////////
			if ( gap_size <= 1 ){
				if ( !Check_chain_closable(screening_pose, five_prime_chain_break_res, gap_size ) ) continue;
				count_data_.chain_closable_count++;
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			/////////////// Full_atom_van_der_Waals_screening                        /////////////////
			//////////////////////////////////////////////////////////////////////////////////////////
			//////*In Parin's old code, fa_rep_cut was much more lenient for during chain_break_closure (pass screen if delta_rep<10) Parin Jan 28, 2010
			bool apply_fa_atr_cut = ( gap_size > 0 );
			if ( !Full_atom_van_der_Waals_screening( screening_pose, base_rep_score, base_atr_score, delta_rep_score, delta_atr_score, apply_fa_atr_cut ) ) continue;

			//////////////////////////////////////////////////////////////////////////////////////////
			// Almost ready to actually score pose.
			//////////////////////////////////////////////////////////////////////////////////////////
			apply_rotamer( pose, current_rotamer );

			//////////////////////////////////////////////////////////////////////////////////////////
			///////////////Chain_break_screening -- CCD closure /////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////
			bool 	bulge_added( false );
			if ( gap_size == 0 /*really need to close it!*/ ){

				//Problem is that chain_break backbone torsions is currently to initialized to value after CCD of last pose July 27, 2009
				//Need to fix this!! But this maybe beneficial since chain break torsion value of the
				// different sampling pose might be similar??
				apply_rotamer(chain_break_screening_pose, current_rotamer );
				if( ! Chain_break_screening( chain_break_screening_pose, chainbreak_scorefxn_) ) continue;

				// Need to be very careful here -- do CCD torsions ever overlap with pose torsions?
				Copy_CCD_torsions( pose, chain_break_screening_pose);
				if (!found_a_centroid_interaction_partner && moving_positions.size() < 2) apply_bulge_variant( pose, bulge_added, delta_atr_score ); /*further cut on atr, inside*/
				}
			/////////////////////////////////////////////////////////////////////////////////////////////

			//				if(fast_) screening_pose.dump_pdb( tag + ".pdb" );

			////////////////Add pose to pose_data_list if pose have good score////////////////////////////////////////////

			if ( o2star_screen_ ) sample_o2star_hydrogen( pose );

			Real current_score=Pose_selection_by_full_score(pose_data_list, pose, tag, sampling_scorefxn_);

			if(verbose_){
				std::cout << tag <<  std::endl;
				Output_data(silent_file_data, silent_file_, tag, true, pose, get_native_pose(), actually_moving_res, Is_prepend);

			}
		} //while( rotamer_generator->has_another_rotamer() )

		// just make sure something gets output.
		if ( pose_data_list.size() == 0 )  Pose_selection_by_full_score(pose_data_list, pose, "U_000", sampling_scorefxn_);


		Output_title_text("Final sort and clustering");
		std::sort(pose_data_list.begin(), pose_data_list.end(), sort_criteria);
		if ( cluster_rmsd_ > 0.0 ) cluster_pose_data_list(pose_data_list);
		if( pose_data_list.size()>num_pose_kept_ ) pose_data_list.erase(pose_data_list.begin()+num_pose_kept_, pose_data_list.end());
		std::cout<< "after erasing.. pose_data_list= " << pose_data_list.size() << std::endl;

		pose_data_list_=pose_data_list;

		pose = pose_save;
		outfile.close();

		std::cout << "FINAL COUNTS" << std::endl;
		if( gap_size <= 1) std::cout << " chain_closable_count= " << count_data_.chain_closable_count << std::endl;
		if( gap_size == 0 ) std::cout << " angle_n= " << count_data_.good_angle_count << " dist_n= " << count_data_.good_distance_count << std::endl;
		std::cout << "stack= " << count_data_.base_stack_count << " pair= " << count_data_.base_pairing_count;
		std::cout << " atr= " << count_data_.good_atr_rotamer_count;
		std::cout << " rep= " << count_data_.good_rep_rotamer_count;
		std::cout << " both= " << count_data_.both_count;
		std::cout << " rmsd= " << count_data_.rmsd_count << " tot= " << count_data_.tot_rotamer_count << std::endl;
//		std::cout << "bulge_rotamer_list_size= " << bulge_rotamer_list_size << std::endl;
		std::cout << "Total time in StepWiseRNA_ResidueSampler: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	}

	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::initialize_scorefunctions(){
		using namespace core::scoring;

		///////////////////////////////////////////////////////////////////
		// Bare minimum to check for contact (fa_atr) but not clash (fa_rep)
		atr_rep_screening_scorefxn_ =  new ScoreFunction;
		atr_rep_screening_scorefxn_->set_weight( fa_atr  , 0.23 );
		atr_rep_screening_scorefxn_->set_weight( fa_rep  , 0.12 );

		///////////////////////////////////////////////////////////////////
		chainbreak_scorefxn_ =  new ScoreFunction;
		chainbreak_scorefxn_->set_weight( angle_constraint, 1.0 );
		chainbreak_scorefxn_->set_weight( atom_pair_constraint, 1.0 );

		////////////////////Setup sampling scoring//////////////////////////////////////////////////////////////////////////////
    //1. Want to increase fa_rep during the minimization phase but want to keep it at 0.12 during the sample phase
	  //2. Sugar scoring is always turned off during sampling stage since it screw up pose selection. (TURN IT BACK ON: RD 01/31/2010)
		//3. Harmonic and Linear Chain_break scoring is always turned off during sampling stage
		sampling_scorefxn_ = scorefxn_->clone();

		//		sampling_scorefxn_->set_weight( rna_sugar_close, 0.0 );
		sampling_scorefxn_->set_weight( fa_rep, 0.12 );
		//		sampling_scorefxn_->set_weight( angle_constraint, 0.0 );
		//		sampling_scorefxn_->set_weight( atom_pair_constraint, 0.0 );
		sampling_scorefxn_->set_weight( linear_chainbreak, 0.0);
		if ( scorefxn_->get_weight( rna_bulge ) > 0.0 ) sampling_scorefxn_->set_weight( rna_bulge, 0.5 /*This is totally arbitrary*/);

		///////////////////////////////////////////////////////////////////
		o2star_pack_scorefxn_ = new ScoreFunction;
		// Each of the following terms have been pretty optimized for the packer (trie, etc.)
		o2star_pack_scorefxn_->set_weight( fa_atr, sampling_scorefxn_->get_weight( fa_atr ) );
		o2star_pack_scorefxn_->set_weight( fa_rep, sampling_scorefxn_->get_weight( fa_rep ) );
		o2star_pack_scorefxn_->set_weight( hbond_lr_bb_sc, sampling_scorefxn_->get_weight( hbond_lr_bb_sc ) );
		o2star_pack_scorefxn_->set_weight( hbond_sr_bb_sc, sampling_scorefxn_->get_weight( hbond_sr_bb_sc ) );
		o2star_pack_scorefxn_->set_weight( hbond_sc, sampling_scorefxn_->get_weight( hbond_sc ) );
		o2star_pack_scorefxn_->set_energy_method_options( sampling_scorefxn_->energy_method_options() );
		// note that geom_sol is not optimized well --> replace with lk_sol for now.
		o2star_pack_scorefxn_->set_weight( fa_sol, sampling_scorefxn_->get_weight( lk_nonpolar ) );
		// Note that: rna_torsion, rna_sugar_close, fa_stack not optimized -- also irrelevant for 2'-OH sampling.

		// just a comparison. This is extremely slow. Would need to implement trie
		//  for geom_sol, lk_nonpolar, and hackelec... Not too hard, but I don't feel like doing it now.
		//o2star_pack_scorefxn_ = sampling_scorefxn_->clone();

	}

	////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< pose_data_struct2 > &
	StepWiseRNA_ResidueSampler::get_pose_data_list(){
		return pose_data_list_;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::Copy_CCD_torsions(pose::Pose & pose, pose::Pose const & template_pose) const {

 		using namespace core::chemical;
		using namespace core::conformation;
	  using namespace core::id;

		Size const five_prime_res = job_parameters_->first_chain_break_res();
		Size const three_prime_res = five_prime_res+1;

		conformation::Residue const & lower_res=template_pose.residue(five_prime_res);
		conformation::Residue const & upper_res=template_pose.residue(three_prime_res);

		//Even through there is the chain_break, alpha of 3' and epl and gamma of 5' should be defined due to the existence of the upper and lower variant type atoms.

		for(Size n=1; n<=3; n++){ //alpha, beta, gamma of 3' res
			pose.set_torsion( TorsionID( three_prime_res, id::BB,  n ), upper_res.mainchain_torsion(n) );
		}

		for(Size n=5; n<=6; n++){ //epsilon and zeta of 5' res
			pose.set_torsion( TorsionID( five_prime_res, id::BB,  n ), lower_res.mainchain_torsion(n) );
		}
	}


	////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_ResidueSampler::Chain_break_screening(  pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn ){

		using namespace core::scoring;

 		static protocols::rna::RNA_LoopCloser rna_loop_closer;
		Size const five_prime_res = job_parameters_->first_chain_break_res();
		//		Real const mean_dist_err=rna_loop_closer.apply( chain_break_screening_pose, five_prime_res);
		rna_loop_closer.apply( chain_break_screening_pose, five_prime_res);

		(*chainbreak_scorefxn)(chain_break_screening_pose);

		scoring::EMapVector & energy_map= chain_break_screening_pose.energies().total_energies();
		Real angle_score = energy_map[scoring::angle_constraint];
		Real distance_score = energy_map[scoring::atom_pair_constraint];


		if(angle_score<5) count_data_.good_angle_count++;
		if(distance_score<5) count_data_.good_distance_count++;
		if((angle_score<5) && (distance_score<5)){
			count_data_.both_count++;
			if(verbose_){
				//				std::cout << " C5_O3= " << C5_O3_distance << " C5_O3_n= " << count_data_.C5_O3_distance_count;
				std::cout << "  angle= " << angle_score << " dist= " << distance_score;
				std::cout << " angle_n= " << count_data_.good_angle_count;
				std::cout << " dist_n= " << count_data_.good_distance_count;
				std::cout << " both= " << count_data_.both_count;
				std::cout << " tot= " << count_data_.tot_rotamer_count << std::endl;
			}
			return true;
		} else {
			return false;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_ResidueSampler::apply_bulge_variant( core::pose::Pose & pose, bool & bulge_added, Real const & delta_atr_score ){

		static Real const atr_cutoff_for_bulge( 0.0 );
		// FOLLOWING IS NOT GREAT -- requires there to be only one nucleotide, I think.
		assert( job_parameters_->working_moving_res_list().size() == 1 );
		Size const moving_res(  job_parameters_->working_moving_res_list()[1] ); // corresponds to user input.

		bulge_added = false;

		if ( allow_bulge_at_chainbreak_ ) {
			if ( delta_atr_score >= atr_cutoff_for_bulge ) {

				if (verbose_) std::cout << "delta_atr " << delta_atr_score << " passes cutoff for bulge " << atr_cutoff_for_bulge << std::endl;
				//std::cout << "Is BULGE already there? " << pose.residue( moving_res ).has_variant_type( "BULGE" ) << std::endl;;
				pose::add_variant_type_to_pose_residue( pose, "BULGE", moving_res );
				bulge_added = true;

			} else {

				if (verbose_) std::cout << "delta_atr " << delta_atr_score << " DOES NOT PASS cutoff for bulge " << atr_cutoff_for_bulge << std::endl;

			}
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::Update_pose_data_list(std::string const & tag, utility::vector1< pose_data_struct2 > & pose_data_list, pose::Pose const & current_pose, Real const & current_score) const{

		//The order of evaluation of the two expression in the if statement is important!
		if(pose_data_list.size() < num_pose_kept_ || current_score < pose_data_list[num_pose_kept_].score) {

			if(verbose_){
				std::cout << "tag= " << tag;

				if(pose_data_list.size() >= num_pose_kept_) std::cout << " cutoff score= " << pose_data_list[num_pose_kept_].score;

				std::cout<< " score= " << current_score;
			}

			pose_data_struct2 current_pose_data;
			current_pose_data.pose_OP=new pose::Pose;
			(*current_pose_data.pose_OP)=current_pose;
			current_pose_data.score = current_score;
			current_pose_data.tag=tag;

			if ( get_native_pose())  {
				setPoseExtraScores( *current_pose_data.pose_OP, "all_rms",
														core::scoring::rms_at_corresponding_heavy_atoms( *current_pose_data.pose_OP, *get_native_pose() ) );
			}

			pose_data_list.push_back(current_pose_data);
			if(verbose_) std::cout << " pose_data_list.size= " << pose_data_list.size() << std::endl;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	StepWiseRNA_ResidueSampler::Pose_selection_by_full_score(utility::vector1< pose_data_struct2 >& pose_data_list, pose::Pose & current_pose, /*utility::vector1<Real> const & current_rotamer,*/ std::string const & tag, core::scoring::ScoreFunctionOP const & scorefxn) const{

		using namespace core::scoring;

		Real const current_score=(*scorefxn)(current_pose);

		//Very bad score pose...don't bother to add to list. Before Jan 28,2009 used to be if(current_score>-1) return 99.99;
		//		if(current_score>99.99) return 99.99;

		Update_pose_data_list(tag, pose_data_list, current_pose, current_score);

		if((pose_data_list.size()==num_pose_kept_*multiplier_)){
			std::sort(pose_data_list.begin(), pose_data_list.end(), sort_criteria);
			cluster_pose_data_list(pose_data_list);
			if(pose_data_list.size()>num_pose_kept_){
				pose_data_list.erase(pose_data_list.begin()+num_pose_kept_, pose_data_list.end());
				std::cout<< "after erasing.. pose_data_list.size()= " << pose_data_list.size() << std::endl;
			}else{
				std::cout<< "pose_data_list.size()= " << pose_data_list.size() << std::endl;
			}
		}


		return current_score;

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Dec 18, 2009...took off alot of optimization from this code since it is very fast (not rate limiting) anyways.
	void
	StepWiseRNA_ResidueSampler::cluster_pose_data_list(utility::vector1< pose_data_struct2 >& pose_data_list) const{

		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const actually_moving_res = job_parameters_->actually_moving_res();

		utility::vector1< bool > pose_state_list( pose_data_list.size(), true );

		Size num_clustered_pose=0;

		for(Size i=1; i<=pose_data_list.size(); i++){

			if(pose_state_list[i]==true){
				num_clustered_pose++;
				for(Size j=i+1; j<=pose_data_list.size(); j++){

					Real rmsd = suite_rmsd( (*pose_data_list[i].pose_OP), (*pose_data_list[j].pose_OP), actually_moving_res, Is_prepend );

					if(rmsd < cluster_rmsd_){
						pose_state_list[j]=false;
						if(verbose_) std::cout << "rmsd= " << rmsd << "  pose" << pose_data_list[j].tag << " is a neighbor of pose " << pose_data_list[i].tag << std::endl;
					}
				}
			}
		}


		utility::vector1< pose_data_struct2> clustered_pose_data_list;

		for(Size i=1; i<=pose_data_list.size(); i++) {
			if(pose_state_list[i]==true){
				clustered_pose_data_list.push_back(pose_data_list[i]);
			}
		}

		pose_data_list=clustered_pose_data_list;

	}

	///////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_ResidueSampler::Full_atom_van_der_Waals_screening(
																																pose::Pose & current_pose_screen,
																																Real const & base_rep_score,
																																Real const & base_atr_score,
																																Real & delta_atr_score,  //Added by Rhiju. Parin S Jan 28, 2009
																																Real & delta_rep_score,  //Added by Rhiju. Parin S Jan 28, 2009
																																bool const apply_fa_atr_cut /* = true */ ){

		using namespace core::scoring;

		(*atr_rep_screening_scorefxn_)(current_pose_screen);

		EnergyMap const & energy_map = current_pose_screen.energies().total_energies();

		Real rep_score = atr_rep_screening_scorefxn_->get_weight(fa_rep) * energy_map[scoring::fa_rep];
		Real atr_score = atr_rep_screening_scorefxn_->get_weight(fa_atr) * energy_map[scoring::fa_atr];

		delta_rep_score=rep_score-base_rep_score;
		delta_atr_score=atr_score-base_atr_score;

//		static bool const verbose_( false ); THIS SHOULD NOT BE HERE!!!! Jan 28, 2010 Parin S

		if( delta_rep_score<rep_cutoff_ ) count_data_.good_rep_rotamer_count++;

		if( delta_atr_score<(-1) || !apply_fa_atr_cut ) count_data_.good_atr_rotamer_count++;

		//		std::cout << "base atr score: " << base_atr_score << "     new atr score: " << atr_score << std::endl;

		if( ( delta_atr_score<(-1) || !apply_fa_atr_cut ) && ( (delta_rep_score+delta_atr_score) < 0) ) {
			//	if((delta_atr_score<(-1)) && ((delta_rep_score+delta_atr_score) < 200) ) { //This causes about 5times more pose to pass the screen (50,000 poses vs 10,000 poses)
			count_data_.both_count++;
			if ( verbose_ ) {
				std::cout << " rep= " << delta_rep_score << " atr= " << delta_atr_score;
				std::cout << "  stack_n= " << count_data_.base_stack_count << " pair_n= " << count_data_.base_pairing_count;
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
	// Could also use the GreenPacker --> potentially can do full pack at very little cost.
	void
	StepWiseRNA_ResidueSampler::initialize_o2star_packer_task( core::pose::Pose const & pose ){

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

	}

	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ResidueSampler::initialize_o2star_green_packer( core::pose::Pose & pose )
	{
		using namespace protocols::moves;
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		o2star_green_packer_ = new protocols::simple_moves::GreenPacker;

		ObjexxFCL::FArray1D< bool > const & partition_definition = job_parameters_->partition_definition();
		bool const root_partition = partition_definition( pose.fold_tree().root() );

		Size const nres = pose.total_residue();
		protocols::simple_moves::UserDefinedGroupDiscriminatorOP user_defined_group_discriminator( new protocols::simple_moves::UserDefinedGroupDiscriminator);
		utility::vector1< Size > group_ids;

		Size current_group = 0;
		Size spectator_group = 1;

		for (Size i = 1; i <= nres; i++ ) {

			if ( partition_definition( i ) != root_partition ) {
				current_group = 0;
				std::cout << "GREENPACKER SAMPLER " << i << std::endl;
			} else {
				std::cout << "GREENPACKER SPECTATOR   " << i <<  " --> group " << spectator_group << std::endl;
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
		for (Size i = 1; i <= nres; i++) {
			if ( !pose.residue(i).is_RNA() ) continue;
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
	StepWiseRNA_ResidueSampler::sample_o2star_hydrogen( core::pose::Pose & pose ){

		//std::cout << "Packing 2'-OH ... ";
		if ( use_green_packer_ ) {
			o2star_green_packer_->apply( pose );
		} else {

			//problem with bulge variant -- need to initialize_o2star_packer_task each time.
			initialize_o2star_packer_task( pose );

			pack::rotamer_trials( pose, *o2star_pack_scorefxn_, o2star_pack_task_ );
			//pack::pack_rotamers( pose, *o2star_pack_scorefxn_, o2star_pack_task_ );

		}
		//std::cout << " done. " << std::endl;

	}

	////////////////////////////////////////////////////////////////////////
	std::string
	StepWiseRNA_ResidueSampler::create_tag(std::string const prestring, StepWiseRNA_RotamerGenerator_WrapperOP const & rotamer_generator){

		using namespace ObjexxFCL;

		std::string tag=prestring;

		for(Size list_position=rotamer_generator->rotamer_generator_list_size(); list_position>=2; list_position--){ //For dinucleotide
			tag.append("_" + lead_zero_string_of(rotamer_generator->group_rotamer(list_position), 4));
		}

		tag.append("_" + lead_zero_string_of(rotamer_generator->group_rotamer(1), 3));
		tag.append("_" + lead_zero_string_of(rotamer_generator->subgroup_rotamer(1), 5));

		return tag;
	}

/*
	//Best way to make this robust is to tokenize the tag.
	std::string
	StepWiseRNA_ResidueSampler::create_tag(std::string const & prestring, Size const & group_rotamer, Size const & subgroup_rotamer, std::string const & old_tag){

		using namespace ObjexxFCL;

		std::string tag=old_tag;
		if(tag=="") tag.append(prestring);
		else tag[0]=prestring[0];

		tag.append("_");
		tag.append(lead_zero_string_of(group_rotamer, 3));
		tag.append("_");
		tag.append(lead_zero_string_of(subgroup_rotamer, 5));
		//			tag.append("_");
		//		 	tag.append(lead_zero_string_of(fine_rotamer_count, 4));
		//  tag.append("_");
		//  tag.append(Get_one_letter_name(Current_Rebuild_Unit::rebuild_unit.rebuild_residue.name));
		//  tag.append(lead_zero_string_of(Current_Rebuild_Unit::rebuild_unit.rebuild_residue.seq_num, 2)); //Made this change on Sep 22, 2009.
		//		 	std::cout << "tag= " << tag << std::endl;
		return tag;
	}

	std::string
	StepWiseRNA_ResidueSampler::create_tag(std::string const & prestring, Size const bulge_rotamer_ID, Size const & group_rotamer, Size const & subgroup_rotamer, std::string const & old_tag){

		using namespace ObjexxFCL;

		std::string tag=old_tag;
		if(tag=="") tag.append(prestring);
		else tag[0]=prestring[0];

		tag.append("_");
		tag.append(lead_zero_string_of(bulge_rotamer_ID, 4));
		tag.append("_");
		tag.append(lead_zero_string_of(group_rotamer, 3));
		tag.append("_");
		tag.append(lead_zero_string_of(subgroup_rotamer, 5));

		return tag;
	}
*/

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_fast( bool const & setting ){
    fast_ = setting;
//		if (fast_) num_pose_kept_ = 40;
		if (fast_) num_pose_kept_ = 10;
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_native_rmsd_screen( bool const & setting ){
    native_rmsd_screen_ = setting;
		if (native_rmsd_screen_) num_pose_kept_ = 20;
  }
  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_verbose( bool const & setting ){
    verbose_ = setting;
  }


  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_o2star_screen( bool const & setting ){
    o2star_screen_ = setting;
  }


  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_ResidueSampler::set_allow_bulge_at_chainbreak( bool const & setting ){
    allow_bulge_at_chainbreak_ = setting;
  }


  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseRNA_ResidueSampler::set_output_filename( std::string const & output_filename){
		output_filename_=output_filename;
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
	StepWiseRNA_ResidueSampler::output_pose_data_list( std::string const & silent_file ) const{
		using namespace core::io::silent;
		Output_pose_data_list( pose_data_list_, silent_file );

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
	}

	void
	StepWiseRNA_ResidueSampler::set_num_pose_kept( core::Size const & num_pose_kept ){
		num_pose_kept_=num_pose_kept ;

		// PARIN THIS DOES NOT MAKE ANY SENSE??
		set_native_rmsd_screen(native_rmsd_screen_); // does nothin??
		set_fast(fast_); // does nothing??
	}

}
}
}
