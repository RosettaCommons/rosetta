// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Minimizer
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_Bin_Screener.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <core/conformation/Conformation.hh>
#include <protocols/rna/RNA_LoopCloser.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>
#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.rna_stepwise_rna_minimizer" ) ;

namespace protocols {
namespace swa {
namespace rna {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseRNA_Minimizer::StepWiseRNA_Minimizer(
																							 utility::vector1 <pose_data_struct2> const & pose_data_list,
																							 StepWiseRNA_JobParametersCOP & job_parameters ):
		pose_data_list_( pose_data_list ),
		job_parameters_( job_parameters ),
		silent_file_( "silent_file.txt" ),
		verbose_(false),
		native_screen_(false),
		native_screen_rmsd_cutoff_(3.0),
		centroid_screen_(true),
		perform_o2star_pack_(true), 
		output_before_o2star_pack_(false), 
		skip_minimize_(false), 
		num_pose_minimize_(99999999999), 
		minimize_and_score_sugar_(true) 
  {
		set_native_pose( job_parameters_->working_native_pose() );
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseRNA_Minimizer::~StepWiseRNA_Minimizer()
  {}

	/////////////////////
	std::string
	StepWiseRNA_Minimizer::get_name() const {
		return "StepWiseRNA_Minimizer";
	}

	//////////////////////////////////////////////////////////////////////////


	void
	StepWiseRNA_Minimizer::apply( core::pose::Pose & pose ) {

		using namespace core::scoring;
		using namespace core::scoring::rna;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::optimization;

		//Assume that the both the harmonic constraint and chemical::CUTPOINT_LOWER/UPPER for the CCD is already set up....(actually there is no harmonic constraints May 14, Parin 2010)...
		//Assume that the virtual phosphate is correctly

		Output_title_text("Enter StepWiseRNA_Minimizer::apply");

		bool const perform_minimizer_run=true; //this is for debugging!

		Output_boolean(" verbose_=", verbose_ ); std::cout << std::endl;
		Output_boolean(" native_screen_=", native_screen_ ); std::cout << std::endl;
		std::cout << " native_screen_rmsd_cutoff_=" << native_screen_rmsd_cutoff_ << std::endl;
		Output_boolean(" centroid_screen_=", centroid_screen_ ); std::cout << std::endl;
		Output_boolean(" perform_o2star_pack_=", perform_o2star_pack_); std::cout << std::endl;
		Output_boolean(" output_before_o2star_pack_=", output_before_o2star_pack_ ); std::cout << std::endl;
		std::cout << " (Upper_limit) num_pose_minimize_=" << num_pose_minimize_<< " pose_data_list.size()= " << pose_data_list_.size() <<  std::endl;
		Output_boolean(" minimize_and_score_sugar_=", minimize_and_score_sugar_ ); std::cout << std::endl;
		Output_boolean(" user_inputted_VDW_bin_screener_=", user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ); std::cout << std::endl;
		Output_seq_num_list(" working_global_sample_res_list=", job_parameters_->working_global_sample_res_list() ); 
		Output_boolean(" skip_minimize_=", skip_minimize_); std::cout << std::endl;
		Output_boolean(" perform_minimizer_run=", perform_minimizer_run); std::cout << std::endl;

		clock_t const time_start( clock() );

		Size const moving_res(  job_parameters_->working_moving_res() );
		Size const actually_moving_res(  job_parameters_->actually_moving_res() );
		Size const gap_size(  job_parameters_->gap_size() );
		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
	
		utility::vector1< pose_data_struct2 > minimized_pose_data_list;

		RNA_LoopCloser rna_loop_closer;

		SilentFileData silent_file_data;

	  AtomTreeMinimizer minimizer;
    float const dummy_tol( 0.00000025);
    bool const use_nblist( true );
    MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
    options.nblist_auto_update( true );

		if(pose_data_list_.size()==0){
			 std::cout << "pose_data_list_.size()==0, early exit from StepWiseRNA_Minimizer::apply" << std::endl;
			 return;
		}

		pose::Pose dummy_pose = (*pose_data_list_[1].pose_OP);
		Size const nres =  dummy_pose.total_residue();
		//////////////////////////////May 30 2010//////////////////////////////////////////////////////////////
		utility::vector1< core::Size > const & working_fixed_res= job_parameters_->working_fixed_res();
		utility::vector1< core::Size > working_minimize_res;
		bool o2star_pack_verbose=true;

		for(Size seq_num=1; seq_num<=pose.total_residue(); seq_num++ ){
			if(Contain_seq_num(seq_num, working_fixed_res)) continue;
			working_minimize_res.push_back(seq_num);
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////

		//Check scorefxn
		std::cout << "check scorefxn" << std::endl;
		scorefxn_->show( std::cout, dummy_pose );

		if( move_map_list_.size()==0 ) move_map_list_=Get_default_movemap( dummy_pose );

		if(minimize_and_score_sugar_==false){
			std::cout << "WARNING: minimize_and_score_sugar_ is FALSE, freezing all DELTA, NU_1 and NU_2 torsions." << std::endl;
			for(Size n=1; n<=move_map_list_.size(); n++){
				Freeze_sugar_torsions( move_map_list_[n] , nres );
			}
		}

		for(Size n=1; n<=move_map_list_.size(); n++){
			std::cout << "OUTPUT movemap #" << n << std::endl;
			Output_movemap(move_map_list_[n], nres);
		}

		if( gap_size == 0 ) {
      //			std::cout << "five_prime_chain_break_res= " << five_prime_chain_break_res << std::endl;
			//			Add_harmonic_chainbreak_constraint(pose, five_prime_chain_break_res ); //This should be already added during the sampling phrase
		}

		// put in terms to keep chainbreak closed during minimize.
		// even if there is not chainbreak closure, let's have this term (will be 0.0 ),
		//  to keep the number of scoreterms the same in the final outfile.
    //		scorefxn_->set_weight( linear_chainbreak, 1.0 ); //Mod this on May 10,2010..set linear_chain_break from score weight file instead...

 		for(Size i=1; i<=pose_data_list_.size(); i++){

			if(i>num_pose_minimize_){
				std::cout << "WARNING MAX num_pose_minimize_(" << num_pose_minimize_ << ") EXCEEDED, EARLY BREAK." << std::endl;
				break;
			}

			std::cout << "Minimizing pose_num= " << i << std::endl;

 			std::string tag=pose_data_list_[i].tag;
			pose=(*pose_data_list_[i].pose_OP); //This actually create a hard copy.....need this to get the graphic working..

			if(verbose_ && output_before_o2star_pack_){
			 	tag[0]='B'; 
				std::cout << tag <<  std::endl;
				(*scorefxn_)(pose); //Score pose to ensure that that it has a score to be output, Parin May 31, 2010.
				Output_data(silent_file_data, silent_file_+"_before_o2star_pack", tag, false , pose, get_native_pose(), job_parameters_);
			}

			Remove_virtual_O2Star_hydrogen(pose);

			if(perform_o2star_pack_) o2star_minimize(pose, scorefxn_, get_surrounding_O2star_hydrogen(pose, working_minimize_res, o2star_pack_verbose) );
				
			if(verbose_ && !output_before_o2star_pack_){
			 	tag[0]='B'; //B for before minimize
				std::cout << tag <<  std::endl;
				(*scorefxn_)(pose); //Score pose to ensure that that it has a score to be output
				Output_data(silent_file_data, silent_file_+"_before_minimize", tag, false , pose, get_native_pose(), job_parameters_);
			}

			if(skip_minimize_) continue;


 			///////Minimization/////////////////////////////////////////////////////////////////////////////
			if (gap_size == 0){
				rna_loop_closer.apply( pose, five_prime_chain_break_res ); //This doesn't do anything if rna_loop_closer was already applied during sampling stage...May 10,2010
			
				if(verbose_){ //May 10, 2010....just for consistency check...
				 	tag[0]='C'; //B for before minimize
					std::cout << tag <<  std::endl;
					(*scorefxn_)(pose); //Score pose to ensure that that it has a score to be output
					Output_data(silent_file_data, silent_file_+"_after_loop_closure_before_minimize", tag, false , pose, get_native_pose(), job_parameters_);
				}
			}			

			for(Size round=1; round<=move_map_list_.size(); round++){
				core::kinematics::MoveMap mm=move_map_list_[round];

				if(perform_minimizer_run) minimizer.run( pose, mm, *(scorefxn_), options );
				if(perform_o2star_pack_) o2star_minimize(pose, scorefxn_, get_surrounding_O2star_hydrogen(pose, working_minimize_res, o2star_pack_verbose) );

				//Mod this out on May 10, extra CDD just lead to more randomness...extra CCD can make chain_break score better but also can make torsion score at chain_break as well as 
				//fa_rep worst since not sensitive to these score elements..plus X5 linear_chain_break score should ensure that chain close well by this point...May 10..
				//Umm...change my decision...include an extra minimizer after CCD instead..This decision is base on observation that linear_chain_break  score is significantly 
				//better with extra round of CCD....The extra minimizer and o2star_pack is meant to help remove any new clash introduced by the extra CCD...
				//Checked minimizing time and found that this extra minimizer doesn't significantly slow down code....642.05 sec (extra minimize) vs. 603.01 sec (without extra minimize)...
				//This time is from running 1zih 108 poses...minimize two nucleotides....May 10, 2010...

				if( gap_size == 0 ){ 
					Real mean_dist_err = rna_loop_closer.apply( pose, five_prime_chain_break_res );
					std::cout << "mean_dist_err (round= " << round << " ) = " <<  mean_dist_err << std::endl;
					 //new May 10, 2010, temporary mod out to reproduce May 4th results/////
					if(perform_o2star_pack_) o2star_minimize(pose, scorefxn_, get_surrounding_O2star_hydrogen(pose, working_minimize_res, o2star_pack_verbose) ); 
					if(perform_minimizer_run) minimizer.run( pose, mm, *(scorefxn_), options );
					/////////////////////////////////////////////////////////////////////////
				}

			}

			////////////////////////////////Final screening //////////////////////////////////////////////
			//Sept 19, 2010 screen for weird delta/sugar conformation. Sometime at the chain closure step, minimization will severely screw up the pose sugar
			//if the chain is not a closable position. These pose will have very bad score and will be discarded anyway. The problem is that the clusterer check
			//for weird delta/sugar conformation and call utility_exit_with_message() if one is found (purpose being to detect bad silent file conversion)
			//So to prevent code exit, will just screen for out of range sugar here.
			if (gap_size == 0){		

				//Oct 28, 2010 ...ok check for messed up nu1 and nu2 as well. Note that check_for_messed_up_structure() also check for messed up delta but the range is smaller than below.
				if(check_for_messed_up_structure(pose, tag)==true){
					std::cout << "gap_size == 0, " << tag << " discarded: messed up structure " << std::endl;
					continue;
				}

				conformation::Residue const & five_prime_rsd = pose.residue(five_prime_chain_break_res);
				Real const five_prime_delta = numeric::principal_angle_degrees(five_prime_rsd.mainchain_torsion( DELTA ));

	

				if( (five_prime_rsd.has_variant_type("VIRTUAL_RNA_RESIDUE")==false) && (five_prime_rsd.has_variant_type("VIRTUAL_RIBOSE")==false)){
					if( (five_prime_delta>1.0 && five_prime_delta<179.00)==false ){

						std::cout << "gap_size == 0, " << tag << " discarded: five_prime_chain_break_res= " << five_prime_chain_break_res << " five_prime_CB_delta= " << five_prime_delta << " is out of range " << std::endl;
						continue;
					}
				}	

				conformation::Residue const & three_prime_rsd = pose.residue(five_prime_chain_break_res+1);
				Real const three_prime_delta = numeric::principal_angle_degrees(three_prime_rsd.mainchain_torsion( DELTA ));

				if( (three_prime_rsd.has_variant_type("VIRTUAL_RNA_RESIDUE")==false) && (three_prime_rsd.has_variant_type("VIRTUAL_RIBOSE")==false)){
					if( (three_prime_delta>1.0 && three_prime_delta<179.00)==false ){

						std::cout << "gap_size == 0, " << tag << " discarded: three_prime_chain_break_res= " << (five_prime_chain_break_res+1) << " three_prime_CB_delta= " << three_prime_delta << " is out of range " << std::endl;
						continue;
					}
				}
			}

			
			if(centroid_screen_){
				if ( !base_centroid_screener_->Update_base_stub_list_and_Check_that_terminal_res_are_unstacked( pose, true /*reinitialize*/ ) ){
					std::cout << tag << " discarded: fail Check_that_terminal_res_are_unstacked	" << std::endl;	 
					continue;
				}
				if ( !Pose_screening(pose, tag, silent_file_data) ){
					std::cout << tag << " discarded: fail Pose_screening " << std::endl;	 
					continue;
				}
			}

			//March 20, 2011..This is neccessary though to be consistent with FARFAR. Cannot have false low energy state that are due to empty holes in the structure.
			//Feb 22, 2011. Warning this is inconsistent with the code in SAMPLERER:
			//In Both standard and floating base sampling, user_input_VDW_bin_screener_ is not used for gap_size==0 or internal case.
			//This code is buggy for gap_size==0 case IF VDW_screener_pose contains the residue right at 3' or 5' of the building_res
			//Internal case should be OK.
			if(user_input_VDW_bin_screener_->user_inputted_VDW_screen_pose() ){ 
				if(i==1) std::cout << "user_inputted_VDW_screen_pose=true" << std::endl;

				//Convert to using physical_pose instead of bin for screening (as default) on June 04, 2011
				utility::vector1 < core::Size > const & working_global_sample_res_list=job_parameters_->working_global_sample_res_list();
				bool const pass_physical_pose_VDW_rep_screen=user_input_VDW_bin_screener_->VDW_rep_screen_with_act_pose( pose, working_global_sample_res_list, true /*local verbose*/);

				if( pass_physical_pose_VDW_rep_screen==false){
					std::cout << tag << " discarded: fail physical_pose_VDW_rep_screen" << std::endl;	 
				 	continue;
				}

			}
			//////////////////////////////////////////////////////////////////////////////////////////////

			pose_data_struct2 pose_data;

			pose_data.tag=tag;
			pose_data.score=(*scorefxn_)(pose);
			pose_data.pose_OP=new pose::Pose;
			(*pose_data.pose_OP)=pose;
			minimized_pose_data_list.push_back(pose_data);

			// might as well output as we go along -- no longer clustering in between.

			tag[ 0 ] = 'M';
			(*scorefxn_)(pose); //Score pose to ensure that that it has a score to be output
			Output_data(silent_file_data, silent_file_, tag, false , pose, get_native_pose(), job_parameters_);

			std::cout << "Total time in StepWiseRNA_Minimizer: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

  	}

		Output_title_text("Exit StepWiseRNA_Minimizer::apply");

//Clustering the pose will reduces the work load of the master node but is not necessary, Dec 20, 2009 Parin.
//		Output_title_text("Final_Sort_and_Clustering");
//		std::sort(pose_data_list.begin(), pose_data_list.end(), sort_citeria);
//		cluster_pose_data_list(pose_data_list);
//		std::cout << "minimized_pose_data_list.size() = " << minimized_pose_data_list.size() << std::endl;


// 		for(Size i=1; i<=minimized_pose_data_list.size(); i++){
// 			pose_data_struct2 & pose_data=minimized_pose_data_list[i];
// 			pose::Pose & pose=*pose_data.pose_OP;

// 			pose_data.tag[0]='M';  //M for minimized, these are the poses that will be clustered by master node.
// 			std::cout << pose_data.tag;
// 			Output_data(silent_file_data, silent_file_, pose_data.tag, false, pose, get_native_pose(), actually_moving_res, Is_prepend);

// 		}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}


	void
	StepWiseRNA_Minimizer::Freeze_sugar_torsions(core::kinematics::MoveMap & mm, Size const nres) const {

		 using namespace core::id;

		 std::cout << "Freeze pose sugar torsions, nres=" << nres << std::endl;
		 	 	 	 
		 for(Size i=1; i<=nres; i++){	 
			
				mm.set( TorsionID( i  , id::BB,  4 ), false ); //delta_i
				mm.set( TorsionID( i  , id::CHI, 2 ), false ); //nu2_i	
				mm.set( TorsionID( i  , id::CHI, 3 ), false );	//nu1_i	 
	
		 }
	}


	//Cannot pass pose in as constant due to the setPoseExtraScores function
	bool
	StepWiseRNA_Minimizer::Pose_screening(pose::Pose & pose, std::string tag, core::io::silent::SilentFileData & silent_file_data) const{

		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::optimization;


		Size const moving_res(  job_parameters_->working_moving_res() );
		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const gap_size( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();

		bool pass_screen=true;

		tag[0]='S'; //S for screening

		if( gap_size == 1){
			if( !Check_chain_closable(pose, five_prime_chain_break_res, gap_size) ){
				tag[0]='F'; //Mark as unclosable loop failure so that will not be use to rebuild next (last) residue
				pass_screen=false;
			}
		}

		if(native_screen_ && get_native_pose()){	//Before have the (&& !Is_chain_break condition). Parin Dec 21, 2009
	
			Real const rmsd= suite_rmsd(*get_native_pose(), pose,  moving_res, Is_prepend, true /*ignore_virtual_atom*/); 
			Real const loop_rmsd=	 rmsd_over_residue_list( *get_native_pose(), pose, job_parameters_, true /*ignore_virtual_atom*/);

			if(rmsd > native_screen_rmsd_cutoff_ || loop_rmsd > native_screen_rmsd_cutoff_){ 
				tag[0]='R'; //Mark as RMSD failure
				std::cout << tag << " discarded: fail native_rmsd_screen. rmsd= " << rmsd << " loop_rmsd= " << loop_rmsd << " native_screen_rmsd_cutoff_= " << native_screen_rmsd_cutoff_ << std::endl;
				pass_screen=false;
			}
		}

		if(verbose_){
			std::cout << tag <<  std::endl;
			(*scorefxn_)(pose); //Score pose to ensure that that it has a score to be output
			Output_data(silent_file_data, silent_file_ + "_screen", tag, false, pose, get_native_pose(), job_parameters_);
  		}


		return pass_screen;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1 <core::kinematics::MoveMap>
	StepWiseRNA_Minimizer::Get_default_movemap( core::pose::Pose const & pose ) const{

		Size const nres( pose.total_residue() );

		utility::vector1 <core::kinematics::MoveMap> move_map_list;

		core::kinematics::MoveMap mm;
		Figure_out_moving_residues( mm, pose );

		// Allow sugar torsions to move again (RD 01/31/2010), now
		// that rotamers have been pre-optimized for ribose closure, and
		// sugar_close is turned back on.
		//		Freeze_sugar_torsions(mm, nres); //Freeze the sugar_torsions!

		move_map_list.push_back(mm);
		return move_map_list;
	}

	////////////////////////////////////////////////////////////////////////////////////
	// This is similar to code in RNA_Minimizer.cc
	void
	StepWiseRNA_Minimizer::Figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose ) const
  {

		using namespace core::id;
		using namespace core::scoring::rna;
		using namespace ObjexxFCL;

		Size const nres( pose.total_residue() );

		utility::vector1< core::Size > const & fixed_res( job_parameters_->working_fixed_res() );

		ObjexxFCL::FArray1D< bool > allow_insert( nres, true );
		for (Size i = 1; i <= fixed_res.size(); i++ ) allow_insert( fixed_res[ i ] ) = false;

		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );

		for  (Size i = 1; i <= nres; i++ )  {

			utility::vector1< TorsionID > torsion_ids;

			for ( Size rna_torsion_number = 1; rna_torsion_number <= NUM_RNA_MAINCHAIN_TORSIONS; rna_torsion_number++ ) {
				torsion_ids.push_back( TorsionID( i, id::BB, rna_torsion_number ) );
			}
			for ( Size rna_torsion_number = 1; rna_torsion_number <= NUM_RNA_CHI_TORSIONS; rna_torsion_number++ ) {
				torsion_ids.push_back( TorsionID( i, id::CHI, rna_torsion_number ) );
			}


			for ( Size n = 1; n <= torsion_ids.size(); n++ ) {

				TorsionID const & torsion_id  = torsion_ids[ n ];

				id::AtomID id1,id2,id3,id4;
				bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
				if (fail) continue; //This part is risky, should also rewrite...

				// Dec 19, 2010..Crap there is a mistake here..should have realize this earlier...
				//Should allow torsions at the edges to minimize...will have to rewrite this. This effect the gamma and beta angle of the 3' fix res.

				// If there's any atom that is in a moving residue by this torsion, let the torsion move.
				//  should we handle a special case for cutpoint atoms? I kind of want those all to move.
				if ( !allow_insert( id1.rsd() ) && !allow_insert( id2.rsd() ) && !allow_insert( id3.rsd() )  && !allow_insert( id4.rsd() ) ) continue;
				mm.set(  torsion_id, true );

			}

		}

		std::cout << "pose.fold_tree().num_jump()= " << pose.fold_tree().num_jump() << std::endl;

		for (Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
			Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

			if ( allow_insert( jump_pos1 ) || allow_insert( jump_pos2 ) ) 	 mm.set_jump( n, true );
			std::cout << "jump_pos1= " << jump_pos1 << " jump_pos2= " << jump_pos2 << " mm.jump= "; Output_boolean(allow_insert( jump_pos1 ) || allow_insert( jump_pos2 ) );  std::cout << std::endl;
		}

	}


	////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Minimizer::set_move_map_list(utility::vector1 <core::kinematics::MoveMap> const & move_map_list){
		move_map_list_ = move_map_list;
	}


  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseRNA_Minimizer::set_silent_file( std::string const & silent_file){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Minimizer::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		scorefxn_ = scorefxn;
	}

  //////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP &
	StepWiseRNA_Minimizer::silent_file_data(){
		return sfd_;
	}

	//////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Minimizer::set_base_centroid_screener( StepWiseRNA_BaseCentroidScreenerOP & screener ){
		base_centroid_screener_ = screener;
	}



}
}
}
