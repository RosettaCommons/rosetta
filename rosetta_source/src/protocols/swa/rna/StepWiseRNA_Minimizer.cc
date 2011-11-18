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
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseJobParameters.hh>
#include <protocols/swa/StepWiseJobParameters.fwd.hh>
#include <protocols/swa/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
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

#ifdef WIN32
	#include <ctime>
#endif


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

static basic::Tracer TR( "protocols.swa.stepwise_rna_minimizer" ) ;

namespace protocols {
namespace swa {
namespace rna {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseRNA_Minimizer::StepWiseRNA_Minimizer(
																							 utility::vector1 <pose_data_struct2> const & pose_data_list,
																							 StepWiseJobParametersCOP & job_parameters ):
		pose_data_list_( pose_data_list ),
		job_parameters_( job_parameters ),
		silent_file_( "silent_file.txt" ),
		verbose_(false),
		native_screen_(false),
		rmsd_cutoff_(3.0),
		minimize_rounds_( 1 )
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
		using namespace core::scoring::constraints;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::optimization;

		//Assume that the both the harmonic constraint and chemical::CUTPOINT_LOWER/UPPER for the CCD is already set up.
		//Assume that the virtual phosphate is correctly setup

		clock_t const time_start( clock() );

		ConstraintSetOP cst_set = pose.constraint_set()->clone();

		Size const actually_moving_res(  job_parameters_->actually_moving_res() );
		Size const gap_size(  job_parameters_->gap_size() );
		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const five_prime_chain_break_res = job_parameters_->first_chain_break_res();

		utility::vector1< pose_data_struct2 > minimized_pose_data_list;

		RNA_LoopCloser rna_loop_closer;

	  AtomTreeMinimizer minimizer;
    float const dummy_tol( 0.00000025);
    bool const use_nblist( true );
    MinimizerOptions options( "dfpmin_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
    options.nblist_auto_update( true );

		if(pose_data_list_.size()==0) return;

		pose = (*pose_data_list_[1].pose_OP);

		//Check scorefxn
		std::cout << "check scorefxn" << std::endl;
		scorefxn_->show( std::cout, pose );

		SilentFileData sfd_dummy; // probably should not be necessary
		sfd_ = new SilentFileData;
		std::string silent_file_for_output = verbose_ ?  "after_minimize_"+silent_file_  :  "" ;

		if( move_map_list_.size()==0 ) move_map_list_=Get_default_movemap( pose );

		if( gap_size == 0 ) {
			std::cout << "five_prime_chain_break_res= " << five_prime_chain_break_res << std::endl;
			//			Add_harmonic_chainbreak_constraint(pose, five_prime_chain_break_res ); //This should be already added during the sampling phrase
		}

		// put in terms to keep chainbreak closed during minimize.
		// even if there is not chainbreak closure, let's have this term (will be 0.0 ),
		//  to keep the number of scoreterms the same in the final outfile.
		scorefxn_->set_weight( linear_chainbreak, 10.0 );


 		for(Size i=1; i<=pose_data_list_.size(); i++){

			std::cout << "Minimizing pose_num= " << i << " of " << pose_data_list_.size() << std::endl;

 			std::string tag=pose_data_list_[i].tag;
			pose=(*pose_data_list_[i].pose_OP); //This actually create a hard copy.....need this to get the graphic working..

			pose.constraint_set( cst_set ); // this is necessary because poses from clustering went thorugh silent struct and lost their constraints

			Remove_virtual_O2Star_hydrogen(pose);
			//			bool const o2star_hydrogens_were_virtual = Remove_virtual_O2Star_hydrogen(pose);
			// if ( o2star_hydrogens_were_virtual )
			o2star_minimize(pose, scorefxn_);

			if(verbose_){
			 	tag[0]='B'; //B for before minimize
				std::cout << tag <<  std::endl;
				//Output_data(sfd_dummy, "before_minimize_"+silent_file_, tag, false, pose, get_native_pose(), actually_moving_res, Is_prepend);
			}

 			///////Minimization/////////////////////////////////////////////////////////////////////////////
			if (gap_size == 0) rna_loop_closer.apply( pose, five_prime_chain_break_res );

			for(Size round=1; round<=move_map_list_.size(); round++){
				if (verbose_) std::cout << "Minimizing round: " << round << " of " << move_map_list_.size() << std::endl;
				core::kinematics::MoveMap mm=move_map_list_[round];

				minimizer.run( pose, mm, *(scorefxn_), options );
				o2star_minimize(pose, scorefxn_);
				if( gap_size == 0 ){
					Real mean_dist_err = rna_loop_closer.apply( pose, five_prime_chain_break_res );
					std::cout << "mean_dist_err (round= " << round << " ) = " <<  mean_dist_err << std::endl;
				}
			}

			////////////////////////////////Final screening //////////////////////////////////////////////
			if ( base_centroid_screener_ && !base_centroid_screener_->Update_base_stub_list_and_Check_that_terminal_res_are_unstacked( pose, true /*reinitialize*/ ) )  continue;

			if ( !Pose_screening(pose, tag, sfd_dummy) ) continue;


			//////////////////////////////////////////////////////////////////////////////////////////////

			pose_data_struct2 pose_data;

			pose_data.tag=tag;
			pose_data.score=(*scorefxn_)(pose);
			pose_data.pose_OP=new pose::Pose;
			(*pose_data.pose_OP)=pose;
			minimized_pose_data_list.push_back(pose_data);

			// might as well output as we go along -- no longer clustering in between.

			tag[ 0 ] = 'M';
			Output_data( *sfd_, silent_file_for_output, tag, false, pose, get_native_pose(), actually_moving_res, Is_prepend);

			std::cout << "Total time in StepWiseRNA_Minimizer: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

  	}


		if(  minimized_pose_data_list.size() == 0 ){ // output something or we'll get weird errors in the workflow...
			Output_data( *sfd_, silent_file_for_output, "DUMMY", false, pose, get_native_pose(), actually_moving_res, Is_prepend);
		}

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
// 			Output_data(sfd_dummy, silent_file_, pose_data.tag, false, pose, get_native_pose(), actually_moving_res, Is_prepend);

// 		}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	//Cannot pass pose in as constant due to the setPoseExtraScores function
	bool
	StepWiseRNA_Minimizer::Pose_screening(pose::Pose & pose, std::string tag, core::io::silent::SilentFileData & silent_file_data) const{

		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::optimization;

		Size const actually_moving_res( job_parameters_->actually_moving_res() ); // used for rmsd calcs
		bool const Is_prepend(  job_parameters_->Is_prepend() );
		Size const gap_size( job_parameters_->gap_size() ); /* If this is zero or one, need to screen or closable chain break */
		Size const five_prime_chain_break_res = job_parameters_->first_chain_break_res();

		bool pass_screen=true;

		tag[0]='S'; //S for screening

		if( gap_size == 1) {
			if( !Check_chain_closable(pose, five_prime_chain_break_res) ){
				tag[0]='F'; //Mark as unclosable loop failure so that will not be use to rebuild next (last) residue
				pass_screen=false;
			}
		}

		if(native_screen_ && !get_native_pose()){	//Before have the (&& !Is_chain_break condition). Parin Dec 21, 2009
			Real rmsd= suite_rmsd(pose, *get_native_pose(), actually_moving_res, Is_prepend);

			if(rmsd>rmsd_cutoff_){ //rmsd_cutoff here is previously define as rmsd_cutoff_of_sampling_phase+1. Parin Dec 21, 2009
				tag[0]='R'; //Mark as RMSD failure
				pass_screen=false;
			}
		}

		if(verbose_){
			std::cout << tag <<  std::endl;
			Output_data(silent_file_data, "screen_"+silent_file_, tag, false, pose, get_native_pose(), actually_moving_res, Is_prepend);
  	}


		return pass_screen;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1 <core::kinematics::MoveMap>
	StepWiseRNA_Minimizer::Get_default_movemap( core::pose::Pose const & pose ) const{

		Size const nres( pose.total_residue() );

		utility::vector1 <core::kinematics::MoveMap> move_map_list;

		core::kinematics::MoveMap mm;
		protocols::swa::Figure_out_moving_residues( mm, pose, job_parameters_->working_fixed_res() );

		// Allow sugar torsions to move again (RD 01/31/2010), now
		// that rotamers have been pre-optimized for ribose closure, and
		// sugar_close is turned back on.
		//		Freeze_sugar_torsions(mm, nres); //Freeze the sugar_torsions!

		if(verbose_) Output_movemap(mm, nres);

		for ( Size n = 1; n <= minimize_rounds_; n++ ) move_map_list.push_back(mm);

		return move_map_list;
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
