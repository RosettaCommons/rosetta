// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author shilei

#include <iostream>
#include <iomanip>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>

//pose
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

//score
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//job distribution
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/tools/make_vector1.hh>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

//MPI
#include <mpi.h>
#include <utility/mpi_util.hh>

//docking
#include <protocols/docking/DockingHighResLegacy.hh>
#include <protocols/docking/metrics.hh>

//jumps
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/Conformation.hh>

//matrix
#include <numeric/xyzMatrix.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>

//id
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>


//silent
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

//options
#include <basic/options/option_macros.hh>

using namespace core;
using namespace core::scoring;
using namespace std;
using utility::vector1;
using basic::options::option;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;
using std::string;

///work/shilei/rosetta/rosetta_source/src/apps/pilot/wendao/test_bbmc.cc
OPT_1GRP_KEY(Integer,docking_parallel,nstruct)
OPT_1GRP_KEY(Boolean,docking_parallel,rtmin)
OPT_1GRP_KEY(Real,docking_parallel,scaling)
OPT_1GRP_KEY(Real,docking_parallel,trans_pert)
OPT_1GRP_KEY(Real,docking_parallel,rot_pert)
OPT_1GRP_KEY(Real,docking_parallel,min_trans_pert)
OPT_1GRP_KEY(Real,docking_parallel,min_rot_pert)
OPT_1GRP_KEY(String,docking_parallel,outtag)
OPT_1GRP_KEY(Real,docking_parallel,temp)
OPT_1GRP_KEY(Integer,docking_parallel,coor_freq)
OPT_1GRP_KEY(Boolean,docking_parallel,swap)
OPT_1GRP_KEY(Real,docking_parallel,swap_factor)
OPT_1GRP_KEY(Boolean,docking_parallel,restart)


///////////////////////////////////////////////////////////////////////////////
int dump_pose_diff(core::pose::Pose const & pose,core::pose::Pose const & ref_pose, utility::vector1< double > & message)  {
	using core::Size;
	using core::Real;
	using basic::Warning;
	using namespace core::id;
	using namespace core::scoring;
	using core::kinematics::Jump;


	int ii=1;
	//std::cout << "pose.num_jump(): " << pose.num_jump() << std::endl;
	for(int jump = 1, jump_end = pose.num_jump(); jump <= jump_end; ++jump) {
		//rotaiton
		for ( int i = 1; i <= 3; ++i ) {
		      for ( int j = 1; j <= 3; ++j ) {
			 message[ii]=pose.jump(jump).get_rotation()(i,j);
			 ++ii;
			}
		}
		//translation
	    	for ( int i = 1; i <= 3; ++i ) {
    			   message[ii]=pose.jump(jump).get_translation()(i);
			   ++ii;
    		}
	}

	//add changes of torsion angles

  	for(Size rsd = 1, rsd_end = pose.total_residue(); rsd <= rsd_end; ++rsd) {
    		bool const is_jump_residue = pose.fold_tree().is_jump_point(rsd);
		int bb_precision = 6;

		Real bb_tol = 1.0, sc_tol = 1.0;
		for(int i = 0; i < bb_precision; ++i) bb_tol /= 10.0;

		if (!is_jump_residue) {
    		for(Size atom = 1, atom_end = pose.residue(rsd).natoms(); atom <= atom_end; ++atom) {
      			core::id::AtomID aid(atom, rsd);
      			core::id::DOF_ID dof_phi(aid, PHI);
        		int precision = bb_precision;
        		Real tol = bb_tol;

        		Real before_phi = ref_pose.dof(dof_phi);
		        Real after_phi = pose.dof(dof_phi);
        		bool const changed_phi = (std::abs(numeric::nearest_angle_radians(after_phi, before_phi) - before_phi) > tol);
        		if( changed_phi ) {
				message[ii]=rsd;
				++ii;
				message[ii]=atom;
				++ii;
				message[ii]=after_phi;
				++ii;
          			//std::cout.precision( precision );
          			//std::cout  << "Changes: " << rsd << " " << atom << " " << before_phi << " " << after_phi << " " << after_phi-before_phi;
          			}
          			//std::cout << '\n';
    		}// end loop over atoms
		} //non-jumps
  	}// end loop over residues
	return ii-1;
}

///////////////////////////////////////////////////////////////////////////////
void change_docking_pose(core::pose::Pose & to_change_pose, utility::vector1< double > & message, int message_size)  {
        using core::Size;
        using core::Real;
        using basic::Warning;
        using namespace core::id;
        using namespace core::scoring;
        using core::kinematics::Jump;

	//add rb changes
        numeric::xyzMatrix<core::Real> rotation_matrix=numeric::xyzMatrix<Real>::rows(message[1],message[2],message[3],message[4],message[5],message[6],message[7],message[8],message[9]);
        numeric::xyzVector<core::Real> translation_vector=numeric::xyzVector<core::Real>::xyzVector(message[10],message[11],message[12]);
        core::kinematics::Jump tmpJump=to_change_pose.jump(1);
        tmpJump.set_rotation(rotation_matrix);
        tmpJump.set_translation(translation_vector);
        to_change_pose.conformation().set_jump(1,tmpJump);

	Real tmp;

	//torsion angle changes
	for (int ii=13; ii<=message_size; ii=ii+3) {
	to_change_pose.set_dof(core::id::DOF_ID(core::id::AtomID(int(message[ii+1]),int(message[ii])),PHI),message[ii+2]);
	}

}

///////////////////////////////////////////////////////////////////////////////
void run_parallel_docking() {

	int my_rank( 0 ), nprocs( 1 ), tag_( 1 );

	MPI_Status stat_;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);/* get number of processes */

	//read input pose on the every node, rather than communicating
	Pose pose,original_pose,native_pose,previous_pose,best_pose;

	//read in initial structures
	Size input_pdb_size=basic::options::option[ basic::options::OptionKeys::in::file::s ]().size();
	Size random_index=numeric::random::random_range(1,input_pdb_size);
	std::string pdbnamei=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[random_index];
	core::import_pose::pose_from_pdb( pose, pdbnamei.c_str() );
	//std::cout << " node: " << my_rank << " random_index: "  << random_index << std::endl;

	//keep in the same reference poses
	std::string pdbname1=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
	core::import_pose::pose_from_pdb( original_pose, pdbname1.c_str() );

	//read in the native pose
        if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
                core::import_pose::pose_from_pdb( native_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]() );
        } else {
                throw( utility::excn::EXCN_BadInput("native expected for this app") );
        }


	//silent io setup
	core::io::silent::SilentFileData sfd;
	std::string outtag= basic::options::option[basic::options::OptionKeys::docking_parallel::outtag];
	std::string silentfilename=outtag+"highres_dock_" + string_of( my_rank ) + ".silent";
	sfd.set_filename(silentfilename);

	core::Size  total_atoms(0);
	for(Size rsd = 1, rsd_end = pose.total_residue(); rsd <= rsd_end; ++rsd) {
    		for(Size atom = 1, atom_end = pose.residue(rsd).natoms(); atom <= atom_end; ++atom) {
			++total_atoms;
		}
	}

	//initialize docking
	utility::vector1< core::Size > movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	core::scoring::ScoreFunctionOP docking_scorefxn_high_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );

	//add cst to the scoring
	if (basic::options::option[basic::options::OptionKeys::constraints::cst_file].user()) {
		protocols::simple_moves::ConstraintSetMoverOP docking_constraint_ = new protocols::simple_moves::ConstraintSetMover();
		Real cst_weight_=basic::options::option[basic::options::OptionKeys::constraints::cst_weight];
		docking_constraint_->apply(pose);
		docking_scorefxn_high_->set_weight( core::scoring::atom_pair_constraint, cst_weight_ );
	}

	protocols::docking::DockingHighResLegacyOP docking_highres_mover_ = new protocols::docking::DockingHighResLegacy( movable_jumps_, docking_scorefxn_high_, docking_scorefxn_high_);

	//need to hook up the flags for rtmin, unboundrot, ex1, ex2
	Size nstruct=basic::options::option[basic::options::OptionKeys::docking_parallel::nstruct];
	Size curr_struct=1;

	float initial_trans_magnitude=basic::options::option[basic::options::OptionKeys::docking_parallel::trans_pert];
	float initial_rot_magnitude=basic::options::option[basic::options::OptionKeys::docking_parallel::rot_pert];
	float min_trans_magnitude=basic::options::option[basic::options::OptionKeys::docking_parallel::min_trans_pert];
	float min_rot_magnitude=basic::options::option[basic::options::OptionKeys::docking_parallel::min_rot_pert];
	float decay_rate=basic::options::option[basic::options::OptionKeys::docking_parallel::scaling];
	if (basic::options::option[basic::options::OptionKeys::docking_parallel::rtmin]) {
		docking_highres_mover_->set_rt_min(true);
	}
	float temp=basic::options::option[basic::options::OptionKeys::docking_parallel::temp];
	Size coor_freq=basic::options::option[basic::options::OptionKeys::docking_parallel::coor_freq];
	float swap_factor=basic::options::option[basic::options::OptionKeys::docking_parallel::swap_factor];
	Size swap_freq= std::max(int(coor_freq*swap_factor),5);

	//previous_pose;
	previous_pose=pose;
	float current_trans_magnitude=initial_trans_magnitude;
	float current_rot_magnitude=initial_rot_magnitude;
	core::Real best_Isc=999;

	//used to save the difference relative to the first pose
	utility::vector1< core::Real > jump_transformation(12+total_atoms*3);
	int non_zero_size=0;

	utility::vector0< bool >  communicate_all_process_flags(nprocs,false);
	Size num_comm=0;

	do {

	//std::cout << "Run Rosetta on " << my_rank << std::endl;
	if ( curr_struct % coor_freq ==0 ) {
		random_index=numeric::random::random_range(1,input_pdb_size);
		//std::cout << "curr_struct: " << curr_struct << " node: " << my_rank << " random_index: "  << random_index << std::endl;
		pdbnamei=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[random_index];
		core::import_pose::pose_from_pdb( pose, pdbnamei.c_str() );
		previous_pose=pose;
		current_trans_magnitude=initial_trans_magnitude;
		current_rot_magnitude=initial_rot_magnitude;
		best_Isc=999;
	}
	//For every period of structures, re-read the pdbs from input (select at random)

	//set up  perturbation
	docking_highres_mover_->set_trans_magnitude(current_trans_magnitude);
	docking_highres_mover_->set_rot_magnitude(current_rot_magnitude);

	//run docking
	docking_highres_mover_->apply(pose);

	//compute the energy
	core::Real Isc=protocols::docking::calc_interaction_energy(pose,docking_scorefxn_high_,movable_jumps_);

	core::io::silent::BinarySilentStruct ss(pose,(outtag+"rank_"+string_of( my_rank )+"_"+pdbnamei+"_"+string_of(curr_struct)+".pdb"));
	Real rms = core::scoring::CA_rmsd(pose,native_pose);
	ss.add_energy("Isc",Isc);
	ss.add_energy("rms",rms);
	//core::io::silent::BinarySilentStruct ss(pose,(my_rank+"_"+round+".pdb"));
      	sfd.write_silent_struct(ss,silentfilename,false);
	//Write to silent files and compute rms, score etc

	//store the best pose Isc and their difference from the first one
	if (Isc<best_Isc) {
			best_Isc=Isc;
			jump_transformation.assign(12+total_atoms*3,0);
			non_zero_size=dump_pose_diff(pose, original_pose,jump_transformation);
			best_pose=pose;
	}

	//figure out the rb transformation and packer stuff for easily transformation
	//utility::vector1< core::Real > jump_transformation(12+total_atoms*3);
	//int non_zero_size=dump_pose_diff(pose, original_pose,jump_transformation);
	//find out the first non-zero ones and resize it
	//jump_transformation.resize(non_zero_size);

        //for (int ii=1; ii<=jump_transformation.size(); ++ii)
          //      std::cout <<"Rank "<< my_rank << " jump_transformation[" << ii << "]= " << jump_transformation[ii] << std::endl;

	//save pdbs
	//pose.dump_pdb( "dock_from_rank_" + string_of( my_rank ) + ".pdb" );

	//Pose pose_copy;
	//pose_copy=original_pose;
	//change_docking_pose(pose_copy,jump_transformation,non_zero_size);

	//things seem to work fine of converting changes based on a vector
	//pose_copy.dump_pdb( "dock_from_rank_jump_" + string_of( my_rank ) + ".pdb" );

	//compute the rmsd/transformation/packer information
	std::cout << "Finish " << curr_struct<<  " on " << my_rank << " rms " << rms << " Isc " << Isc << " " << current_trans_magnitude << " " << current_rot_magnitude << std::endl;

	//Sync processes in all calculations
	MPI_Barrier(MPI_COMM_WORLD);

	if (basic::options::option[basic::options::OptionKeys::docking_parallel::swap]) {
	//std::cout << "Node: "<< my_rank<< " swap_freq: " << swap_freq << std::endl ;

	if ( curr_struct % swap_freq ==0 ) {
		swap_freq=std::max(int(swap_freq*swap_factor),2);
		//std::cout << "Node: "<< my_rank<< " curr_struct: " << curr_struct << " best_Isc:" << best_Isc << std::endl ;

	int min_indx=0,max_indx=0;

	//reset communicate_all_process_flags if all communication has happened
	if (num_comm==nprocs) {
		     //std::cout << "Node "<< my_rank << " reset flags and all processors have communicated " << num_comm << std::endl ;
		     communicate_all_process_flags.assign(nprocs,false);
	             num_comm=0;
		}
	//send/receive all Isc to rank 0
	if ( my_rank == 0 ) {

		utility::vector0< core::Real > Isc_all_process(nprocs);
		Isc_all_process[my_rank]=best_Isc;
		for ( int ii = 1; ii <= nprocs-1; ++ii ) {
	//		std::cout << "Node 0 receive double: " << " jump_transformation from node " << ii << std::endl ;
			Isc_all_process[ii]=utility::receive_double_from_node(ii);
		}

		//Time to find out the index of the max and mini value
		Real min_value=999,max_value=-999;
		for ( int ii = 0; ii < Isc_all_process.size(); ++ii ) {
			if ( Isc_all_process[ii] > max_value && communicate_all_process_flags[ii]!=true ){
				max_value=Isc_all_process[ii];
				max_indx=ii;
			}

			if ( Isc_all_process[ii] < min_value && communicate_all_process_flags[ii]!=true ){
				min_value=Isc_all_process[ii];
				min_indx=ii;
			}
		}


		//for (int ii=0; ii<=nprocs-1; ++ii)
			//std::cout <<" "<<ii<<" " << Isc_all_process[ii] << std::endl;

		//std::cout <<" min_indx "<< min_indx <<" min_value " << min_value << std::endl;
		//std::cout <<" max_indx "<< max_indx <<" max_value " << max_value << std::endl;

		//Isc_all_process

	} else {
	//	std::cout << "Node " << my_rank << " send double: " << Isc << std::endl ;
		utility::send_double_to_node( 0, best_Isc);
	}

	//MPI_Bcast
	MPI_Bcast(&max_indx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&min_indx, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	//Label the process has been used
	communicate_all_process_flags[min_indx]=true;
	communicate_all_process_flags[max_indx]=true;
	num_comm=num_comm+2;

	if ( my_rank == min_indx ) {
		//Generate the difference in Pose
		utility::vector1< core::Real > tmp_jump_transformation=jump_transformation;
		//utility::vector1< core::Real > tmp_jump_transformation(12+total_atoms*3);
		//int non_zero_size=dump_pose_diff(pose, original_pose,jump_transformation);
		tmp_jump_transformation.resize(non_zero_size);
		utility::send_integer_to_node( max_indx, non_zero_size);
		//std::cout << "Node min " << my_rank << " Size:" << non_zero_size << std::endl ;
		utility::send_doubles_to_node(max_indx,tmp_jump_transformation);
		//pose.dump_pdb( "dock_from_original_" + string_of( my_rank ) + ".pdb" );

        	//for (int ii=1; ii<=jump_transformation.size(); ++ii)
          	  //    std::cout <<"Rank "<< my_rank << " jump_transformation[" << ii << "]= " << jump_transformation[ii] << std::endl;
		//This should send the data to max_indx

	}

	if ( my_rank == max_indx ) {
		non_zero_size= utility::receive_integer_from_node( min_indx);
		//std::cout << "Node max " << my_rank << " Size:" << non_zero_size << std::endl ;
		utility::vector1< core::Real > tmp_jump_transformation=utility::receive_doubles_from_node(min_indx);
        	//for (int ii=1; ii<=jump_transformation.size(); ++ii)
          	  //    std::cout <<"Rank "<< my_rank << " jump_transformation[" << ii << "]= " << jump_transformation[ii] << std::endl;
		//This should receive the data from min_indx

		//pose.dump_pdb( "dock_from_original_" + string_of( my_rank ) + ".pdb" );
		change_docking_pose(pose,tmp_jump_transformation,non_zero_size);
		//pose.dump_pdb( "dock_from_update_" + string_of( my_rank ) + ".pdb" );

		//save the pose in the previous_pose for restart
		previous_pose=pose;

	}

	//Everyone should work on the bestPose
	if ( my_rank != max_indx ) {
		previous_pose=best_pose;
	}

	//Sync processes in all calculations
	MPI_Barrier(MPI_COMM_WORLD);

	//reduce the search radius after each swap step
	current_trans_magnitude=std::max(current_trans_magnitude*decay_rate,min_trans_magnitude);
	current_rot_magnitude=std::max(current_rot_magnitude*decay_rate,min_rot_magnitude);

	}//end of swap_factor
        }//end of swap

	//each run should start from the previous_pose if restart is true
	if (basic::options::option[basic::options::OptionKeys::docking_parallel::restart]) {
			pose=previous_pose;
	}//use the starting pose
	else{ //Add a metropolis rule for using new poses or original poses
		core::Real Isc_now=protocols::docking::calc_interaction_energy(pose,docking_scorefxn_high_,movable_jumps_);
		core::Real Isc_previous=protocols::docking::calc_interaction_energy(previous_pose,docking_scorefxn_high_,movable_jumps_);
		if ( Isc_now-Isc_previous <= 0 || std::exp( -1 * (Isc_now-Isc_previous) / temp ) > numeric::random::rg().uniform() ) {
			//accept and update previous_pose
          		//std::cout <<"Rank "<< my_rank << " accepting currentE: " << Isc_now << " previous E: " << Isc_previous << " ratio: " << std::exp( -1 * (Isc_now-Isc_previous)) << " random " << numeric::random::rg().uniform() << std::endl;
			previous_pose=pose;
		} else {
			//not accept and don't update previous_pose
          		//std::cout <<"Rank "<< my_rank << " NOT accepting currentE: " << Isc_now << " previous E: " << Isc_previous << " ratio: " << std::exp( -1 * (Isc_now-Isc_previous)) <<" random " << numeric::random::rg().uniform() << std::endl;
			pose=previous_pose;
		}
	} //use the current pose

	//change the perturbation radius
	//current_trans_magnitude=std::max(current_trans_magnitude*decay_rate,min_trans_magnitude);
	//current_rot_magnitude=std::max(current_rot_magnitude*decay_rate,min_rot_magnitude);
	curr_struct++;

	} //done compute one structures
	while (curr_struct<=nstruct);

  	MPI::Finalize();

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {

        // initialize option and random number system
	NEW_OPT(docking_parallel::nstruct, "nstruct for each processor", 250);
	NEW_OPT(docking_parallel::rtmin, "perform rtmin now only in the prepack stage", false);
	NEW_OPT(docking_parallel::scaling, "scaling factor for trans_pert and rot_pert for after each step",0.75);
	NEW_OPT(docking_parallel::trans_pert, "random translational pert",3);
	NEW_OPT(docking_parallel::rot_pert, "random rotational pert",8);
	NEW_OPT(docking_parallel::min_trans_pert, "minimum translational pert",0.5);
	NEW_OPT(docking_parallel::min_rot_pert, "minimum rotational pert",1);
	NEW_OPT(docking_parallel::outtag, "output silent tag","");
	NEW_OPT(docking_parallel::temp, "temperature for accepting Isc moves",1);
	NEW_OPT(docking_parallel::coor_freq, "coor_freq read in new random structures from input", 50);
	NEW_OPT(docking_parallel::swap, "whether to swap poses with low Isc to high Isc", true);
	NEW_OPT(docking_parallel::swap_factor, "factor in terms of steps in exchanging pose with low Isc", 0.5 );
	NEW_OPT(docking_parallel::restart, "whether to start from initial pose after each run", true);
	//NEW_OPT come before devel::init

        devel::init( argc, argv );

        run_parallel_docking();
}
