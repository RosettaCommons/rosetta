// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>

//pose
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>

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
#include <utility/io/izstream.hh>
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
#include <protocols/docking/DockingProtocol.hh>
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

#include <protocols/scoring/Interface.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <boost/timer.hpp>
#include <protocols/docking/util.hh>

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
OPT_1GRP_KEY(Integer,docking_parallel,restart_freq)
OPT_1GRP_KEY(String,docking_parallel,outtag)
OPT_1GRP_KEY(Boolean,docking_parallel,use_cst)
OPT_1GRP_KEY(String,docking_parallel,patchdock_file_prefix)
OPT_1GRP_KEY(Integer,docking_parallel,num_master)


void read_patchdock_entry( utility::vector1< utility::vector1 <double> > &saved_transformations_, string patchdock_fname_ )
{
  utility::io::izstream data( patchdock_fname_ );
  if ( !data )
    utility_exit_with_message( "Cannot open patchdock file: " + patchdock_fname_ );

  std::string line;
  bool entries_found( false );
  while ( getline( data, line ) ) {
    using namespace std;

    utility::vector1< double > transformation(6);

    std::istringstream line_stream( line );
    string first_field;
    line_stream >> first_field;

    if( first_field == "#" ) { entries_found = true; continue; }
    if( !entries_found ) continue;
    core::Size const wheres_pipe( line.find_first_of( "|" ) );
    if( wheres_pipe == string::npos ) break;; // no longer reading entries
    core::Size const transformation_begin( line.find_last_of( "||" ) + 2 );
    std::istringstream transData( line.substr( transformation_begin, 10000) );
    double alpha,beta,gamma,x,y,z;
    if( !transData.fail() ) {
    	transData >> alpha >> beta >> gamma >> x >> y >> z;
        transformation[1]=alpha;
        transformation[2]=beta;
        transformation[3]=gamma;
        transformation[4]=x;
        transformation[5]=y;
        transformation[6]=z;
    }
    saved_transformations_.push_back( transformation );
  }
}

void transform_pose( core::pose::Pose & pose, core::Size const chain, utility::vector1<double> const & t )
{
  core::Size const chain_begin( pose.conformation().chain_begin( chain ) );
  core::Size const chain_end( pose.conformation().chain_end( chain ) );

  numeric::xyzMatrix< double > rotation;
  { //compute rotation matrix (taken from RBSMover), but here expecting radian rather than degrees
    double const sa ( std::sin( t[1] ));
    double const ca ( std::cos( t[1] ));
    double const sb ( std::sin( t[2]  ));
    double const cb ( std::cos( t[2]  ));
    double const sg ( std::sin( t[3] ));
    double const cg ( std::cos( t[3] ));
// Adapted from code sent by Dina Schneidman of the Wolfson lab (Tel-Aviv U)
    rotation.xx( cg * cb ); rotation.xy( -sb*sa*cg - sg * ca ); rotation.xz(  -sb*ca*cg + sg * sa );
    rotation.yx( sg * cb ); rotation.yy(  -sb*sa*sg + ca*cg ); rotation.yz( -sb*ca*sg - sa*cg );
    rotation.zx( sb );            rotation.zy( cb*sa );            rotation.zz(  cb*ca );
  }//compute rotation

//rotate each atom around the geometric centre of the chain
  for( core::Size residue=chain_begin; residue<=chain_end; ++residue ) {
    core::Size const atom_begin( 1 );
    core::Size const atom_end( pose.residue( residue ).natoms() );

    numeric::xyzVector< double > localX, localRX;
    for( core::Size atom=atom_begin; atom<=atom_end; ++atom ) {
      id::AtomID const id( atom, residue );

      localX = pose.xyz( id );
      localRX = rotation * localX;
      pose.set_xyz( id, localRX );
    }
  }

//translate
  numeric::xyzVector< Real > translation(t[4],t[5],t[6]);
  for( core::Size residue=chain_begin; residue<=chain_end; ++residue ) {
    core::Size const atom_begin( 1 );
    core::Size const atom_end( pose.residue( residue ).natoms() );

    for( core::Size atom=atom_begin; atom<=atom_end; ++atom ) {
      id::AtomID const id( atom, residue );

      numeric::xyzVector< double > const new_pos( pose.xyz( id ) + translation );
      pose.set_xyz( id, new_pos );
    }
  }
  // detect disulfides
  if ( option[ in::detect_disulf ].user() ?
      option[ in::detect_disulf ]() : // detect_disulf true
      pose.is_fullatom() // detect_disulf default but fa pose
    )
  {
    pose.conformation().detect_disulfides();
  }
}

///////////////////////////////////////////////////////////////////////////////
/*
int dump_pose_diff(core::pose::Pose const & pose,core::pose::Pose const & ref_pose, utility::vector1< double > & message, core::Size &shift)  {
        using core::Size;
        using core::Real;
        using basic::Warning;
        using namespace core::id;
        using namespace core::scoring;
        using core::kinematics::Jump;


        int ii=shift;
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

        for(Size rsd = 1, rsd_end = pose.size(); rsd <= rsd_end; ++rsd) {
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
        to_change_pose.conformation().set_jump_now(1,tmpJump);

        Real tmp;

        //torsion angle changes
        for (int ii=13; ii<=message_size; ii=ii+3) {
        to_change_pose.set_dof(core::id::DOF_ID(core::id::AtomID(int(message[ii+1]),int(message[ii])),PHI),message[ii+2]);
        }

}
*/

///////////////////////////////////////////////////////////////////////////////
void run_parallel_docking() {

	int my_rank( 0 ), nprocs( 1 )/*, tag_( 1 )*/;

	//MPI_Status stat_;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);/* get number of processes */

	//read input pose on the every node, rather than communicating
	Pose pose,native_pose,original_pose;
        std::string pdbname;
        if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
              pdbname=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
              core::import_pose::pose_from_file( pose, pdbname.c_str() , core::import_pose::PDB_file);
              original_pose=pose;
        } else {
              throw( utility::excn::EXCN_BadInput("expected -s for this app") );
        }

	//read in the native pose
        if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
                core::import_pose::pose_from_file( native_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
        } else {
                throw( utility::excn::EXCN_BadInput("native expected for this app") );
        }

        core::Size  total_atoms(0);
        for(Size rsd = 1, rsd_end = pose.size(); rsd <= rsd_end; ++rsd) {
                for(Size atom = 1, atom_end = pose.residue(rsd).natoms(); atom <= atom_end; ++atom) {
                        ++total_atoms;
                }
        }


	//my_rank 0 reads all patchdock orientations
	//send random one over to different rank
	utility::vector1< utility::vector1<double> > saved_transformations_;

	//check if patchdock_file exist
        if (!basic::options::option[basic::options::OptionKeys::docking_parallel::patchdock_file_prefix].user()) {
                throw( utility::excn::EXCN_BadInput("patchdock should be given for this app") );
        }
	std::string filename=basic::options::option[basic::options::OptionKeys::docking_parallel::patchdock_file_prefix]+utility::to_string(my_rank%10);
	read_patchdock_entry(saved_transformations_,filename);

	//read_patchdock_entry(saved_transformations_, basic::options::option[basic::options::OptionKeys::docking_parallel::patchdock_file]);
	Size total=saved_transformations_.size();
	//Size part=total/nprocs;
	Size random_start=numeric::random::random_range(1,total);
	utility::vector1< double > transformation=saved_transformations_[random_start];
	//pose.dump_pdb( "dump_start_"+string_of(my_rank)+".pdb");
	transform_pose( pose, 2, transformation );

       	original_pose=pose;
	//check if in contact, otherwise slide to contact
	//core::Size const rb_move_jump = 1;
	//protocols::scoring::Interface interface( rb_move_jump );
	//interface.calculate( pose );
	//if (interface.interface_nres()==0) {
	//	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	//	core::scoring::ScoreFunctionOP scorefxn_cen( ScoreFunctionFactory::create_score_function("interchain_cen") );
	//	protocols::docking::DockingLowResOP docking_lowres_mover = new protocols::docking::DockingLowRes( scorefxn_cen, rb_move_jump );
	//	docking_lowres_mover->apply(pose);
	//	protocols::forge::methods::restore_residues( original_pose, pose );
	//	core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);
	//}

       	//original_pose=pose;
	//pose.dump_pdb( "dump_start_transform_"+string_of(my_rank)+".pdb");

	//initialize docking
	utility::vector1< core::Size > movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
        core::scoring::ScoreFunctionOP scorefxn_cen( ScoreFunctionFactory::create_score_function("interchain_cen") );
	core::scoring::ScoreFunctionOP docking_scorefxn_high_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	core::scoring::ScoreFunctionOP cst_score_( core::scoring::ScoreFunctionFactory::create_score_function("empty") );

	//add cst to the scoring
	if (basic::options::option[basic::options::OptionKeys::constraints::cst_file].user()) {
		protocols::simple_moves::ConstraintSetMoverOP docking_constraint_ = new protocols::simple_moves::ConstraintSetMover();
		Real cst_weight_=basic::options::option[basic::options::OptionKeys::constraints::cst_weight];
		docking_constraint_->apply(pose);
		if (basic::options::option[basic::options::OptionKeys::docking_parallel::use_cst]) {
			docking_scorefxn_high_->set_weight( core::scoring::atom_pair_constraint, cst_weight_ );
			cst_score_->set_weight( core::scoring::atom_pair_constraint, cst_weight_ );
		} else {
			cst_score_->set_weight( core::scoring::atom_pair_constraint, cst_weight_ );
		}
	}

	//core::kinematics::FoldTree original_ft = pose.fold_tree();
	//protocols::docking::setup_foldtree( pose, utility::to_string("_"), movable_jumps_);
	//protocols::docking::DockMCMProtocolOP docking_highres_mover_ = new protocols::docking::DockMCMProtocol(movable_jumps_,docking_scorefxn_high_,docking_scorefxn_high_);
	//protocols::docking::DockingInitialPerturbationOP perturber = new protocols::docking::DockingInitialPerturbation( movable_jumps_, true);
	protocols::docking::DockingProtocolOP docking_highres_mover_ = new protocols::docking::DockingProtocol( movable_jumps_, false, true, true, scorefxn_cen, docking_scorefxn_high_);

	//need to hook up the flags for rtmin, unboundrot, ex1, ex2
	Size nstruct=basic::options::option[basic::options::OptionKeys::docking_parallel::nstruct];
	Size curr_struct=1;

	//float initial_trans_magnitude=basic::options::option[basic::options::OptionKeys::docking_parallel::trans_pert];
	//float initial_rot_magnitude=basic::options::option[basic::options::OptionKeys::docking_parallel::rot_pert];
	Size restart_freq=basic::options::option[basic::options::OptionKeys::docking_parallel::restart_freq];
	Size num_master=basic::options::option[basic::options::OptionKeys::docking_parallel::num_master];

	utility::vector1< double > results(18);
	if ( my_rank%num_master== 0 ) {
		utility::vector1< utility::vector1<double> > cached_results_;

		//silent io setup
		std::string outtag= basic::options::option[basic::options::OptionKeys::docking_parallel::outtag];
		std::stringstream silentfilestreamname;
		silentfilestreamname << outtag << "highres_dock_" << my_rank << ".sc";
		std::string silentfilename=silentfilestreamname.str();
		std::ofstream o(silentfilename.c_str());
		Size n=0;
		while ( n < nstruct*(nprocs-1) ) {
		//utility::vector1< double > transformation=utility::receive_doubles_from_node(my_rank/num_master);
		for (int i=my_rank+1; i<std::min(int(my_rank+num_master),int(nprocs)); ++i) {
			//std::cout << "receive results from " << i << " n is " << n << std::endl;
			cached_results_.push_back( utility::receive_doubles_from_node(i) );
			n++;
		}
        //	MPI_Barrier(MPI_COMM_WORLD);

		//std::cout << "Supposingly received result and n is " << n << std::endl;

		if (cached_results_.size()>num_master) {
		for(utility::vector1< utility::vector1<double> >::const_iterator cit = cached_results_.begin(); cit!= cached_results_.end(); ++cit){
		//	std::cout << "Results: " << (*cit)[1]<<" "<<(*cit)[2]<<" "<<(*cit)[3]<<" "<<(*cit)[4]<< " "<<(*cit)[5]<< " "<<(*cit)[6] << std::endl;
			o<<(*cit)[1]<<" "<<(*cit)[2]<<" "<<(*cit)[3]<<" "<<(*cit)[4]<< " "<<(*cit)[5]<< " "<<(*cit)[6]<<" "<< (*cit)[7]<<" "<<(*cit)[8]<<" "<<(*cit)[9]<<" "<<(*cit)[10]<< " "<<(*cit)[11]<< " "<<(*cit)[12] << " "<<(*cit)[13]<< " "<<(*cit)[14]<< " "<<(*cit)[15] << " "<<(*cit)[16]<< " "<<(*cit)[17]<< " "<<(*cit)[18]  << std::endl;
		}
			cached_results_.clear();
		}

		}
		o.close();

	} else {
		do {

		boost::timer timer;
		Size total=saved_transformations_.size();
               	Size random_start=numeric::random::random_range(1,total);

		//std::cout << "Run Rosetta on " << my_rank << " processor " << curr_struct << std::endl;
		if ( curr_struct % restart_freq ==0 ) {
                        utility::vector1< double > transformation=saved_transformations_[random_start];
                        transform_pose( pose, 2, transformation );

			original_pose=pose;
                        //protocols::scoring::Interface interface( rb_move_jump );
                        //interface.calculate( pose );
                        //if (interface.interface_nres()==0) {
                        //        core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
                        //        protocols::docking::DockingLowResOP docking_lowres_mover = new protocols::docking::DockingLowRes( scorefxn_cen, rb_move_jump );
                        //        docking_lowres_mover->apply(pose);
                        //        protocols::forge::methods::restore_residues( original_pose, pose );
                         //       core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);
                        //}
			//original_pose=pose;
			//original_ft = pose.fold_tree();
			//protocols::docking::setup_foldtree( pose, utility::to_string("_"), movable_jumps_);
		}

		//run docking
		//if ( basic::options::option[basic::options::OptionKeys::docking_parallel::pert] ) perturber->apply(pose);
		docking_highres_mover_->apply(pose);
		//pose.fold_tree(original_ft);

		//compute the energy
		double CstScore=0.0;
		if (basic::options::option[basic::options::OptionKeys::constraints::cst_file].user()) {
			(*cst_score_)(pose);
			CstScore= pose.energies().total_energies().dot( cst_score_->weights() );
		} else {
			CstScore=-999;
		}
		(*docking_scorefxn_high_)(pose);
		double totalScore= pose.energies().total_energies().dot( docking_scorefxn_high_->weights() );
		double Isc=protocols::docking::calc_interaction_energy(pose,docking_scorefxn_high_,movable_jumps_);
		Real rms = core::scoring::CA_rmsd(pose,native_pose);

		//std::cout << "start fold-tree afterdocking: "<< pose.fold_tree() << std::endl;
		core::kinematics::FoldTree fold_tree_;
		fold_tree_.clear();
		core::Size nres = pose.size();
		core::Size chain2start(pose.conformation().chain_begin( 2 ));
		core::Size chain1end(pose.conformation().chain_end( 1 ));
		//core::Size jumppoint=chain2start;
		fold_tree_.add_edge( 1, chain1end, core::kinematics::Edge::PEPTIDE );
		fold_tree_.add_edge( chain2start, nres, core::kinematics::Edge::PEPTIDE );
		fold_tree_.add_edge( 1, chain2start, 1 );
		fold_tree_.reorder( 1 );
		pose.fold_tree(fold_tree_);
		//std::cout << "clear fold-tree afterdocking: "<< pose.fold_tree() << std::endl;
		//protocols::docking::DockJumps movable_jumps_out(1);
		//protocols::docking::setup_foldtree( pose, utility::to_string("_"), movable_jumps_out);
		//std::cout << "reset fold-tree afterdocking: "<< pose.fold_tree() << std::endl;

		//utility::vector1< double > results(6+12+total_atoms*3);
		//results.assign(6+12+total_atoms*3,0);
		//int non_zero_size=dump_pose_diff(pose, original_pose,results,7);
		//jump_transformation.resize(non_zero_size);

		utility::vector1< double > results(18);
		results[1]=my_rank;
		results[2]=rms;
		results[3]=Isc;
		results[4]=CstScore;
		results[5]=totalScore;
		results[6]=timer.elapsed() ;
		//results[6]=random_start;
		results[7]=pose.jump(1).get_rotation()(1,1);
		results[8]=pose.jump(1).get_rotation()(1,2);
		results[9]=pose.jump(1).get_rotation()(1,3);
		results[10]=pose.jump(1).get_rotation()(2,1);
		results[11]=pose.jump(1).get_rotation()(2,2);
		results[12]=pose.jump(1).get_rotation()(2,3);
		results[13]=pose.jump(1).get_rotation()(3,1);
		results[14]=pose.jump(1).get_rotation()(3,2);
		results[15]=pose.jump(1).get_rotation()(3,3);
		results[16]=pose.jump(1).get_translation()(1);
		results[17]=pose.jump(1).get_translation()(2);
		results[18]=pose.jump(1).get_translation()(3);
		//pose.dump_pdb( "dump_"+string_of(my_rank)+"_"+string_of(curr_struct)+".pdb");
		//send data to writing node
		//std::cout << "Finish " << curr_struct<<  " on " << my_rank << " rms " << rms << " Isc " << Isc << " random_start " << random_start <<  " send to " << my_rank/num_master<< " transformation " << results[7] << " " << results[8] << " " << results[9] << " " << results[10] << " " << results[11] << " " << results[12] << " " << results[13] << " " << results[14] << " " << results[15] << " " << results[16] << " " << results[17] << " " << results[18] << " fold_tree: "<<  pose.fold_tree() << std::endl;
                utility::send_doubles_to_node((my_rank/num_master)*num_master,results);

		//compute the rmsd/results/packer information

		pose=original_pose;

		curr_struct++;

		} //done compute one structures
		while (curr_struct<=nstruct);
	}//end of worker

  	MPI::Finalize();

}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] ) {

        // initialize option and random number system
	NEW_OPT(docking_parallel::nstruct, "nstruct for each processor", 500);
	NEW_OPT(docking_parallel::outtag, "output silent tag","");
	NEW_OPT(docking_parallel::patchdock_file_prefix, "file containing patchdock output","");
	NEW_OPT(docking_parallel::use_cst, "use cst docking",false);
        NEW_OPT(docking_parallel::restart_freq, "restart_freq read in new random structures from input", 100);
        NEW_OPT(docking_parallel::num_master, "number of masters for writing", 100);
	//NEW_OPT come before devel::init

        devel::init( argc, argv );

        run_parallel_docking();
}
