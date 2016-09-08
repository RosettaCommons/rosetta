//#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <stdio.h>

namespace myspace {

using namespace core;

utility::vector1< Size >
get_touched_res( pose::Pose const pose, 
		 utility::vector1< Size > const loopres
		 )
{

  utility::vector1< Size > touched_res = loopres;
  Real const dist_cut( 6.0 );
  
  for( Size i = 1; i <= loopres.size(); ++i ){
    Size const ires( loopres[i] );
    Vector const Cb_crd_i = pose.residue( ires ).nbr_atom_xyz();

    for( Size jres = 1; jres <= pose.size(); ++jres ){
      if( touched_res.contains( jres ) ) continue;

      Vector const Cb_crd_j = pose.residue( jres ).nbr_atom_xyz();

      if( Cb_crd_i.distance( Cb_crd_j ) < dist_cut )
	touched_res.push_back( jres );

    }
  }

  return touched_res;
}

  //inline
kinematics::MoveMapOP
get_movemap( pose::Pose const pose,
	     utility::vector1< Size > loopres,
	     bool const nonideal
	     )
{
  kinematics::MoveMapOP mm = new core::kinematics::MoveMap;

  // Expand loopres into its neighbors
  utility::vector1< Size > touched_residue = get_touched_res( pose, loopres );

  mm->set_bb( false );
  mm->set_chi( false );
  mm->set_jump( true );
  for( Size ires = 1; ires <= touched_residue.size(); ++ires ){
    Size const resno( touched_residue[ ires ] );
    
    if( loopres.contains( ires ) ) mm->set_bb  ( resno, true );
    mm->set_chi ( resno, true );
    if( nonideal ){
      for( Size j=1; j<=pose.residue_type(resno).natoms(); ++j ){
	mm->set( id::DOF_ID(id::AtomID(j,resno), id::THETA ), true );
	mm->set( id::DOF_ID(id::AtomID(j,resno), id::D ), true );
      }
    }
  }
  return mm;
}

inline
kinematics::MoveMapOP
get_movemap( bool const nonideal )
{
  kinematics::MoveMapOP mm = new kinematics::MoveMap;

  mm->set_bb( true );
  mm->set_chi( true );
  mm->set_jump( true );
  if( nonideal ){
    mm->set( id::D, true );
    mm->set( id::THETA, true );
  }
  return mm;
}

  //inline
protocols::simple_moves::PackRotamersMoverOP
setup_packer( pose::Pose const &pose,
	      kinematics::MoveMap const mm,
	      scoring::ScoreFunctionCOP sfxn )
{
  using namespace pack;
  using namespace pack::task;
  using namespace pack::task::operation;

  // Setup TaskOperation from movemap
  TaskFactoryOP local_tf = new TaskFactory();
  local_tf->push_back( new RestrictToRepacking() );
  PreventRepackingOP turn_off_packing = new PreventRepacking();

  for( Size ires = 1; ires <= pose.size(); ++ires ) {
    if (!mm.get_chi(ires) ){
      turn_off_packing->include_residue(ires);
    }
  }
  local_tf->push_back( turn_off_packing );
  
  //Include current rotamer by default
  local_tf->push_back( new IncludeCurrent() );

  protocols::simple_moves::PackRotamersMoverOP packer
    = new protocols::simple_moves::PackRotamersMover( sfxn );

  packer->task_factory(local_tf);

  return packer;
}

inline
void
ramp_repack_min( pose::Pose &pose,
		 kinematics::MoveMap const mm,
		 scoring::ScoreFunctionCOP sfxn,
		 optimization::MinimizerOptionsCOP minoptions,
		 bool const cartesian
		 ){


  // Movemap, Packer: full repack/min
  protocols::simple_moves::PackRotamersMoverOP packer = setup_packer( pose, mm, sfxn );

  packer->apply( pose );
  
  scoring::ScoreFunction sfxn_loc = *sfxn;
  optimization::MinimizerOptions minoptions_loc = *minoptions;

  core::scoring::EnergyMap full_weights = sfxn_loc.weights();

  // Ramp schedule
  float w_ramp[] = { 0.02, 0.25, 0.55, 1.0 };
  int max_iter[] = { 50, 50, 100, 200 };

  for( int i = 0; i< 4; ++i ){
    minoptions_loc.max_iter( (Size)(max_iter[i]) );
    sfxn_loc.set_weight( scoring::fa_rep, 
			 full_weights[ scoring::fa_rep ] * (Real)(w_ramp[i]) );
    
    if( cartesian ){
      optimization::CartesianMinimizer minimizer;
      minimizer.run( pose, mm, sfxn_loc, minoptions_loc );
    } else {
      optimization::AtomTreeMinimizer minimizer;
      minimizer.run( pose, mm, sfxn_loc, minoptions_loc );
    }
  }
}

inline
std::map< Size, Size>
get_resmap( pose::Pose const &pose,
	    pose::Pose const &ref_pose
	    )
{
  std::map< Size, Size > resmap;

  for ( Size ii = 1; ii <= pose.size(); ++ii ) {
    Size ii_pdb( pose.pdb_info()->number( ii ) );
    
    for ( Size jj = 1; jj <= ref_pose.size(); ++jj ) {
      Size jj_pdb( ref_pose.pdb_info()->number( jj ) );
      
      if( ii_pdb == jj_pdb ){
	resmap[ii] = jj;
	break;
      }
    }
    
  }
  return resmap;
}

inline
void
evaluate_and_write( scoring::ScoreFunctionCOP sfxn,
		    pose::Pose const & native_pose,
		    pose::Pose const & pose0,
		    utility::vector1< pose::Pose > const &poses )
{

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  std::map< Size, Size > resmap = get_resmap( pose0, native_pose );

  std::cout << std::setw(4) << "I ";
  std::cout << std::setw(8) << "IRMSD " << std::setw(8) << "RMSD ";
  std::cout << std::setw(8) << "GDTMM " << std::setw(8) << "GDTHA ";
  std::cout << std::endl;

  Real gdtha, gdtmm, rmsd, irmsd;

  gdtha = scoring::gdtha   ( native_pose, pose0, resmap );
  gdtmm = scoring::CA_gdtmm( native_pose, pose0, resmap );
  rmsd  = scoring::CA_rmsd ( native_pose, pose0, resmap );
  printf( "%4d %8.3f %8.3f %8.4f %8.4f\n", 0, 0.0, rmsd, gdtmm, gdtha );

  for( Size i_pose = 1; i_pose <= poses.size(); ++i_pose ){
    pose::Pose pose_work( poses[i_pose] );

    gdtha = scoring::gdtha   ( native_pose, pose_work, resmap );
    gdtmm = scoring::CA_gdtmm( native_pose, pose_work, resmap );
    rmsd  = scoring::CA_rmsd ( native_pose, pose_work, resmap );

    irmsd = scoring::CA_rmsd( pose0, pose_work );
    
    printf( "%4d %8.3f %8.3f %8.4f %8.4f\n", i_pose, irmsd, rmsd, gdtmm, gdtha );

    std::stringstream pdbname("");
    std::string prefix;
    if( option[ out::prefix ].user() ){
      prefix = option[ out::prefix ]();
    } else {
      prefix = "out";
    }
    pdbname << prefix << "_" << i_pose << ".pdb";
    pose_work.dump_pdb( pdbname.str() );
  }
}

} // myspace
