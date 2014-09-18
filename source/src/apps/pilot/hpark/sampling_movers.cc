#include <apps/pilot/hpark/sampling_utils.hh>
#include <apps/pilot/hpark/sampling_movers.hh>

#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/md/CartesianMD.hh>
#include <protocols/simple_moves/CombinePoseMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/LoopHashLibrary.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/Minimizer.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

#include <numeric/random/random.hh>
#include <utility/vector1.hh>

namespace myspace {


using namespace core;

utility::vector1< pose::Pose >
test_NMrelaxer( pose::Pose const pose,
		scoring::ScoreFunctionCOP sfxn,
		optimization::MinimizerOptionsCOP minoptions,
		Size const nstruct,
		bool const cartesian,
		bool const repack
		)
{
  utility::vector1< pose::Pose > poses;

  Size const NMODE( 5 );
  Size const NSCALE( 2 );
  
  // Modes setup
  utility::vector1< Size > modes;
  for( Size i = 1; i <= NMODE; ++i ) modes.push_back( i );
  
  // Run CartNM & TorsionMin
  protocols::normalmode::NormalModeRelaxMoverOP NMMover;
  kinematics::MoveMapOP mm = get_movemap( false ); // allow only torsion move

  if( cartesian ){
    NMMover = 
      new protocols::normalmode::CartesianNormalModeMover( pose, sfxn, mm, "CA", 10.0, "min" );
      
  } else {
    NMMover = 
      new protocols::normalmode::TorsionNormalModeMover( pose, sfxn, mm, "CA", 10.0, "min" );
  }

  NMMover->set_cartesian_minimize( false );
  NMMover->set_minoption( minoptions );
  
  // Iter over different scales
  for(Size i_mode = 1; i_mode <= modes.size(); ++i_mode )
    NMMover->set_mode( modes[i_mode] );
    
  // Positive direction
  for( Size i = 1; i <= NSCALE; ++i ){
    pose::Pose pose_tmp( pose );
    Real const scale = i;
    
    NMMover->set_extrapolate_scale( scale );
    NMMover->apply( pose_tmp );
    poses.push_back( pose_tmp );
  }
  
  // Negative direction
  NMMover->invert_direction();
  for( Size i = 1; i <= NSCALE; ++i ){
    pose::Pose pose_tmp( pose );
    Real const scale = i;
    
    NMMover->set_extrapolate_scale( scale );
    NMMover->apply( pose_tmp );
    poses.push_back( pose_tmp );
  }

  return poses;
}

utility::vector1< pose::Pose >
test_MD( pose::Pose const pose,
	 scoring::ScoreFunctionCOP sfxn,
	 scoring::ScoreFunctionCOP sfxn_obj,
	 optimization::MinimizerOptionsCOP minoptions,
	 Size const nstruct,
	 Size const nper_trj
	 )
{
  utility::vector1< pose::Pose > poses;
  Size const nrun( nstruct/nper_trj );

  kinematics::MoveMapOP mm = get_movemap( true ); //nonideal

  protocols::md::CartesianMD MD( pose, sfxn, mm );
  MD.set_store_trj( true );
  MD.set_scorefxn_obj( sfxn_obj );

  // default
  MD.set_nstep( 50000 );
  MD.set_temperature( 300.0 );

  for( Size irun = 1; irun <= nrun; ++irun ){
    pose::Pose pose_work( pose );
    MD.apply( pose_work );

    utility::vector1< pose::Pose > poses_out = MD.dump_poses( pose );
    for( Size i = 1; i <= poses_out.size(); ++i ){
      pose_work = poses_out[i];
      // Minimize?
      // optimization::CartesianMinimizer minimizer;
      // minimizer.run( pose_work, mm, *sfxn, minoptions );
      poses.push_back( pose_work );
    }
  }

  return poses;
}

utility::vector1< pose::Pose >
test_relax( pose::Pose const pose,
	    scoring::ScoreFunctionCOP sfxn,
	    Size const nstruct
	    )
{

  std::string const relax_script =
    basic::options::option[ basic::options::OptionKeys::relax::script ]();

  protocols::relax::FastRelax relax( sfxn->clone(), relax_script );

  utility::vector1< pose::Pose > poses;
  for( Size istr = 1; istr <= nstruct; ++istr ){
    pose::Pose pose_work( pose );
    relax.apply( pose_work );
    poses.push_back( pose_work );
  }

  return poses;
}

// full Repack & Cartmin
utility::vector1< pose::Pose >
test_recombine( pose::Pose const pose,
		pose::Pose const pose2,
		scoring::ScoreFunctionCOP sfxn,
		optimization::MinimizerOptionsCOP minoptions,
		Size const nstruct,
		bool const repack
		)
{

  chemical::ResidueTypeSetCAP fa_rsdset = 
    chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  // Right now let's just test torsion combiner
  protocols::simple_moves::CombinePoseMover comb( sfxn, pose );

  comb.set_store_silents( true );
  comb.set_max_struct( nstruct );
  comb.set_max_try( nstruct*2 );

  // defaults
  comb.set_minfrac_crossover( 0.1 );
  comb.set_maxfrac_crossover( 0.4 );

  // Run and return decoys
  pose::Pose pose_tmp( pose2 );
  comb.apply( pose_tmp );
  std::vector< io::silent::SilentStructOP > decoys = comb.return_silent();

  assert( decoys.size() == nstruct );

  // Run 
  kinematics::MoveMapOP mm = get_movemap( false ); 
  utility::vector1< pose::Pose > poses;

  for( Size i_decoy = 0; i_decoy < nstruct-1 ; ++i_decoy ){
    decoys[i_decoy]->fill_pose( pose_tmp, *fa_rsdset );
    // Apply minimization
    ramp_repack_min( pose_tmp, *mm, sfxn, minoptions, true );
    poses.push_back( pose_tmp );
  }

  return poses;
}

// by default run torsion minimization
utility::vector1< pose::Pose >
test_bbgauss( pose::Pose const pose,
	      scoring::ScoreFunctionCOP sfxn,
	      optimization::MinimizerOptionsCOP minoptions,
	      Size const nstruct
	      )
{
  utility::vector1< pose::Pose > poses;
  for( Size istr = 1; istr <= nstruct; ++istr ){
    
  }

  return poses;
}

// Repack & Cartmin
utility::vector1< pose::Pose >
test_loophash( pose::Pose const pose,
	       scoring::ScoreFunctionCOP sfxn,
	       optimization::MinimizerOptionsCOP minoptions,
	       Size const nstruct,
	       bool const local
	       )
{

  using namespace protocols::loophash;

  utility::vector1< pose::Pose > poses;

  chemical::ResidueTypeSetCAP cen_rsdset = 
    chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
  chemical::ResidueTypeSetCAP fa_rsdset = 
    chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
 
  //protocols::moves::MoverOP tocen 
  //  = new protocols::simple_moves::SwitchResidueTypeSetMover( chemical::CENTROID );

  // Setup LoopHash
  // just use 10 res
  utility::vector1<Size> loop_sizes;
  loop_sizes.push_back( 10 );

  LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes, 0, 1 );
  loop_hash_library->load_db();

  LocalInserterOP inserter;
  if( local ){
    inserter = new LocalInserter_FixFrame();
  } else {
    inserter = new LocalInserter_SimpleMin();
  }

  LoopHashSamplerOP LHsampler = new LoopHashSampler( loop_hash_library, inserter );
  
  //LHsampler->set_conserve_frame( true );
  LHsampler->set_min_bbrms( 20.0 );
  LHsampler->set_max_bbrms( 1400.0 );
  LHsampler->set_min_rms( 0.5 );
  LHsampler->set_max_rms( 4.0 );
  //LHsampler->set_nonideal( true );
  LHsampler->set_max_struct( 10 ); // try 10 per each window

  std::vector< io::silent::SilentStructOP > decoys;
  Size const nres( pose.total_residue() - 10 );

  // Iter for random windows
  for( Size istr = 1; istr <= nstruct; ++istr ){
    // Randomly pick resno
    Size const ires = (Size)( nres*RG.uniform() );

    LHsampler->set_start_res( ires );
    LHsampler->set_stop_res( ires );

    std::vector< io::silent::SilentStructOP > decoys;

    // Should input be centroid?
    LHsampler->build_structures( pose, decoys );

    kinematics::MoveMapOP mm;
    if( local ){
      utility::vector1< Size > LHres;
      for( Size i = ires; i < ires; ++i ) LHres.push_back( i );
      mm = get_movemap( pose, LHres, true );

    } else {
      mm = get_movemap( false );
    }

    protocols::simple_moves::PackRotamersMoverOP packer = setup_packer( pose, *mm, sfxn );

    // Store
    pose::Pose pose_tmp;
    optimization::CartesianMinimizer minimizer;

    for( Size i_decoy = 0; i_decoy < decoys.size(); ++i_decoy ){
      //decoys[i_decoy]->fill_pose( pose_tmp, *cen_rsdset );
      decoys[i_decoy]->fill_pose( pose_tmp, *fa_rsdset );

      // Apply minimization at restricted region
      packer->apply( pose_tmp );
      minimizer.run( pose_tmp, *mm, *sfxn, *minoptions );
      poses.push_back( pose_tmp );

      if( poses.size() >= nstruct ) return poses;
    }
  }

  return poses;
}

}
