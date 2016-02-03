// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <apps/pilot/hpark/sampling_utils.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/md/CartesianMD.hh>
#include <protocols/simple_moves/CombinePoseMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/Minimizer.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/io/silent/SilentStruct.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <numeric/random/random.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <sys/time.h>

OPT_1GRP_KEY(String, fpd, method)
OPT_1GRP_KEY(Boolean, fpd, coord_cst)
OPT_1GRP_KEY(Integer, fpd, niter)
OPT_1GRP_KEY(Boolean, fpd, repack)
OPT_1GRP_KEY(IntegerVector, fpd, loopres1)

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
  Size const NSCALE( 3 );
  
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
  for(Size i_mode = 1; i_mode <= modes.size(); ++i_mode ){
    NMMover->set_mode( modes[i_mode] );

		std::cout << "Normal mode for mode no " << i_mode << std::endl;

		// Positive direction
		for( Size i = 1; i <= NSCALE; ++i ){
			pose::Pose pose_tmp( pose );
			Real const scale = i;
			
			std::cout << "Positive, Nscale " << i << " " << std::endl;
			
			NMMover->set_extrapolate_scale( scale );
			NMMover->apply( pose_tmp );
			poses.push_back( pose_tmp );
		}
		
		// Negative direction
		NMMover->invert_direction();
		for( Size i = 1; i <= NSCALE; ++i ){
			pose::Pose pose_tmp( pose );
			Real const scale = i;
			
			std::cout << "Negative, Nscale " << i << " " << std::endl;
			
			NMMover->set_extrapolate_scale( scale );
			NMMover->apply( pose_tmp );
			poses.push_back( pose_tmp );
		}
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
  Size const nrun( std::min( 1, (int)(nstruct/nper_trj)) );

  kinematics::MoveMapOP mm = get_movemap( true ); //nonideal

  protocols::md::CartesianMD MD( pose, sfxn, mm );
  MD.set_store_trj( true );
  MD.set_scorefxn_obj( sfxn_obj );

  // default
  MD.set_nstep( 50000 );
  MD.set_temperature( 300.0 );

	// Starting restraint
	MD.set_constraint( 2.0 );

  for( Size irun = 1; irun <= nrun; ++irun ){
		std::cout << "irun: " << irun << std::endl;

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
							Size const nstruct,
							Size const report_step
							)
{
  utility::vector1< pose::Pose > poses;

	Real const kT( 0.56 );
	//Size const nstep_per_trj( 100000 );
	Size const nstep_per_trj( 1000 );
	//Size const ntrj( std::min( 1, (int)(nstruct / (nstep_per_trj/report_step))));
	Size const ntrj = 1;

	//setup SidechainMC mover
	protocols::simple_moves::sidechain_moves::SidechainMCMover scmc;
	pack::task::PackerTaskOP pt = pack::task::TaskFactory::create_packer_task( pose );
	scmc.set_task( pt );
	pt->restrict_to_repacking();

	scmc.init_task( pose );
	scmc.set_ntrials( 100 );
	scmc.set_prob_uniform( 0.0 );
	scmc.set_prob_withinrot( 0.0 );
	scmc.set_prob_random_pert_current( 0.1 ); 
	scmc.set_preserve_detailed_balance( false );
	//scmc.set_accept_according_to_dunbrack( false );
	scmc.set_temperature( kT );
	scmc.set_scorefunction( *sfxn );
	scmc.setup( sfxn );

	// Setup BBgaussian mover
	pose::Pose pose_work( pose );

	kinematics::MoveMapOP mm = get_movemap( false );
	protocols::simple_moves::BBG8T3AMover bbgmover;
	bbgmover.movemap( mm );

	// Run MC
	protocols::moves::MonteCarloOP mc = 
		new protocols::moves::MonteCarlo( pose_work, *sfxn, kT );

	std::cout << "ntrj/nstep per trj? " << ntrj << " " << nstep_per_trj << std::endl;

  for( Size itrj = 1; itrj <= ntrj; ++itrj ){
		// Restart from initial
		pose_work = pose;

		for( Size istep = 1; istep <= nstep_per_trj; ++istep ){
			Real prob = numeric::random::rg().uniform();
			Real proposal_density_ratio( 1.0 );
			std::string movetype;

			//if( istep%100 == 0 ){
			//	std::cout << "i_trj/istep: " << itrj << " " << istep << std::endl;
			//}

			if ( prob > 0.75 ){
				bbgmover.apply( pose_work );
				movetype = bbgmover.type();
				proposal_density_ratio = bbgmover.last_proposal_density_ratio();
			} else {
				scmc.apply( pose_work );
				movetype = bbgmover.type();
			}
			//
			//mc->set_last_accepted_pose( pose_work );

			mc->boltzmann( pose_work, movetype, proposal_density_ratio);

			if( istep%report_step == 0 )
				poses.push_back( pose_work );

			if( poses.size() >= nstruct ) return poses;
		}
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
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  utility::vector1< pose::Pose > poses;

  chemical::ResidueTypeSetCAP cen_rsdset = 
    chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
  chemical::ResidueTypeSetCAP fa_rsdset = 
    chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
 
  protocols::moves::MoverOP tocen 
    = new protocols::simple_moves::SwitchResidueTypeSetMover( chemical::CENTROID );
  protocols::moves::MoverOP tofa
    = new protocols::simple_moves::SwitchResidueTypeSetMover( chemical::FA_STANDARD );

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
  LHsampler->set_max_struct( 10 ); // try 10 per each window
	LHsampler->set_nonideal( true ); // always because starting is not usually idealized

  std::vector< io::silent::SilentStructOP > decoys;
  Size const nres( pose.total_residue() - 10 );

  // Iter for random windows
  for( Size istr = 1; istr <= nstruct; ++istr ){
		Size ires( 0 );

		if( option[ fpd::loopres1 ].user() ){
			utility::vector1< int > loopres1 = option[ fpd::loopres1 ]();
			Size const nloopres( loopres1.size() );
			Size const i = (Size)( nloopres*RG.uniform() );
			ires = loopres1[i];

		} else {
			// Randomly pick resno
			ires = (Size)( nres*RG.uniform() );
		}

		std::cout << "Run LH at: " << ires << " for " << istr << "th structure..." << std::endl;

    LHsampler->set_start_res( ires );
    LHsampler->set_stop_res( ires );

    std::vector< io::silent::SilentStructOP > decoys;

    // Should input be centroid?
		if( local ){
			LHsampler->build_structures( pose, decoys, true );
		} else {
			pose::Pose pose_in( pose );
			tocen->apply( pose_in );
			LHsampler->build_structures( pose_in, decoys, false );
		}

		std::cout << "LH done! generated: " << decoys.size() << std::endl;

    kinematics::MoveMapOP mm;
    if( local ){
      utility::vector1< Size > LHres;
      for( Size i = ires; i <= ires + 10; ++i ) LHres.push_back( i );

      mm = get_movemap( pose, LHres, true );
      //mm = get_movemap( false );

    } else {
      mm = get_movemap( false );
    }

    protocols::simple_moves::PackRotamersMoverOP packer = setup_packer( pose, *mm, sfxn );

    // Store
    pose::Pose pose_tmp;
    optimization::CartesianMinimizer minimizer;

    for( Size i_decoy = 0; i_decoy < decoys.size(); ++i_decoy ){
			if( local ){
				decoys[i_decoy]->fill_pose( pose_tmp, *fa_rsdset );
			} else {
				decoys[i_decoy]->fill_pose( pose_tmp, *cen_rsdset );
				tofa->apply( pose_tmp );
			}

			/*
			std::stringstream pdbname;
			pdbname << i_decoy << ".init.pdb";
			pose_tmp.dump_pdb( pdbname.str() );
			*/

      // Apply minimization at restricted region
      packer->apply( pose_tmp );
      minimizer.run( pose_tmp, *mm, *sfxn, *minoptions );
      poses.push_back( pose_tmp );

			/*
			std::stringstream pdbname2;
			pdbname2 << i_decoy << ".minpack.pdb";
			pose_tmp.dump_pdb( pdbname2.str() );
			*/

      if( poses.size() >= nstruct ) return poses;
    }
  }

  return poses;
}

} // myspace

int main( int argc, char *argv [] ){

	using namespace core;
	using namespace myspace;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	NEW_OPT( fpd::method, "", "" );
	NEW_OPT( fpd::coord_cst, "", false );
	NEW_OPT( fpd::niter, "", 1 );
	NEW_OPT( fpd::repack, "", false );
	utility::vector1< int > emptyv( 0 );
	NEW_OPT( fpd::loopres1, "", emptyv );

	devel::init(argc, argv);

	//protocols::moves::MoverOP tocen 
	//	= new protocols::simple_moves::SwitchResidueTypeSetMover( chemical::CENTROID );
  chemical::ResidueTypeSetCAP rsd_set
		= chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	// Read and store init pose
  pose::Pose pose0, native_pose;
	import_pose::pose_from_file( pose0, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file);

	import_pose::pose_from_file( native_pose, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);

	scoring::ScoreFunctionCOP sfxn 
    = scoring::ScoreFunctionFactory::create_score_function( option[ score::weights] );

	optimization::MinimizerOptionsOP minoption
		= new optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone",
																					0.000001, true, false, false );

	minoption->max_iter( 50 );

	if( !option[ out::nstruct ].user() || option[ out::nstruct ]() == 0 ){
		std::cout << "!out::nstruct > 0" << std::endl;
		return 0;
	}

	std::string const method( option[ fpd::method ]() );
	Size const nstruct( option[ out::nstruct ]() );	
	Size const niter( option[ fpd::niter ]() );
	bool const repack ( option[ fpd::repack ]() ); // not used???

	utility::vector1< pose::Pose > poses;

	for( Size iter = 1; iter <= niter; ++iter ){
		if( method == "cartNM" ){
			poses = test_NMrelaxer( pose0, sfxn, minoption, nstruct, true, repack );

		} else if( method == "torNM" ){
			poses = test_NMrelaxer( pose0, sfxn, minoption, nstruct, false, repack );

		} else if ( method == "MD" ){
			Size const nper_trj( 10 );
			poses = test_MD( pose0, sfxn, sfxn->clone(), minoption, nstruct, nper_trj ); // just use sfxn for obj function

		} else if ( method == "relax" ){ // Cart/Tor is controlled by inputs
			poses = test_relax( pose0, sfxn, nstruct );

		} else if ( method == "recombine" ){
			pose::Pose pose2;
			import_pose::pose_from_file( pose2, *rsd_set, option[ in::file::s ](2) , core::import_pose::PDB_file);
			poses = test_recombine( pose0, pose2, sfxn, minoption, nstruct, repack );

		} else if ( method == "bbgauss" ){
			poses = test_bbgauss( pose0, sfxn, minoption, nstruct, 10 );

		} else if ( method == "localLH" ){
			poses = test_loophash( pose0, sfxn, minoption, nstruct, true );

		} else if ( method == "globalLH" ){
			poses = test_loophash( pose0, sfxn, minoption, nstruct, false );

		} else {
			std::cout << "No matching methods in [torNM/cartNM/MD/relax/bbgauss/loophash]" << std::endl;
			break;
		}

		// write
		evaluate_and_write( sfxn, native_pose, pose0, poses );
	}

	return 0;
}

