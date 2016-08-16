// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project
#include <apps/pilot/hpark/NMsearch.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

// Pose
#include <core/id/AtomID.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Scoring
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/dssp/Dssp.hh>

// Normal Mode
#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>

// LoopHash
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <core/io/silent/SilentStruct.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// MD
#include <protocols/md/CartesianMD.hh>

// etc.
#include <numeric/random/random.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <sys/time.h>

// Silent
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

OPT_1GRP_KEY(Integer, fpd, looplen)
OPT_1GRP_KEY(Integer, fpd, loopcut)
OPT_1GRP_KEY(String, fpd, mode)
OPT_1GRP_KEY(Boolean, fpd, torsion)
OPT_1GRP_KEY(Boolean, fpd, cut_insertion)

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::loophash;

namespace myspace{


///  Evaluator
Evaluator::Evaluator(){}
Evaluator::~Evaluator(){}

Evaluator::Evaluator( pose::Pose const &pose, 
											pose::Pose const &native )
{

	native_ = native;
	resmap_.clear();

  for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
    Size ii_pdb( pose.pdb_info()->number( ii ) );
    for ( Size jj = 1; jj <= native_.total_residue(); ++jj ) {
      Size jj_pdb( native_.pdb_info()->number( jj ) );
      if( ii_pdb == jj_pdb ){
				resmap_[ii] = jj;
				break;
      }
    }
  }
}

//// Scheduler
Scheduler::Scheduler( pose::Pose &pose,
											pose::Pose const &native,
											std::string const mode,
											Size const maxiter )
{

  set_default_parameters();
  maxiter_ = maxiter;

	// just dummy
  scoring::ScoreFunctionCOP sfxn 
    = scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

	if( mode.compare("NMscan") == 0 || mode.compare("NMmix") == 0 ||
			mode.compare("NMeval") == 0 || mode.compare("combine") == 0 ){
		if( option[ fpd::torsion ]() ){
			NMmover_ = new protocols::normalmode::TorsionNormalModeMover( sfxn, "CA", 9.0, "min" );
			// Set max pert as 2.0 Angstrom from initial structure
			set_pertscale( 2.0 );
		} else {
			NMmover_ = new protocols::normalmode::CartesianNormalModeMover( sfxn, "CA", 9.0, "min" );
			// Set max pert as 2.0 Angstrom from initial structure
			// Assuming that Cartesian movement is more unnatural than torsion 
			set_pertscale( 2.0 );
		}
		NMmover_->set_cst_sdev( 0.5 ); // strict
	}

	if( mode.compare("LHfull") == 0 || mode.compare("LHloop") == 0 ||
			mode.compare("combine") == 0 ){
		setup_LHsampler( pose );
	} else if ( mode.compare("LHeval") == 0 ){
		setup_LHsampler( pose, 20 );
	}


  // Evaluator
  evaluator_ = Evaluator( pose, native );
  evaluator_.store_gdt0( pose );
  gdtmax_ = evaluator_.gdt0;
}

void
Scheduler::set_default_parameters()
{
  iter_updated_last_ = 0;
  iter_added_last_ = 0;
  maxmode_ = 10;
  max_scale_ = 2.0; // Apprx. pert to bring 2.0 Ang
  min_scale_ = 0.5;
  gdt_refresh_cut_ = 0.8;
  looplen_ = option[ fpd::looplen ]();
	loopcut_ = option[ fpd::loopcut ]();
}

void
Scheduler::set_pertscale( Real const scale )
{
  pert_scales_.resize( 0 );

  //pert_scales_.push_back( -1.5 ); 
  pert_scales_.push_back( -1.0 ); 
  pert_scales_.push_back( -0.8 ); 
  pert_scales_.push_back( -0.6 );
  pert_scales_.push_back( -0.4 ); 
  pert_scales_.push_back( -0.2 ); 
  pert_scales_.push_back( -0.1 );
  //pert_scales_.push_back( -0.05 );
  //pert_scales_.push_back(  0.05 );
  pert_scales_.push_back(  0.1 ); 
  pert_scales_.push_back(  0.2 ); 
  pert_scales_.push_back(  0.4 );
  pert_scales_.push_back(  0.6 ); 
  pert_scales_.push_back(  0.8 );
  pert_scales_.push_back(  1.0 );
  //pert_scales_.push_back(  1.5 );
  for( Size i = 1; i <= pert_scales_.size(); ++i ) pert_scales_[i] *= scale;
}

/// NormalMode methods

pose::Pose
Scheduler::NM_linesearch( pose::Pose const &pose,
													Real &gdtmin, Real &pert_best,
													std::vector< io::silent::SilentStructOP > &ss )
{
	protocols::moves::MoverOP tocen 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );

	// line search
	gdtmin = 0.0;
	pert_best = 0.0;

	pose::Pose pose_return( pose );
	pose::Pose pose_cen( pose );
	tocen->apply( pose_cen );

	Real gdt_init = evaluator_.apply( pose_return );

	for( Size k = 1; k <= pert_scales_.size(); ++k ){
		// Get copy of starting pose in centroid
		pose::Pose pose_work( pose_cen );
		Real pert_scale = pert_scales_[k];

		NMmover_->set_extrapolate_scale( pert_scale );
		NMmover_->apply( pose_work );
		Real gdtmin_k = evaluator_.apply( pose_work );

		io::silent::SilentStructOP silent =
      io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
 
    silent->fill_struct( pose_work );
		ss.push_back( silent );

		if( gdtmin_k > gdtmin ){
			pose_return = pose_work;
			gdtmin = gdtmin_k;
			pert_best = pert_scale;
		}
	}
  // Convertor
	protocols::moves::MoverOP tofa 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );

	if( pose_return.is_centroid() ) tofa->apply( pose_return );
	return pose_return;
}

void 
Scheduler::control_schedule( pose::Pose const &pose,
														 pose::Pose &pose_ref,
														 Real const gdt_to_ref )
{
	// Normal mode refreshing
	if( gdt_to_ref < gdt_refresh_cut_ ){
		NMmover_->refresh_normalmode();
		std::cout << "GDT to ref: " << gdt_to_ref << ", refreshing NormalMode based on gdt to ref cut..." << std::endl;
		pose_ref = pose;
		iter_updated_last_ = 0;

	} else if (iter_updated_last_ >= 20){
		NMmover_->refresh_normalmode();
		std::cout << "Refreshing NormalMode based on update schedule..." << std::endl;
		pose_ref = pose;
		iter_updated_last_ = 0;
	} else {
		iter_updated_last_ ++;
	}
 
	// Normal mode adding
	if( iter_added_last_ >= 50 ){
		std::cout << "Adding NormalMode based on update schedule..." << std::endl;
		maxmode_ += 10;
		iter_added_last_ = 0;
	} else {
		iter_added_last_ ++;
	}
}

void
Scheduler::run_NM( pose::Pose &pose, std::string const mode )
{
	pose::Pose pose_ref( pose );

  for( Size iter = 1; iter <= maxiter_; ++iter ){
		Real gdtmin, gdt_to_ref, pert_best;
    pose::Pose pose_tmp;

		std::vector< io::silent::SilentStructOP > ss;

		if( mode.compare("mix") == 0) {
			NMmover_->set_random_mode( 10, "probabilistic", 0.3 );
			Real pert_scale = ( (max_scale_ - min_scale_)*RG.uniform() + min_scale_ ); // -5.0 ~ 5.0, abs > 1.0
			if( numeric::random::rg().uniform() < 0.5 ) pert_best *= -1.0;

			NMmover_->set_extrapolate_scale( pert_best );
			NMmover_->apply( pose_tmp );

			Real const gdtmin = evaluator_.apply( pose_tmp );
			Real const gdt_to_ref = evaluator_.apply( pose_ref, pose_tmp );

			accept_and_report( pose_tmp, pose, iter, 0, pert_best, gdtmin, gdt_to_ref, false );

		} else if ( mode.compare("scan") == 0 ) {
			// Randomly pick mode
			Size i_mode = (Size)(maxmode_*RG.uniform())+1;
			NMmover_->set_mode( i_mode );
			// Line search
			pose_tmp = NM_linesearch( pose, gdtmin, pert_best, ss );
			Real const gdt_to_ref = evaluator_.apply( pose_ref, pose_tmp );

			accept_and_report( pose_tmp, pose, iter, i_mode, pert_best, gdtmin, gdt_to_ref, false );

			control_schedule( pose, pose_ref, gdt_to_ref );
		}
  }

}

void
Scheduler::run_NM_evalmodes( pose::Pose &pose, 
														 Size const maxmode )
{
	Real const &gdt0 = evaluator_.gdt0;

	protocols::moves::MoverOP tocen 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
	protocols::moves::MoverOP tofa 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );

 	pose::Pose pose_cen( pose );
	tocen->apply( pose_cen );
	pose::Pose pose_ref( pose_cen );
	
	printf("GDT0: %6.2f\n", gdt0*100.0);

  for( Size i_mode = 1; i_mode <= maxmode; ++i_mode ){
		NMmover_->set_mode( i_mode );
		utility::vector1< Size > accepted( pert_scales_.size()/2 );

		Real gdt_top( 0.0 );

		pose::Pose pose_top( pose_cen );
		for( Size k = 1; k <= pert_scales_.size(); ++k ){
			// Get copy of starting pose in centroid
			pose::Pose pose_work( pose_cen );
			Real pert_scale = pert_scales_[k];
			NMmover_->set_extrapolate_scale( pert_scale );
			NMmover_->apply( pose_work );

			Real gdtmin_k = evaluator_.apply( pose_work );
			Real gdt_to_ref = evaluator_.apply( pose_work, pose_ref );

			// Check acceptance
			Real dgdt = gdtmin_k-gdt0;
			if( dgdt > 0.001 ){
				Size i_acc = k;
				if( k > pert_scales_.size()/2 ) i_acc = pert_scales_.size() - k + 1;
				accepted[i_acc]++;
			}

			printf("Mode/Pert/GDTref/GDTmin/dGDTmin: %2d %5.2f %6.2f %6.2f %6.2f\n",
						 i_mode, pert_scale, gdt_to_ref*100.0, gdtmin_k*100.0, dgdt*100.0);

			if( gdtmin_k > gdt_top ){
				gdt_top = gdtmin_k;
				pose_top = pose_work;
			}

		}
		Real const w_dynamic = NMmover_->get_dynamic_scale();

		printf("Mode %2d, w_dynamic %6.2f, Best dGDT: %6.2f, Acceptance: ",
					 i_mode, w_dynamic, (gdt_top-gdt0)*100.0 );
		for( Size i_acc = 1; i_acc <= accepted.size(); ++i_acc ) printf(" %1d", accepted[i_acc] );
		printf("\n");

		if( gdt_top - gdt0 > 0.03 ){
			std::stringstream pdbname("");
			pdbname << option[ out::prefix ]() << i_mode << ".pdb";
			tofa->apply( pose_top );
			pose_top.dump_pdb( pdbname.str() );
		}
  }
}

////////////////////
// Loop Hash Methods
////////////////////

void
Scheduler::setup_LHsampler( pose::Pose const &pose,
														Size const max_struct )
{
  utility::vector1 < core::Size > loop_sizes;
  loop_sizes.push_back( looplen_ );
	
  LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes, 0, 1 );
  loop_hash_library->load_db();

	if( option[ fpd::cut_insertion ]() ){
		LHsampler_->set_conserve_frame( true );
		LocalInserter_CartesianOP simple_inserter( new LocalInserter_Cartesian() );
		LHsampler_ = new LoopHashSampler::LoopHashSampler( loop_hash_library, simple_inserter );
	} else {
		LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
		LHsampler_ = new LoopHashSampler::LoopHashSampler( loop_hash_library, simple_inserter );
	}

  LHsampler_->set_min_bbrms( 20.0 );
  LHsampler_->set_max_bbrms( 1400.0 );
  LHsampler_->set_min_rms( 0.5 );
  LHsampler_->set_max_rms( 4.0 );
	LHsampler_->set_nonideal( true );

	LHsampler_->set_max_struct( max_struct );
}

pose::Pose
Scheduler::loophash_search( pose::Pose const &pose,
														Real &gdt, Size const ires,
														std::vector< io::silent::SilentStructOP > &ss
														)
{

	chemical::ResidueTypeSetCAP cen_rsdset = 
		chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	protocols::moves::MoverOP tocen 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );

	gdt = 0.0;

	// Pick replacing resid
	LHsampler_->set_start_res( ires );
	LHsampler_->set_stop_res( ires );

	// Get centroid pose
	pose::Pose pose_cen( pose );
	tocen->apply( pose_cen );

	LHsampler_->build_structures( pose_cen, ss );

	// Get best among LoopHash results
	pose::Pose pose_work;
	pose::Pose pose_tmp( pose );
	for( Size i_decoy = 0; i_decoy < ss.size(); ++i_decoy ){
		// Convert decoy to pose_work
		ss[i_decoy]->fill_pose( pose_work, *cen_rsdset );
		Real gdt_work = evaluator_.apply( pose_work );

		if( gdt_work > gdt ){ 
			pose_tmp = pose_work;
			gdt = gdt_work;
		}
	}
	return pose_tmp;
}

Size
Scheduler::pick_hashing_res( pose::Pose const &pose,
														 std::string const mode,
														 Size &Lcount )
{
	Size const nres( pose.total_residue() );
	
	if( mode.compare("random") == 0 ){
		return (Size)(numeric::random::rg().uniform()*(nres-looplen_)) + 1;

	} else if ( mode.compare("looponly") == 0 ){
		Size const maxiter = 10;

		for( int Lcount_cut = looplen_-4; Lcount_cut != loopcut_; --Lcount_cut ){
			for( Size iter = 1; iter <= maxiter; ++iter ){
				Size ires = (Size)(numeric::random::rg().uniform()*(nres-looplen_)) + 1;
				Lcount = 0;
				for( Size shift = 0; shift < looplen_; ++shift ){
					if( pose.secstruct( ires+shift ) == 'L' ) Lcount++;
				}

				// Pick only if the region contains at least 5 coils
				if( Lcount > Lcount_cut )	return ires;

			}
		}
	}
}

void
Scheduler::run_LH( pose::Pose &pose, std::string const mode )
{
	Real gdt;
	pose::Pose pose_ref( pose );

	std::vector< io::silent::SilentStructOP > ss;

  for( Size iter = 1; iter <= maxiter_; ++iter ){
		// Pick replacing resid / Search using LH
		Size Lcount;
		Size ires = pick_hashing_res( pose, mode, Lcount );

		pose::Pose pose_tmp = loophash_search( pose, gdt, ires, ss );
		Real gdt_to_ref = evaluator_.apply( pose_ref, pose_tmp );

		// Accept and report
		accept_and_report( pose_tmp, pose, iter, ires, (Real)(Lcount), gdt, gdt_to_ref, true );
  }
}

void
Scheduler::run_LHeval( pose::Pose &pose )
{
	chemical::ResidueTypeSetCAP cen_rsdset = 
		chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	protocols::moves::MoverOP tocen 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
	protocols::moves::MoverOP tofa 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );

	pose::Pose pose_ref( pose );

	Real const &gdt0 = evaluator_.gdt0;
	printf("GDT0: %6.2f\n", gdt0*100.0);

  for( Size ires = 1; ires <= pose.total_residue()-looplen_+1; ++ires ){
		// Pick replacing resid
		LHsampler_->set_start_res( ires );
		LHsampler_->set_stop_res( ires );

		pose::Pose pose_cen( pose_ref );
		tocen->apply( pose_cen );

		Real gdtchk = evaluator_.apply( pose_cen );

		// Count 2ndary residue num
		Size Hcount( 0 );
		Size Ecount( 0 );
		Size Lcount( 0 );
		for( Size shift = 1; shift <= looplen_; ++shift ){
			if( pose.secstruct( ires+shift ) == 'H' ){
				Hcount++;
			} else if ( pose.secstruct( ires+shift ) == 'E' ){
				Ecount++;
			} else {
				Lcount++;
			}
		}

		// Build
		std::vector< core::io::silent::SilentStructOP > decoys;
		LHsampler_->build_structures( pose_cen, decoys );

		pose::Pose pose_work;
		pose::Pose pose_top( pose_cen );

		Size nacc( 0 );
		Real gdt_top( 0.0 );
		for( Size i_decoy = 0; i_decoy < decoys.size(); ++i_decoy ){
			// Convert decoy to pose_work
			decoys[i_decoy]->fill_pose( pose_work, *cen_rsdset );
			Real gdt_k = evaluator_.apply( pose_work );
			Real gdt_to_ref = evaluator_.apply( pose_ref, pose_work );
			Real dgdt = gdt_k - gdt0;

			printf("Res/i_decoy/GDTref/GDT/dGDT: %2d %3d %6.2f %6.2f %6.2f\n",
						 ires, i_decoy, gdt_to_ref*100.0, gdt_k*100.0, dgdt*100.0);

			if( gdt_k > gdt_top ){
				gdt_top = gdt_k;
				pose_top = pose_work;
			}

			if( dgdt > 0.001 ) nacc++;
		}

		printf("Res %3d, H/E/L %d/%d/%d, Best dGDT: %6.2f, Acceptance: %6.2f\n",
					 ires, Hcount, Ecount, Lcount, 
					 (gdt_top-gdt0)*100.0, nacc*100.0/decoys.size() );

		if( gdt_top - gdt0 > 0.03 ){
			std::stringstream pdbname("");
			pdbname << option[ out::prefix ]() << ires << ".pdb";
			tofa->apply( pose_top );
			pose_top.dump_pdb( pdbname.str() );
		}

  }
}

////////////
/// Combined
////////////

void
Scheduler::run_combine( pose::Pose &pose, 
												bool const report_full )
{
  // Set starting...
	pose::Pose pose_ref( pose );
	Real gdtmin, pert_best;

	Real lh_portion = 0.7;

	io::silent::SilentFileData sfd;
	std::stringstream ssname;
	ssname << option[ out::prefix ]() << ".out";

	// Run
  for( Size iter = 1; iter <= maxiter_; ++iter ){
		pose::Pose pose_tmp;

		std::vector< io::silent::SilentStructOP > ss;

		// Select randomly either LH or NM
		bool run_lh( true );
		if( numeric::random::rg().uniform() > lh_portion ) run_lh = false;

		if( run_lh ){
			Size Lcount;
			//Size ires = pick_hashing_res( pose, "looponly", Lcount );
			Size ires = pick_hashing_res( pose, "random", Lcount );

			pose_tmp = loophash_search( pose, gdtmin, ires, ss );
			Real const gdt_to_ref = evaluator_.apply( pose_ref, pose_tmp );

			accept_and_report( pose_tmp, pose, iter, ires, (Real)(Lcount), gdtmin, gdt_to_ref, true );

			// NM with "scan" method
		} else {

			// Randomly pick mode
			Size i_mode = (Size)(maxmode_*RG.uniform())+1;
			NMmover_->set_mode( i_mode );

			// Line search along mode
			//pose_tmp = NM_linesearch( pose, gdtmin, pert_best );
			pose_tmp = NM_linesearch( pose, gdtmin, pert_best, ss );

			Real const gdt_to_ref = evaluator_.apply( pose_ref, pose_tmp );

			accept_and_report( pose_tmp, pose, iter, i_mode, pert_best, gdtmin, gdt_to_ref, false );

			control_schedule( pose, pose_ref, gdt_to_ref );
		}

		if( report_full ){
			for( Size i_ss = 0; i_ss < ss.size(); ++i_ss ){
				std::stringstream tag;
				tag << "search." << iter << "." << i_ss;
				ss[i_ss]->set_decoy_tag( tag.str() );
				sfd.write_silent_struct( *ss[i_ss], ssname.str() );
			}
		}
  }

}

void 
Scheduler::accept_and_report( pose::Pose &pose_tmp,
															pose::Pose &pose,
															Size const iter,
															Size const arg1,
															Real const arg2,
															Real const gdt,
															Real const gdt_to_ref,
															bool const is_lh,
															Real const pdb_dump_cut )
{

	protocols::moves::MoverOP tofa 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );

	Real const &gdt0 = evaluator_.gdt0;

	Size acceptance( 0 );
	if( gdt > gdtmax_ ) acceptance = 1;

	std::string argstr;
	if( is_lh ){
		argstr = "LH/res1/Nloop";
	} else {
		argstr = "NM/mode/pert_";
	}

	printf("Iter/%s/GDTcurr/GDTMax/dGDT: %4d %3d %4.1f %6.2f %6.2f %6.2f %6.2f\n", 
				 argstr.c_str(), iter, arg1, arg2, gdt_to_ref*100.0,
				 gdt*100.0, gdtmax_*100.0, (gdtmax_-gdt0)*100.0 );
	
	if( gdt > gdtmax_ ){
		if( gdt - gdtmax_ > 0.005 ){
			std::stringstream pdbname("");
			pdbname << option[ out::prefix ]() << iter << ".pdb";
			if( pose_tmp.is_centroid() ) tofa->apply( pose_tmp );
			pose_tmp.dump_pdb( pdbname.str() );
		}
		pose = pose_tmp;
		gdtmax_ = gdt;
	}
}
// class Scheduler

void
Scheduler::run_md( pose::Pose &pose )
{
	// First report the initial GDT
	printf( "GDT0: %8.3f\n", evaluator_.gdt0 );

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb( true );
	movemap->set_chi( true );
	movemap->set_jump( true );
	movemap->set( core::id::THETA, true );
	movemap->set( core::id::D, true );

	core::scoring::ScoreFunctionCOP sfxn
		= core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

	protocols::md::CartesianMD md( pose, sfxn, movemap );
	md.use_rattle( true );
	md.set_nstep( 5000 );
	md.set_reportstep( 10 );
	md.set_premin( 200 );

	if( (*sfxn)[ core::scoring::facts_elec ] > 1e-3 ) md.set_context_update_step( 1 );
	md.set_constraint( 2.0 );

	// Save initial pose
	pose::Pose const pose0( pose );

	// Run
	for( Size i_run = 1; i_run <= 1; ++i_run ){
		std::cout << "MD run: " << i_run << std::endl;
		pose = pose0;
		md.apply( pose );
	}
}

} // namespace myspace


int main( int argc, char *argv [] ){

  NEW_OPT( fpd::looplen, "looplen", 10 );
  NEW_OPT( fpd::loopcut, "loopcut", 3 );
  NEW_OPT( fpd::mode, "mode", "" );
  NEW_OPT( fpd::torsion, "torsion", true );
  NEW_OPT( fpd::cut_insertion, "cut_insertion", false );

  devel::init(argc, argv);

	using namespace myspace;

	// Additional options
	std::string mode = option[ fpd::mode ]();

  core::pose::Pose pose, native;
  core::chemical::ResidueTypeSetCAP rsd_set
    = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  
  core::import_pose::pose_from_file( pose, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file); 
  core::import_pose::pose_from_file( native, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file); 

	Scheduler scheduler( pose, native, mode );

	// feed in 2ndary structure info
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	std::cout << "Running with mode " << mode << std::endl;
	if( mode.compare("NMscan") == 0 ){
		scheduler.run_NM( pose, "scan" );
	} else if ( mode.compare("NMmix") == 0 ){
		scheduler.run_NM( pose, "mix" );
	} else if ( mode.compare("NMeval") == 0 ){
		scheduler.run_NM_evalmodes( pose, 50 ); // N modes to evaluate

	} else if ( mode.compare("LHfull") == 0 ){
		scheduler.run_LH( pose, "random" );
	} else if ( mode.compare("LHloop") == 0 ){
		scheduler.run_LH( pose, "looponly" );

	} else if ( mode.compare("LHeval") == 0 ){
		scheduler.run_LHeval( pose );

	} else if ( mode.compare("combine") == 0 ){
		scheduler.run_combine( pose, true );

	} else if ( mode.compare("MD") == 0 ){
		scheduler.run_md( pose );

	} else {
		std::cout << "Empty or inappropriate mode." << std::endl;
		return 0;
	}
 
	std::stringstream pdbname("");
	pdbname << option[ out::prefix ]() << "final.pdb";
	pose.dump_pdb( pdbname.str() );
  return 0;
}

