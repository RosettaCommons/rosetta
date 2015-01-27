// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/WorkUnit_Sampler.cc
/// @brief
/// @author Hahnbeom Park: Generalized as a "Sampler" from "Loop Hasher"

#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <protocols/mpi_refinement/WorkUnit_Loop.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/constraints/util.hh>

#include <core/fragment/FragmentIO.hh>
#include <protocols/loops/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// for LoopHash
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/LoopHashSampler.hh>

// for FragInsert
#include <protocols/hybridization/CartesianSampler.hh>
#include <core/pose/selection.hh>

// for PartialAbinitio
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>
#include <protocols/hybridization/HybridizeFoldtreeDynamic.hh>
#include <protocols/hybridization/AllResiduesChanged.hh>

// for KICCloser
#include <protocols/comparative_modeling/LoopRelaxMover.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <ObjexxFCL/string.functions.hh>
#include <utility/excn/Exceptions.hh>
#include <utility> //for std::pair

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("WorkUnit_Sampler.LoopType");

/// Currently registered samplers:
// LoopHash, FragInsert, KicCloser, PartialAbinitio

////////////////////////////////////////////
//////// WorkUnit LoopHash
WorkUnit_LoopHash::WorkUnit_LoopHash( core::Size start_ir, core::Size end_ir, 
																			core::Size ssid, 
																			core::Size is_global )
{
	set_defaults();
	set_start(start_ir);
	set_end(end_ir);
	set_ssid( ssid );
	set_global( is_global );
}

void
WorkUnit_LoopHash::set_defaults(){}

// Override
void
WorkUnit_LoopHash::init_from_cmd( const core::Size mpi_rank )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	utility::vector1 < core::Size > loop_sizes = option[ lh::loopsizes]();
	core::Size num_partitions = option[ OptionKeys::wum::n_slaves_per_master](); 
    if( option[ lh::num_partitions].user() )
        num_partitions = option[ lh::num_partitions]();
    core::Size assigned_num = mpi_rank % num_partitions;
	try{
		library_ = protocols::loophash::LoopHashLibraryOP(
  	 new protocols::loophash::LoopHashLibrary( loop_sizes, num_partitions, assigned_num ) );
		// load initial library from disk
		library_->load_mergeddb();
	}
	catch( utility::excn::EXCN_Msg_Exception e ){
		e.show( std::cout );
		e.show( std::cerr );
		throw;
	}
	library_->mem_foot_print();
}

void
WorkUnit_LoopHash::run()
{
  using namespace core::pose;
	using namespace protocols::loops;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;

	if( decoys().size() == 0 ){
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	bool global( true );
	if( get_global() == 0 ) global = false;

	core::io::silent::SilentStructCOP start_struct = decoys().get_struct(0);
	core::pose::Pose pose;
	decoys().get_pose( 0, pose );

	TR << "Executing WorkUnit_LoopHash_Mover on ssid " << start_struct->get_energy("ssid") << std::endl;

	// clear the sotre of structures
	decoys().clear();

	runtime_assert( library_ );

	LocalInserterOP inserter;
	//if( global ){
	inserter = LocalInserterOP( new LocalInserter_SimpleMin() );
	//} else {
	//inserter = LocalInserterOP( new LocalInserter_FixFrame() );
	//}
  LoopHashSampler  lsampler( library_, inserter );

  lsampler.set_start_res( get_start()  );
  lsampler.set_stop_res ( get_end() );
  lsampler.set_min_bbrms( 20.0   );
  lsampler.set_max_bbrms( 1400.0 );
  lsampler.set_min_rms( 0.5 );
  lsampler.set_max_rms( 4.0 );

	if( option[ OptionKeys::lh::bss]() || global ) lsampler.set_nonideal(true);

	// convert pose to centroid pose:
	if( global ){
		if( pose.is_fullatom() ){
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
		}
		core::pose::set_ss_from_phipsi( pose );
	}

	TR.Debug << "Running loophash function: Start: " << get_start() << " End: " << get_end() << std::endl;
	core::Size starttime = time(NULL);

	core::scoring::ScoreFunctionCOP sfxn_loc = get_energy( "talaris2013_cart" );

	std::vector < core::io::silent::SilentStructOP > decoys_tmp;
	//if( global ){
	lsampler.build_structures( pose, decoys_tmp ); //, false );
	//} else {
	//lsampler.build_structures( pose, decoys_tmp, true );
	//}

	// Get loop res: will bring unexpected behavior if using more than 1 loop len
	utility::vector1< core::Size > loopres;
	loopres.resize( 0 );
	core::Size const looplen = library_->hash_sizes()[0]; 
	for( core::Size ires = get_start(); 
			 ires <= std::min( pose.total_residue(), get_end()+looplen ); ++ires ) 
		loopres.push_back(ires);

	// Run minpack before storing
	for( core::Size i_str = 0; i_str < decoys_tmp.size(); ++i_str ){
		core::pose::Pose pose_tmp;
		decoys_tmp[i_str]->fill_pose( pose_tmp );

		// Final boolean option: ramp (4steps) or not?
		ramp_minpack_loop2( pose_tmp, loopres, sfxn_loc, true, false, false, 6.0 );
		store_to_decoys( start_struct, pose_tmp );
	}

	core::Size endtime = time(NULL);
	TR.Debug << "Build " << decoys().size() << " structures in " << endtime - starttime << " s " << std::endl;

}

////////////////////////////////////////////
//////// WorkUnit FragInsert
WorkUnit_FragInsert::WorkUnit_FragInsert( core::Size const nsteps,
																					core::Size const scoretype,
																					core::Size const res1,
																					core::Size const res2,
																					bool const fullatom )
{
	set_defaults();
	set_nsteps( nsteps );
	set_scoretype( scoretype );
	set_res1( res1 );
	set_res2( res2 );
	fullatom_ = fullatom;
}

void
WorkUnit_FragInsert::set_defaults(){}

void
WorkUnit_FragInsert::run()
{
	using namespace core::io::silent;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	if( decoys().size() == 0 ){
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	SilentStructCOP start_struct = decoys().get_struct(0);

	core::pose::Pose pose;
	decoys().get_pose( 0, pose );
	runtime_assert( pose.is_fullatom() );

	// clear the store of structures
	decoys().clear();

	// get number of runs
	core::Size const n( 1 ); //let's just run once

	TR << "Executing WorkUnit_FragInsert on " << get_res1() << " " << get_res2() << std::endl;
	core::Size starttime = time(NULL);

	// Starts here
	// 1. Score function setup
	std::string obj_name = get_scoretype() == 2 ? "goap" : "talaris2013_cart_softrep";

	core::scoring::ScoreFunctionOP sfxn_sampling = get_energy( "talaris2013_cart_softrep" );
	core::scoring::ScoreFunctionOP sfxn_cen      = get_energy( "cen_loop" );
	core::scoring::ScoreFunctionOP sfxn_obj      = get_energy( obj_name );

	// 2. Loop region
	std::set<core::Size> user_pos;
	utility::vector1< core::Size > loopres;
	for( core::Size i = get_res1(); i <= get_res2(); ++i ){
		user_pos.insert( i );
		loopres.push_back( i );
	}

	if( user_pos.size() == 0 ){
		TR << "Empty loop region! nothing to execute." << std::endl;
		return;
	}

	// get options
	core::Real mctemp = option[ lh::mpi_metropolis_temp ]();
	core::Size looplen( get_res2() - get_res1() + 1 );
	core::Size n_frags = option[ frags::n_frags ]();

	// WARNING: loops::frag_size ? frags::frag_size?
	utility::vector1< core::Size > frag_sizes_all = option[ frags::frag_sizes ](); // use this instead which doesn't require addtional cmd line
	utility::vector1< core::Size > frag_sizes;
	for( core::Size i = 1; i <= frag_sizes_all.size(); ++i ){
		if( frag_sizes_all[i] <= 2*(looplen-1) ) frag_sizes.push_back( frag_sizes_all[i] );
	} 

	// 3. Setup sampler
	// dynamic fragment
	protocols::hybridization::CartesianSampler sampler;

	// cmd fragment
	//utility::vector1< core::fragment::FragSetOP > fragsets;
	//protocols::hybridization::CartesianSampler sampler( fragsets );

	sampler.set_userpos( user_pos );
	sampler.set_strategy( "user" );
	sampler.set_rms_cutoff( 3.0 );
	sampler.set_restore_csts( false );
	sampler.set_frag_sizes( frag_sizes );
	if( fullatom_ ){
		sampler.set_fa_scorefunction( sfxn_sampling );
		sampler.set_mc_scorefunction( sfxn_obj );
	} else {
		sampler.set_scorefunction( sfxn_cen );
		sampler.set_mc_scorefunction( sfxn_cen );
	}
	sampler.set_nfrags( n_frags );

	// to be tested
	sampler.set_temp( mctemp ); // temperature, default 2.0
	sampler.set_fullatom( fullatom_ ); //check
	sampler.set_bbmove( true ); //
	sampler.set_ncycles( get_nsteps() );
	sampler.set_nminsteps( 25 );

	// SS biasing for extended
	// TODO: balancing SS - let's think about this later...
	/*
	std::map< core::Size, std::pair< char, core::Real > > bias_frac;
	if( option[ in::file::psipred_ss2 ].user() ){
		std::string ssfile = option[ in::file::psipred_ss2 ]();
		std::map< core::Size, utility::vector1< core::Real > > SSprob = read_ss2( ssfile );

		for( core::Size ires = get_res1(); ires <= get_res2(); ++ires ){
			utility::vector1< core::Real > resprob = SSprob.at(ires);
			core::Real probC = resprob[1];
			core::Real probH = resprob[2];
			core::Real probE = resprob[3];

			core::Real confidence = std::abs( probH - probE );

			core::Real probE_boost = probE + (1.0 - confidence)*probC;
			bias_frac[ires] = std::make_pair( 'E', probE_boost );
		}
	} else {
		for( core::Size ires = get_res1(); ires <= get_res2(); ++ires ){
			bias_frac[ires] = std::make_pair( 'E', 0.33 ); // at least make 1/3 extended by default
		}
	}
	*/

	//sampler.set_bias_fraction( bias_frac );

  protocols::moves::MoverOP tofa 
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );

	// 4. Run!
	for( core::Size i = 1; i <= n; ++i ){
		core::pose::Pose pose_work( pose );
		sampler.apply( pose_work );

		//core::io::silent::SilentStructOP ss = 
		//	core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

		if( pose_work.is_centroid() ){
			tofa->apply( pose_work );

			// reconstruct non-loop region
			for( core::Size ires = 1; ires <= pose_work.total_residue(); ++ires ){
				if( ires >= get_res1() && ires <= get_res2() ) continue;

				for( core::Size j=1; j<=pose.residue_type(ires).natoms(); ++j ){
					core::Vector const &xyz_ref = pose.xyz( core::id::AtomID(j,ires) );
					pose_work.set_xyz( core::id::AtomID(j,ires), xyz_ref );
				}
			}
		}
		ramp_minpack_loop2( pose_work, loopres, sfxn_sampling, true, false, false, 0.0 );
		//store_to_decoys( start_struct, ss, "_"+ObjexxFCL::string_of( i ) );
		store_to_decoys( start_struct, pose_work ); 
	}

	core::Size endtime = time(NULL);
	TR.Debug << "Build " << decoys().size() << " structures in ";
	TR.Debug << endtime - starttime << " s " << std::endl;

} // WorkUnit FragInsert

////////////////////////////////////////////
//////// WorkUnit KicCloser
WorkUnit_KicCloser::WorkUnit_KicCloser( core::Size const nsteps,
																				core::Size const scoretype,
																				core::Size const res1,
																				core::Size const res2,
																				bool const kicclose )
{
	set_defaults();
	set_nsteps( nsteps );
	set_scoretype( scoretype );
	set_res1( res1 );
	set_res2( res2 );
	if( kicclose ){
		set_options( "refine_kic" );
	} else {
		set_options( "no" );
	}
}

void
WorkUnit_KicCloser::set_defaults(){}

void
WorkUnit_KicCloser::run()
{
	using namespace core::io::silent;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	if( decoys().size() == 0 ){
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	SilentStructCOP start_struct = decoys().get_struct(0);

	protocols::moves::MoverOP tofa 
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );

	core::pose::Pose pose;
	decoys().get_pose( 0, pose );
	runtime_assert( pose.is_fullatom() );

	// clear the store of structures
	decoys().clear();

	// get number of runs
	core::Size const n( get_nsteps() );

	core::Size starttime = time(NULL);

	utility::vector1< core::Size > loopres;
	loopres.resize( 0 );
	for( core::Size i = get_res1(); i <= get_res2(); ++i ){
		loopres.push_back( i );
		bool is_res_terminus = pose.residue(i).is_terminus();
		if( is_res_terminus ){
			TR << "Detected terminus within loop region. cannot do this." << std::endl;
			return;
		}
	}

	if( loopres.size() == 0 ){
		TR << "Empty loop region! nothing to execute." << std::endl;
		return;
	}

	// extended loopres
	utility::vector1< core::Size > loopresext( loopres );
	if( get_res1() > 1 ) loopresext.push_back( get_res1() - 1 );
	if( get_res2() < pose.total_residue() ) loopresext.push_back( get_res2() + 1 );

	// Starts here
	// 1. Score function setup
	std::string obj_name = ( get_scoretype() == 2 )? "goap" : "talaris2013";

	core::scoring::ScoreFunctionOP sfxn_sampling = get_energy( "talaris2013_cart" );
	core::scoring::ScoreFunctionOP sfxn_obj = get_energy( obj_name );
	core::scoring::ScoreFunctionOP sfxn_cen = get_energy( "cen_loop" );

	if( sfxn_sampling->get_weight( core::scoring::elec_dens_fast ) > 0.0 )
		TR << "Sampling with elec_dens_fast : " << sfxn_sampling->get_weight( core::scoring::elec_dens_fast ) << std::endl;

	// 2. Setup sampler
	protocols::loops::Loops loops;
	loops.add_loop( get_res1(), get_res2() );
	std::string pert, refine;

	pert = "perturb_kic"; //this is just perturbing w/o closure
	refine = get_options(); //this is "THE CLOSURE"

	comparative_modeling::LoopRelaxMover sampler( pert, "no", refine, "no", loops );
	sampler.scorefxns( sfxn_cen, sfxn_sampling );

	TR << "Executing WorkUnit_KicCloser on " << get_res1() << " " << get_res2() << ", using method ";
	TR << pert << " and " << refine << std::endl;

	// 3. Run!
	core::Real emin( 999999.0 );
	core::pose::Pose pose_min( pose );

	for( core::Size i = 1; i <= n; ++i ){ // try multiple times and pick best sfxn_obj pose
		//continue from pose_min - changed 
		core::pose::Pose pose_work( pose_min );

		sampler.apply( pose_work );

		// could be centroid w/o refine call
		if( pose_work.is_centroid() ){
			tofa->apply( pose_work );
			// reconstruct non-loop region
			for( core::Size ires = 1; ires <= pose_work.total_residue(); ++ires ){
				if( ires >= get_res1() && ires <= get_res2() ) continue;

				for( core::Size j=1; j<=pose.residue_type(ires).natoms(); ++j ){
					core::Vector const &xyz_ref = pose.xyz( core::id::AtomID(j,ires) );
					pose_work.set_xyz( core::id::AtomID(j,ires), xyz_ref );
				}
			}
		}

		core::Real score = sfxn_obj->score( pose_work );
		// reject if too close to the starting
		// this happens if KIC recovers the starting torsions
		core::Real looprmsd = core::scoring::CA_rmsd( pose, pose_work, get_res1(), get_res2() );
		TR << "loop sample " << i << ", irmsd/e/emin: " << std::setw(8) << looprmsd << " ";
		TR << std::setw(8) << score << " " << std::setw(8) << emin;

		if( looprmsd < 1.0 ){
			TR << ", rejected by similarity" << std::endl;
		} else if ( score < emin ){
			emin = score;
			pose_min = pose_work;
			TR << ", accepted" << std::endl;
		} else {
			TR << ", rejected by score" << std::endl;
		}
	}

	// ramping relax (necessary if refine_kic does everything)
	// two rounds: tors-min followed by cart-min
	//ramp_minpack_loop2( pose_min, loopresext, sfxn_sampling, false, true, false, 0.0 );
	ramp_minpack_loop2( pose_min, loopresext, sfxn_sampling, true, false, false, 6.0 ); // only cartmin

	store_to_decoys( start_struct, pose_min );

	core::Size endtime = time(NULL);
	TR.Debug << "Build " << decoys().size() << " structures in " << endtime - starttime << " s " << std::endl;

}

////////////////////////////////////////////
//////// WorkUnit PartialAbinitio
WorkUnit_PartialAbinitio::WorkUnit_PartialAbinitio( core::Size const nsteps,
																										bool const reconstruct
																										)
{
  set_defaults();
  set_nsteps( nsteps );
  if( reconstruct ){
    set_reconstruct( 1 );
  } else {
    set_reconstruct( 0 );
  }
}

void
WorkUnit_PartialAbinitio::set_defaults(){}

core::scoring::ScoreFunctionOP
WorkUnit_PartialAbinitio::setup_score( std::string sfxn_name ) const
{
	core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "score0" );

	if( sfxn_name.compare("score0") == 0 ){
		// stage1_1

	} else if( sfxn_name.compare("score1") == 0 ){
		// stage1_2
		sfxn->set_weight( core::scoring::linear_chainbreak, 0.2 );
		sfxn->set_weight( core::scoring::coordinate_constraint, 0.2 );
		sfxn->set_weight( core::scoring::atom_pair_constraint, 0.1 );
		sfxn->set_weight( core::scoring::vdw, 0.2 );
		sfxn->set_weight( core::scoring::env, 1.0 );
		sfxn->set_weight( core::scoring::pair, 1.0 );
		sfxn->set_weight( core::scoring::hs_pair, 2.0 );
		sfxn->set_weight( core::scoring::ss_pair, 0.6 );
		sfxn->set_weight( core::scoring::sheet, 2.0 );

		//STRAND_STRAND_WEIGHTS 1 11
		core::scoring::methods::EnergyMethodOptions score1_options(sfxn->energy_method_options());
		score1_options.set_strand_strand_weights(1,11);
		sfxn->set_energy_method_options(score1_options);

	} else if( sfxn_name.compare("score2") == 0 ||
						 sfxn_name.compare("score5") == 0) {

  // stage1_3
  sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "score2" );
  sfxn->set_weight( core::scoring::vdw, 0.2 );
  sfxn->set_weight( core::scoring::hs_pair, 2.0 );
  sfxn->set_weight( core::scoring::ss_pair, 0.6 );
  sfxn->set_weight( core::scoring::sheet, 2.0 );
  sfxn->set_weight( core::scoring::linear_chainbreak, 0.5 );
  sfxn->set_weight( core::scoring::atom_pair_constraint, 0.1 );

	if( sfxn_name.compare("score5") == 0 ){
			core::scoring::methods::EnergyMethodOptions score5_options(sfxn->energy_method_options());
			score5_options.set_strand_strand_weights(1,11);
			sfxn->set_energy_method_options(score5_options);
		}

	} else if( sfxn_name.compare("score3") == 0 ){
		// stage1_4
		sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
		sfxn->set_weight( core::scoring::linear_chainbreak, 2.0 );
		sfxn->set_weight( core::scoring::vdw, 0.2 );
		sfxn->set_weight( core::scoring::hs_pair, 2.0 );
		sfxn->set_weight( core::scoring::ss_pair, 2.0 );
		sfxn->set_weight( core::scoring::sheet, 2.0 );
		sfxn->set_weight( core::scoring::rsigma, 2.0 );
		sfxn->set_weight( core::scoring::rama, 0.3 );
		sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	}

	return sfxn;
}

void
WorkUnit_PartialAbinitio::run()
{
  using namespace core::io::silent;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  if( decoys().size() == 0 ){
    TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
    return;
  }
  
  SilentStructCOP start_struct = decoys().get_struct(0);
  
  core::pose::Pose pose0;
  decoys().get_pose( 0, pose0 );
  runtime_assert( pose0.is_fullatom() );

  protocols::moves::MoverOP tofa 
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
  protocols::moves::MoverOP tocen 
    ( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

  // clear the store of structures
  decoys().clear();

  // get number of runs
  core::Size starttime = time(NULL);

  core::pose::Pose pose_work( pose0 );
  tocen->apply( pose_work );

  // Ab inition Setup - simplified from original
  // difference from original (or from hybridize):
  //  - no chunk, just big or small fragment
  //  - no smooth_small frag
  //  - no convergence check in stage3
  //  - no stage1, assuming loop is already at reasonable position
  //  - reduced repeats 2000->1000, considering small space to search

  // 1. loop region / foldtree setup
  utility::vector1< core::Real > residue_weights_3mer( pose0.total_residue(), 0.0 );
  utility::vector1< core::Real > residue_weights_9mer( pose0.total_residue(), 0.0 );

  bool contains_loop( false );

  std::string loopstr( option[ lh::loop_string ]() );
  utility::vector1< utility::vector1< core::Size > > loopregions = loopstring_to_loopregions( loopstr );
  utility::vector1< core::Size > loopres = loopstring_to_loopvector( loopstr );

  protocols::loops::Loops chunk;
  core::Size res2_prv( 0 );
  for( core::Size iloop = 1; iloop <= loopregions.size(); ++iloop ){
    utility::vector1< core::Size > &loopreg = loopregions[iloop];

    core::Size res1 = loopreg[1];
    core::Size res2 = loopreg[loopreg.size()];

    for( core::Size ires = res1; ires <= res2; ++ires ) residue_weights_3mer[ ires ] = 1.0;
    for( core::Size ires = res1; ires <= res2; ++ires ) residue_weights_9mer[ ires ] = 1.0;

    if( !pose0.residue( res1 ).is_terminus() && !pose0.residue( res2 ).is_terminus() )
      contains_loop = true;

    if( res1 > 1 ) chunk.add_loop( res2_prv+1, res1 - 1 );

    // why is this happening - virtual residue at the Cterm??
    // -1 : to avoid weird thing on virtual res at the Cterm...
    if( iloop == loopregions.size() && res2 < pose_work.total_residue()-1 )
			chunk.add_loop( res2 + 1, pose_work.total_residue()-1 );

    res2_prv = res2;
  }

  TR << "initialize foldtree" << std::endl;
  utility::vector1< std::pair< core::Size, core::Size > > pairs;
  std::set< core::Size > std_pos;
  protocols::hybridization::HybridizeFoldtreeDynamic foldtree_mover;
  foldtree_mover.initialize( pose_work, chunk, pairs, std_pos );

  for( core::Size ires = 1; ires <= pose_work.total_residue(); ++ires ){
    if( !loopres.contains( ires ) ) constrain_residue( pose_work, ires, 
																											 loopres //"coordinate"
																											 );
  }

  // 2. fragments / scores setup
  core::fragment::FragSetOP frag3
    = core::fragment::FragmentIO().read_data( option[ in::file::frag3 ]() );
  core::fragment::FragSetOP frag9
    = core::fragment::FragmentIO().read_data( option[ in::file::frag9 ]() );
  utility::vector1< core::fragment::FragSetOP > frag3s; frag3s.push_back( frag3 );
  utility::vector1< core::fragment::FragSetOP > frag9s; frag9s.push_back( frag9 );

	// score function
	core::scoring::ScoreFunctionOP score0 = setup_score( "score0" );
	core::scoring::ScoreFunctionOP score1 = setup_score( "score1" );
	core::scoring::ScoreFunctionOP score2 = setup_score( "score2" );
	core::scoring::ScoreFunctionOP score3 = setup_score( "score3" );
	core::scoring::ScoreFunctionOP score5 = setup_score( "score5" );

	//score1->score( pose_work );
	//score1->show( pose_work );

  // cst setup
  if( option[ constraints::cst_file ].user() ){
    TR << "Add input cst for PartialAbinitio: " << option[ constraints::cst_file ]() << std::endl;
    core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose_work );
  }

  // 3. Mover setup
  // steps
  core::Size stage1_1_cycles, stage1_2_cycles, stage1_3_cycles, stage1_4_cycles;
  if( get_reconstruct() ){
    stage1_1_cycles = 1000*get_nsteps();
    stage1_2_cycles = 1000*get_nsteps();
    stage1_3_cycles = 200*get_nsteps(); // 3 to 4: for closure
    stage1_4_cycles = 200*get_nsteps();
  } else {
    stage1_1_cycles = 0;
    stage1_2_cycles = 0;
    stage1_3_cycles = 1000*get_nsteps();
    stage1_4_cycles = 200*get_nsteps();
  }

  utility::vector1< core::Size > jump_anchors; 
  protocols::moves::RandomMoverOP random_big_frag_mover( new protocols::moves::RandomMover() );
  protocols::moves::RandomMoverOP random_small_frag_mover( new protocols::moves::RandomMover() );
  protocols::hybridization::WeightedFragmentTrialMoverOP top_big_frag_mover
		( new protocols::hybridization::WeightedFragmentTrialMover( frag9s, residue_weights_9mer, jump_anchors, 25 ) );

  protocols::hybridization::WeightedFragmentTrialMoverOP small_frag_mover 
    ( new protocols::hybridization::WeightedFragmentTrialMover( frag3s, residue_weights_3mer ) );

  random_big_frag_mover->add_mover( top_big_frag_mover,    1.0 );
  random_small_frag_mover->add_mover( small_frag_mover,    1.0 );

  core::Real const PHI_EXT( -135.0 );
  core::Real const PSI_EXT(  135.0 );

	// Run from here!
  // stage 1-1: fragment moves to get rid of extended chains, up to 2000 cycles until torsions are replaced
  core::Real temp( 2.0 );

  if( get_reconstruct() > 0 ) {
    TR << "Reconstruction called - starting from fully extended at the loop region defined." << std::endl;

    //core::conformation::Conformation &conf = pose_work.conformation();
    for( core::Size ires = 1; ires <= loopres.size(); ++ires ){
      core::Size resno( loopres[ires] );
      //core::conformation::idealize_position( resno, conf );
      pose_work.set_phi( resno, PHI_EXT );
      pose_work.set_psi( resno, PSI_EXT );
      pose_work.set_omega( resno, 180.0 );
    }

    protocols::hybridization::AllResiduesChanged done( pose_work, residue_weights_9mer, jump_anchors );
    protocols::moves::MonteCarloOP mc1( new protocols::moves::MonteCarlo( pose_work, *score0, temp ) );
    mc1->set_autotemp( false, temp );
    (*score0)(pose_work);
    //bool all_res_changed = false;
    for (core::Size i=1; i<=stage1_1_cycles; ++i) {
      //top_big_frag_mover->apply(pose_work);
      random_big_frag_mover->apply( pose_work );
      (*score0)(pose_work);
      mc1->boltzmann(pose_work, "RandomMover");
      if ( done(pose_work) ) {
				//all_res_changed = true;
				break;
      }
    }
    mc1->recover_low( pose_work );
    mc1->reset( pose_work );
  } else {
    TR << "Refinement called - starting from input structure." << std::endl;
  }

  // stage 1-2: fragment and chunk insertion moves with score1, 2000 cycles with autotemp
  if( get_reconstruct() > 0 ) {
    // ramp chainbreak weight as in KinematicAbinitio
    core::Real const setting( 0.25 / 3 * option[ jumps::increase_chainbreak ] );
    score1->set_weight( core::scoring::linear_chainbreak, setting );

    protocols::moves::MonteCarloOP mc2( new protocols::moves::MonteCarlo( pose_work, *score1, temp ) );
    mc2->set_autotemp( true, temp );
    mc2->set_temperature( temp ); // temperature might have changed due to autotemp..
    mc2->reset( pose_work );
    (*score1)(pose_work);

    protocols::moves::TrialMoverOP stage2_trials 
      ( new protocols::moves::TrialMover( random_big_frag_mover, mc2 ) );
    stage2_trials->keep_stats_type( protocols::moves::accept_reject );

    protocols::moves::RepeatMover( stage2_trials, stage1_2_cycles ).apply( pose_work );

    mc2->show_scores();
    mc2->show_counters();
    mc2->recover_low(pose_work);
    mc2->reset( pose_work );
  }

  //pose_work.dump_pdb( option[ out::prefix ]() + ".stage2.pdb" );
  //core::scoring::EnergyMap emap = pose_work.energies().total_energies();
  //{
  //  core::Real cst = emap[ core::scoring::atom_pair_constraint ];
  //  TR.Debug << "atom_pair_constraint after stage2: " << cst << std::endl;
  //}

  // stage 3: fragment insertion moves with alternating score 2 and 5
  if( get_reconstruct() > 0 ) {
    for (int nmacro=1; nmacro<=10; ++nmacro) {
      //TR << "stage1-3, nmacro " << nmacro << std::endl;
      core::Real progress( 1.0* nmacro/10 );
      core::Real const fact(  progress * 1.0/3  * option[ jumps::increase_chainbreak ]);
      core::Real chbrk_weight_stage_3a = 2.5 * fact;
      core::Real chbrk_weight_stage_3b = 0.5 * fact;
      score2->set_weight( core::scoring::linear_chainbreak, chbrk_weight_stage_3a );
      score5->set_weight( core::scoring::linear_chainbreak, chbrk_weight_stage_3b );

      protocols::moves::MonteCarloOP mc3( new protocols::moves::MonteCarlo( pose_work, *score2, temp ) );
      protocols::moves::TrialMoverOP stage3_trials 
				( new protocols::moves::TrialMover( random_big_frag_mover, mc3 ) );
      stage3_trials->keep_stats_type( moves::accept_reject );

      mc3->recover_low(pose_work);
      if ( numeric::mod( nmacro, 2 ) != 0 || nmacro > 6 ) {
				mc3->score_function( *score2 );
				mc3->set_autotemp( true, temp );
				mc3->set_temperature( temp );
				mc3->reset( pose_work );
				(*score2)( pose_work );
      } else {
				mc3->score_function( *score5 );
				mc3->set_autotemp( true, temp );
				mc3->set_temperature( temp );
				mc3->reset( pose_work );
				(*score5)( pose_work );
      }

      protocols::moves::RepeatMover( stage3_trials, stage1_3_cycles ).apply( pose_work );
      mc3->recover_low(pose_work);
      mc3->reset( pose_work );
    }
  }

  // stage 1-4: refinement by inserting 3-mers
  {
    for (int nmacro=1; nmacro<=3; ++nmacro) {
      // ramp chainbreak weight as in KinematicAbinitio
      core::Real progress( 1.0* nmacro/3 );
      core::Real const setting( ( 1.5*progress+2.5 ) * ( 1.0/3) * option[ jumps::increase_chainbreak ]);
      score3->set_weight( core::scoring::linear_chainbreak, setting );
      score3->set_weight( core::scoring::overlap_chainbreak, progress );

      protocols::moves::MonteCarloOP mc4( new protocols::moves::MonteCarlo( pose_work, *score3, temp ) );
      mc4->set_autotemp( true, temp );
      mc4->set_temperature( temp );
      mc4->reset( pose_work );
      (*score3)( pose_work );

      protocols::moves::TrialMoverOP stage4_trials
				( new protocols::moves::TrialMover( random_small_frag_mover,  mc4 ) );

      protocols::moves::RepeatMover( stage4_trials, stage1_4_cycles ).apply( pose_work );
      mc4->recover_low(pose_work);
      mc4->reset( pose_work );
    }
  }
	
  core::scoring::ScoreFunctionOP sfxn_famin = get_energy( "talaris2013_cart" );

  tofa->apply( pose_work );
  superimpose_to_ref( pose0, pose_work, loopres );

  // copy coordinates
  if( contains_loop ) copy_pose_crd( pose0, pose_work, loopregions );

  core::Size chktime = time(NULL);
  TR.Debug << "Abinitio part: " << chktime - starttime << " s " << std::endl;
  TR.Debug << "after pack: " << sfxn_famin->score( pose_work ) << std::endl;

  // pose/loopres/sfxn/nonideal(true)/ramp(true)/eff(false)/envdist(0.0)
  ramp_minpack_loop2( pose_work, loopres, sfxn_famin, false, true, false, 0.0 );
  ramp_minpack_loop2( pose_work, loopres, sfxn_famin, true, false, false, 6.0 );

  // one more round on cart bonded
  core::scoring::EnergyMap emap = pose_work.energies().total_energies();
  core::Real cartbond = emap[ core::scoring::cart_bonded ];
  if( cartbond > pose_work.total_residue() ){
    ramp_minpack_loop2( pose_work, loopres, sfxn_famin, true, false, false, 6.0 );
  }

  TR.Debug << "after min: " << sfxn_famin->score( pose_work ) << std::endl;

	store_to_decoys( start_struct, pose_work );

  core::Size endtime = time(NULL);
  TR.Debug << "Build " << decoys().size() << " structures in " << endtime - starttime << " s " << std::endl;
}

}
}
