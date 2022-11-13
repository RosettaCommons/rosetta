// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/MixedMonteCarlo.cc
/// @brief A test case for running/testing the MixedMC mover
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)


#include <iostream>
#include <core/io/raw_data/ScoreFileData.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>
#include <numeric/random/random.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <protocols/docking/metrics.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MixedMonteCarlo.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <devel/init.hh>
#include <map>

// // C++ headers
#include <fstream>
#include <string>

static basic::Tracer TR( "MixedMonteCarloApp" );
static core::Size trans(1), rot(2);


utility::vector1_int
get_movable_jumps(
	core::pose::Pose& pose,
	std::string partner_chain_info
){
	using namespace core;
	using namespace protocols;
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	char second_chain = '_';
	Size cutpoint = 0;
	utility::vector1_int movable_jumps_;
	kinematics::FoldTree f( pose.fold_tree() );

	// identify second chain from input partner chains
	for ( Size i = 1; i <= partner_chain_info.length()-1; i++ ) {
		if ( partner_chain_info[i-1] == '_' ) { second_chain = partner_chain_info[i]; }
	}
	// identify cutpoints
	for ( Size i = 2; i <= pose.size(); ++i ) {
		if ( pdb_info->chain(i) == second_chain ) {
			cutpoint = i-1;
			break;
		}
	}

	Size jump_pos1 ( residue_center_of_mass( pose, 1, cutpoint ) );
	Size jump_pos2 ( residue_center_of_mass( pose, cutpoint+1, pose.size() ) );

	// Setup foldtree based on cutpoints and jump points
	f.clear();
	f.simple_tree( pose.size() );
	f.new_jump( jump_pos1, jump_pos2, cutpoint );
	movable_jumps_.clear();
	movable_jumps_.push_back( 1 );

	Size chain_begin(0), chain_end(0);

	// rebuild jumps between chains N-terminal to the docking cutpoint
	chain_end = cutpoint;
	chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
	while ( chain_begin != 1 ) {
		chain_end = chain_begin-1;
		f.new_jump( chain_end, chain_begin, chain_end);
		chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
	}

	// rebuild jumps between chains C-terminal to the docking cutpoint
	chain_begin = cutpoint+1;
	chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	while ( chain_end != pose.size() ) {
		chain_begin = chain_end+1;
		f.new_jump( chain_end, chain_begin, chain_end);
		chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	}

	f.reorder( 1 );
	f.check_fold_tree();
	pose.fold_tree( f );
	f.show( std::cout );
	return movable_jumps_;
}


int main(int argc, char ** argv) {

	devel::init( argc, argv );

	std::cout << "######################################################" << std::endl;
	std::cout << "#######      Mixed Resolution Monte Carlo      #######" << std::endl;
	std::cout << "######################################################" << std::endl;

	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
	utility::file::FileName native_file = basic::options::option[ basic::options::OptionKeys::in::file::native ].value();
	int nstruct = basic::options::option[ basic::options::OptionKeys::out::nstruct ].value();
	core::Real tuning_param = basic::options::option[ basic::options::OptionKeys::in::tuning_param ].value();
	bool dock;
	//bool dock =  basic::options::option[ basic::options::OptionKeys::docking::docking ].value();
	core::Real trans_mag = 2;
	core::Real rot_mag = 6;
	std::string scorefile_name = basic::options::option[ basic::options::OptionKeys::out::file::scorefile ].value();
	std::string partner_info_ = basic::options::option[ basic::options::OptionKeys::docking::partners ]();
	bool bb_move = false;
	std::string outpath = basic::options::option[ basic::options::OptionKeys::out::path::all ]().path();
	std::string lowres_sfx = basic::options::option[ basic::options::OptionKeys::score::lowres_weights ].user() ? basic::options::option[ basic::options::OptionKeys::score::lowres_weights ].value() : "motif_dock_score";

	if ( filenames.size() > 0 ) {
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
	} else {
		std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}

	if ( tuning_param > 0 ) {
		std::cout << "The tuning parameter is " << tuning_param << std::endl;
	} else {
		tuning_param = 0;
		std::cout << "The tuning parameter is " << tuning_param << std::endl;
	}

	if ( nstruct > 0 ) {
		std::cout << "The nstruct is " << nstruct << std::endl;
	} else {
		nstruct = 10;
		std::cout << "The nstruct is " << nstruct << std::endl;
	}

	if ( basic::options::option[ basic::options::OptionKeys::docking::docking ].user() ) {
		dock = true;
	} else {
		dock = false;
	}

	if ( basic::options::option[ basic::options::OptionKeys::docking::dock_pert ].user() ) {
		utility::vector1< core::Real > dock_pert = basic::options::option[ basic::options::OptionKeys::docking::dock_pert ].value();
		core::Real trans_mag = dock_pert[trans];
		core::Real rot_mag = dock_pert[rot];
		std::cout << "Rot mag : " << rot_mag << std::endl;
		std::cout << "Trans mag : " << trans_mag << std::endl;
	}

	if ( basic::options::option[ basic::options::OptionKeys::out::file::scorefile ].user() ) {
		scorefile_name = basic::options::option[ basic::options::OptionKeys::out::file::scorefile ].value();
	}

	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
	core::pose::PoseOP ref_pose = core::import_pose::pose_from_file( native_file );
	utility::vector1_int movable_jumps_ = get_movable_jumps( *mypose, partner_info_ );

	//std::cout << "Movable Jumps : " << movable_jumps_ << std::endl;

	// Set-up the scorefunctions
	core::scoring::ScoreFunctionOP highres_sfxn = core::scoring::get_score_function();
	core::scoring::ScoreFunctionOP lowres_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( lowres_sfx );
	// 2. Create a MMC object
	// 3. Switch the pose to centroid and send it to the MMC

	// Set-up the centroid switch
	protocols::simple_moves::SwitchResidueTypeSetMoverOP to_centroid_ = utility::pointer::make_shared< protocols::simple_moves::SwitchResidueTypeSetMover >( core::chemical::CENTROID );
	core::pose::PoseOP low_pose = mypose->clone();
	to_centroid_->apply( *low_pose );

	core::Real highres_score = highres_sfxn->score( *mypose );
	core::Real lowres_score = lowres_sfxn->score( *low_pose );

	std::cout << "HighRes score: " << highres_score << std::endl;
	std::cout << "LowRes score: " << lowres_score << std::endl;
	protocols::moves::MonteCarloOP mc =
		utility::pointer::make_shared<protocols::moves::MonteCarlo>( *mypose, *highres_sfxn, 100);

	protocols::moves::MixedMonteCarloOP mmc =
		utility::pointer::make_shared<protocols::moves::MixedMonteCarlo>( *low_pose, *mypose, *lowres_sfxn, *highres_sfxn, tuning_param, 100 );

	protocols::moves::PyMOLObserverOP the_observer =
		protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
	the_observer->pymol().apply( *mypose);

	core::kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
	mm.set_jump(true);

	// Set up rigid body translation/rotation steps
	core::Size dock_jump = mypose->num_jump();
	protocols::rigid::RigidBodyPerturbNoCenterMoverOP rb_mover( new protocols::rigid::RigidBodyPerturbNoCenterMover( dock_jump, rot_mag, trans_mag ) );


	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	core::optimization::AtomTreeMinimizer atm;

	core::pose::Pose copy_pose;
	core::pose::Pose copy_lowpose;

	std::string filename = "score_up.sc";
	char const *out_file_name = filename.c_str();
	std::fstream out_file( out_file_name, std::ios::out );
	out_file << "highres_total\tlowres_total\tIsc\tIrms\tfnat\tdescription" << std::endl;

	core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
	repack_task->restrict_to_repacking();
	core::pack::pack_rotamers( *mypose, *highres_sfxn, repack_task );

	for ( int i = 0; i < nstruct; ++i ) {
		std::cout << "Iteration count: " << i << std::endl;
		core::Size num_res = mypose->total_residue();
		core::Real rnd_unif = numeric::random::uniform();
		core::Size rnd_res = static_cast< core::Size > ( rnd_unif * num_res + 1 );
		core::Real orig_phi = mypose->phi( rnd_res );
		core::Real orig_psi = mypose->psi( rnd_res );
		core::Real phi_delta = numeric::random::gaussian();
		core::Real psi_delta = numeric::random::gaussian();
		if ( bb_move ) {
			mypose->set_phi( rnd_res, orig_phi + phi_delta);
			mypose->set_psi( rnd_res, orig_psi + psi_delta);
		}
		std::cout << "Docking ? " << dock << std::endl;
		if ( dock ) { rb_mover->apply( *mypose );}

		core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
		repack_task->restrict_to_repacking();
		core::pack::pack_rotamers( *mypose, *highres_sfxn, repack_task );

		copy_pose = *mypose;
		atm.run( copy_pose, mm, *highres_sfxn, min_opts );
		*mypose = copy_pose;
		copy_lowpose = copy_pose;
		to_centroid_->apply( copy_lowpose );
		mmc->boltzmann(copy_lowpose, *mypose);
		mmc->show_scores();
		std::cout << "High Res Score: " << (*highres_sfxn)(*mypose) << std::endl;
		std::cout << "Low Res Score: " << (*lowres_sfxn)(copy_lowpose) << std::endl;

		the_observer->pymol().apply( *mypose );
		//core::Real high_total( mypose->energies().total_energies().dot( highres_sfxn->weights() ) );
		//core::Real low_total( copy_lowpose.energies().total_energies().dot( lowres_sfxn->weights() ) );

		std::ostringstream tag_id;
		tag_id << std::setfill('0') << std::setw(5) << i+1;
		std::string tag( "P_" + tag_id.str() );

		out_file << std::fixed << std::showpoint << std::setprecision(3) << std::setw(10) << mypose->energies().total_energies().dot( highres_sfxn->weights() ) << "\t" ;
		out_file << copy_lowpose.energies().total_energies().dot( lowres_sfxn->weights() ) << "\t" ;
		out_file << protocols::docking::calc_interaction_energy( *mypose, highres_sfxn, movable_jumps_ ) << "\t";
		out_file << protocols::docking::calc_Irmsd(*mypose, *ref_pose, highres_sfxn, movable_jumps_ ) << "\t" ;
		out_file << protocols::docking::calc_Fnat(*mypose, *ref_pose, highres_sfxn, movable_jumps_ ) << "\t" ;
		out_file << tag << std::endl;

		// Dump the PDBs
		mypose->dump_pdb(outpath + tag + ".pdb"  );

		std::cout << "Next iteration... " << std::endl;
	}
	return 0;
}
