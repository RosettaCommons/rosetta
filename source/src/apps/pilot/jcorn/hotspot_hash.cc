// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file hotspot_hash.cc
/// @brief app to find potential hotspot residue placements on the surface of a protein target
///            Method: make a disembodied residue, with option to make residue backbone invisible to the scorefxn (default)
///                    centroid dock residue to the target, then full-atom dock and repack
/// @author Jacob Corn
/// @created May 2008
/// @usage hotspot_hash
/// [-benchmark] <benchmark an existing 2-chain complex>
/// -residue <list of Rosetta Residue name3s || ALL (default)>
/// -run:n_cycles <# stubs to find> (default 10)
/// -hotspot:target <target PDB to hash>
/// -hotspot:target_res <interesting residue on the target> (default none)
/// -hotspot:target_dist <max distance in A from target_res> (default any)
/// -hotspot:angle <Maximum allowed angle between stubCA, target CoM, and stubCB (default off)
/// -hotspot:angle_res <Residue to use for angle calculation from stubCA, <this option>, and stubCB. Used to determine if stub is pointing towards target. 0 uses the default, which is the targets center of mass>
/// -hotspot:length <length of helix as host for hotspot> (default 1)
/// -hotspot:sc_only <Use sidechain only, or also count backbone?> (default T)
/// -hashfile <input existing hash file>
/// -o <output hash file> (default hash.hsh)
/// -scorecut <absolute score or percentage cut> (default off)
/// -score:weights <weights to use for hashing> (default standard, environmental dependence off)
/// -hotspot:envhb <use environmental dependent hbonds? default off)
/// -score:patch <score patch to use for hashing>


#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/scoring/Interface.hh>

#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <core/types.hh>

#include <utility/file/file_sys_util.hh> // file_exists

using basic::T;
using basic::Error;
using basic::Warning;
#include <sstream>
#include <fstream>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>


using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

basic::Tracer TR( "pilot_apps.jcorn.hotspot_hash");

// routine to benchmark score existing contacts as if they were found by hashing
void benchmark_contacts ( pose::Pose const & start_pose, scoring::ScoreFunctionOP scorefxn )
{
	if ( start_pose.conformation().num_chains() < 2 ) {
		TR << "Pose must contain at least two chains!" << std::endl;
		return;
	}
	//TR << "orig # " << start_pose.pdb_numbering() << std::endl;

	std::string name = start_pose.pdb_info()->name();

	// make a copy of the original pose that we can work with
	core::pose::Pose pose = start_pose;

	// minimize sidechains only to get rid of any strain
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->clear();
	mm->set_chi( true );
	protocols::simple_moves::MinMover min_mover( mm, scorefxn, "dfpmin_armijo_nonmonotone", 1e-5, true, false  );
	min_mover.apply( pose );

	// split the pose into separate chains
	utility::vector1< pose::PoseOP > singlechain_poses;
	singlechain_poses = pose.split_by_chain();

	// calculate all interface residues over all jumps
	utility::vector1< Size > interface( pose.total_residue(), false );
	for( Size i=1; i<=pose.num_jump(); ++i ) {
		protocols::scoring::Interface iface( i );
		iface.distance( 10 );
		iface.calculate( pose );
		for( Size n=1; n<=pose.total_residue(); ++n ) {
			if( iface.is_interface( n ) ) {
				interface[n] = true;
			}
		}
	}
	// scoring header
	TR << "SCORE id chain res seqpos fa_atr fa_rep hbond_bb_sc hbond_sc fa_sol total\n";

	// iterate over all single-chain poses
	for ( utility::vector1<pose::PoseOP>::iterator target_chain_it = singlechain_poses.begin();
		  target_chain_it != singlechain_poses.end();
		  ++target_chain_it ) {

		  //for ( conformation::ResidueOPs::iterator res_it = pose.res_begin();
			//res_it != pose.res_end();
			//	++res_it ) {
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

			// make a convenience residue to avoid constant double-dereferencing iterator of OP
			conformation::Residue const & residue =  pose.residue( ii );
			if( !interface[ residue.seqpos() ] ) continue; // skip if we're not in the interface

			// avoid self-interactions
			//if ( res_chain_it == target_chain_it )  continue;

			// make a new residue for each residue in chain1
			// we do this for ALL residues (not just interface ones) to avoid rb_jump confusion in poses with >2 chains
			//for ( conformation::ResidueOPs::iterator res_it = pose->res_begin();
			//	  res_it != (*res_chain_it)->res_end();
			//	  ++res_it ) {

			// don't operate on G or C (since we don't consider them as hotspots, anyway)
			// in the future, can get a HUGE speedup by checking if the residue is at the interface before scoring
			// for now, just working on quick test
			if ( ( residue.name3() == "CYS" ) || ( residue.name3() == "GLY" ) ) {
				continue;
			}

			// append the new residue and set its backbone virtual (emulate de novo hashing)
			(*target_chain_it)->append_residue_by_jump( residue, (*target_chain_it)->total_residue(), "", "", true );
			if( option[hotspot::sc_only]() ) core::pose::add_variant_type_to_pose_residue( **target_chain_it, "SHOVE_BB", (*target_chain_it)->total_residue() );

			//Size const pdb_seqpos = residue->seqpos() + pose.conformation().chain_begin( chain2 ) - 1;
			Size const placed_seqpos = (*target_chain_it)->total_residue();
			Size const pdb_seqpos = residue.seqpos();

			// build up scores
			(*scorefxn)( **target_chain_it );
			// I think accumulate_residue_total_energies is now obsolete
			//scorefxn->accumulate_residue_total_energies( *target_chain_it );
			Real weighted_fa_atr = (*target_chain_it)->energies().residue_total_energies( placed_seqpos )[scoring::fa_atr] * scorefxn->weights()[scoring::fa_atr];
			Real weighted_fa_rep = (*target_chain_it)->energies().residue_total_energies( placed_seqpos )[scoring::fa_rep] * scorefxn->weights()[scoring::fa_rep];
			Real weighted_hbond_bb_sc = (*target_chain_it)->energies().residue_total_energies( placed_seqpos )[scoring::hbond_bb_sc] * scorefxn->weights()[scoring::hbond_bb_sc];
			Real weighted_hbond_sc = (*target_chain_it)->energies().residue_total_energies( placed_seqpos )[scoring::hbond_sc] * scorefxn->weights()[scoring::hbond_sc];
			Real weighted_fa_sol = (*target_chain_it)->energies().residue_total_energies( placed_seqpos )[scoring::fa_sol] * scorefxn->weights()[scoring::fa_sol];
			Real weighted_contact_score = weighted_fa_atr + weighted_fa_rep + weighted_hbond_bb_sc + weighted_hbond_sc + weighted_fa_sol;

			// weighted score output
			// fa_atr fa_rep hbond_bb_sc hbond_sc fa_sol total_score
			if ( weighted_contact_score < -0.001 ) {
				TR  << "WTD_SCORE " << name.substr(0,4) << " " << pose.pdb_info()->chain(pdb_seqpos) << " " << residue.name() << " " <<  pose.pdb_info()->number(pdb_seqpos) << " "
					<< weighted_fa_atr << " "
					<< weighted_fa_rep << " "
					<< weighted_hbond_bb_sc << " "
					<< weighted_hbond_sc << " "
					<< weighted_fa_sol << " "
					<<  weighted_contact_score << "\n";
			}

			// clear the new residue so that we can start over
			(*target_chain_it)->conformation().delete_residue_slow( placed_seqpos );
		} // end residue iteration
	} // end target chain iteration
//	} // end residue source chain iteration
	// flush score output
	TR.flush();
} // end benchmark_contacts


int
main( int argc, char * argv [] )
{
	try {

	using namespace scoring;

	OPT(hotspot::target);
	OPT(hotspot::residue);
	OPT(hotspot::benchmark);
	OPT(hotspot::hashfile);
	OPT(hotspot::target_res);
	OPT(hotspot::target_dist);
	OPT(hotspot::threshold);
	OPT(hotspot::sc_only);
	OPT(out::scorecut);
	OPT(run::n_cycles);
	OPT(out::file::o);
	OPT(out::overwrite);
	OPT(score::weights);
	OPT(score::patch);
	OPT(hotspot::envhb);
	OPT(hotspot::angle);
	OPT(hotspot::angle_res);


	devel::init(argc, argv);

	// turn on extra rotamers for finding good hotspots
	option[ packing::ex1::ex1 ].value( true );
	option[ packing::ex2::ex2 ].value( true );
	option[ packing::extrachi_cutoff ].value( 0 );
	//option[ docking::fake_native ].value( true );  THIS OPTION IS DEPRECATED
	option[ docking::randomize1 ].value( true );
	option[ docking::randomize2 ].value( true );
	option[ in::file::sucker_params ].value( "scoring/sucker/pusher.params" );

	// necessary to make sure NANs in hbonding don't cause an exit
	option[ in::file::fail_on_bad_hbond ].value( false );

	// BENCHMARKING EXISTING COMPLEXES
	if ( option[ hotspot::benchmark ].user() ) {

		core::scoring::ScoreFunctionOP scorefxn;
		if( option[ score::weights ].user() ) {
			scorefxn = core::scoring::getScoreFunction();
		} else {
			scorefxn = core::scoring::getScoreFunctionLegacy( "score13" );
			scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
			scorefxn->set_weight( core::scoring::envsmooth, 0 );
		}
		core::scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		TR << "Setting envinromental hbonds to: " << option[ hotspot::envhb]() << std::endl;
		options.hbond_options().use_hb_env_dep( option[ hotspot::envhb]() );
		scorefxn->set_energy_method_options( options );

		// do benchmarking
		pose::Pose pose;
		if ( option[in::file::l].user() ) {
			utility::vector1< std::string > pdbnames = basic::options::start_files() ;
			for ( utility::vector1<std::string>::iterator filename( pdbnames.begin() );
				  filename != pdbnames.end(); ++filename ) {
				core::import_pose::pose_from_pdb( pose, *filename );
				benchmark_contacts( pose, scorefxn );
			}
		}
		else if ( option[in::file::s].user() ) {
			core::import_pose::pose_from_pdb( pose, basic::options::start_file() );
			benchmark_contacts( pose, scorefxn );

		} else {
			utility_exit_with_message( "No files given: Use either -s or -l to designate a single pdb or a list of pdbs" );
		}
		return 0;
	} // BENCHMARKING


	// main objects
	pose::Pose pose;
	protocols::hotspot_hashing::HotspotStubSet stubset;
	stubset.sc_only( option[hotspot::sc_only]() );

	// SET VARIABLES BASED ON THE COMMAND LINE
	// Residues to use for hashing (defaults to all, sans Gly and Pro)
	utility::vector1< std::string > resnames;
	if (option[ hotspot::residue ].user() ) {
		resnames = option[ hotspot::residue ]();
	}
	else resnames.push_back( "ALL" );

	// set scorefunction.
	core::scoring::ScoreFunctionOP scorefxn;
	if( option[ score::weights ].user() ) {
		scorefxn = core::scoring::getScoreFunction();
	} else {
		scorefxn = core::scoring::getScoreFunctionLegacy( "score13" );
		scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
		scorefxn->set_weight( core::scoring::envsmooth, 0 );
	}
	core::scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	TR << "Setting envinromental hbonds to: " << option[ hotspot::envhb]() << std::endl;
	options.hbond_options().use_hb_env_dep( option[ hotspot::envhb]() );
	scorefxn->set_energy_method_options( options );


	// number of stubs to find (defalts to 1000)
	int n_stubs( 1000 ); // must be int, otherwise subtracting causes baaadness
	if (option[ run::n_cycles ].user() ) {
		n_stubs = option[ run::n_cycles ]();
	}

	// set threshold for finding hotspots (default -1.0)
	core::Real threshold = option[ hotspot::threshold ]();
	stubset.score_threshold( threshold );

	// options for hashing against specific areas on the target
	// defaulting target res to 0 lets us easily test for option setting
	core::Size target_resnum(0);
	core::Real target_distance(20.0);
	if( option[ hotspot::target_res ].user() ) target_resnum = option[ hotspot::target_res ]();
	if( option[ hotspot::target_dist ].user() ) target_distance = option[ hotspot::target_dist ]();

	// read in the pose of our target
	std::string target_fname;
	if ( option[hotspot::target].user() )
	{
		target_fname = option[ hotspot::target ]();
		core::import_pose::pose_from_pdb( pose, target_fname );
	}
	else
	{
		utility_exit_with_message("You must specify a target to hash using -target <filename>");
	}

	// option to read in an existing set
	std::string hashin_fname = "";
	if (option[hotspot::hashfile].user() )
	{
		hashin_fname = option[hotspot::hashfile]();
	}
	// option to write out new set
	std::string hashout_fname;
	if ( option[ out::file::o ].user() )
	{
		hashout_fname = option[ out::file::o ]();
	}
	else
	{
		hashout_fname = "hash.hsh";
	}
	TR << "Writing stubs to " << hashout_fname << std::endl;

	// read existing hashes
	if ( utility::file::file_exists( hashin_fname ) ) {
		if( (option[ out::overwrite ].user()) ) utility::file::file_delete( hashin_fname );
		else {
			stubset.read_data( hashin_fname );
			TR << "Found hash file " << hashin_fname << std::endl;
		}
	}
	if ( utility::file::file_exists( hashout_fname ) ) {
		if( (option[ out::overwrite ].user()) ) utility::file::file_delete( hashout_fname );
		else {
			stubset.read_data( hashout_fname );
			TR << "Found hash file " << hashout_fname << std::endl;
		}
	}

	// fill residues to hash based on convenience name.
	// ALL = everything but G & C
	if( std::find( resnames.begin(), resnames.end(), "ALL" ) != resnames.end() ) {
		resnames.erase( std::find( resnames.begin(), resnames.end(), "ALL" ) );
		resnames.push_back( "ALA" );
		resnames.push_back( "ARG" );
		resnames.push_back( "ASN" );
		resnames.push_back( "ASP" );
		resnames.push_back( "GLU" );
		resnames.push_back( "GLN" );
		resnames.push_back( "HIS" );
		resnames.push_back( "ILE" );
		resnames.push_back( "LEU" );
		resnames.push_back( "LYS" );
		resnames.push_back( "MET" );
		resnames.push_back( "PHE" );
		resnames.push_back( "SER" );
		resnames.push_back( "THR" );
		resnames.push_back( "TRP" );
		resnames.push_back( "TYR" );
		resnames.push_back( "VAL" );
	}

	TR << "Finding hotspots using residues ";
 	for( utility::vector1< std::string >::const_iterator it=resnames.begin() ; it!=resnames.end(); ++it ) {
		TR << *it << " ";
	}
	TR << "." << std::endl;

	// for each residue requested
	for( utility::vector1< std::string >::const_iterator it=resnames.begin() ; it!=resnames.end(); ++it ) {
		std::string resname = *it;

		TR << "Hash contains " << stubset.size(resname) << " " << resname << " stubs." << std::endl;

		// check to see if we've already finished our hash
		int stubs_left = n_stubs;
		stubs_left -= stubset.size( resname );

		// if we're done
		if( stubs_left <= 0 )
		{
			TR << "No more " << resname << " hotspots to be found." << std::endl;
			// perform a scorecut
			if ( basic::options::option[ basic::options::OptionKeys::out::scorecut ].user() )
			{
				Real score_cutoff = option[ out::scorecut ]();
				std::stringstream i;
				i.str("");
				i << score_cutoff;
				stubset.clear();
				stubset.read_data( hashout_fname );
				protocols::hotspot_hashing::HotspotStubSetOP cut_stubs = stubset.subset( score_cutoff );
				std::string newfname = i.str() + "cut_" + hashout_fname;
				cut_stubs->write_all( newfname );
			}
			continue;
		}

    TR << "Finding " << stubs_left << " " << resname << " hotspots." << std::endl;

		// do hashing in 10-stub cycles to minimize file i/o. If we named a target residue, do 1-stub cycles (targets are a rare find)
		Size n_per;
		if( target_resnum ) n_per = 1;
		else n_per = 10;
		Size n_cycles = n_stubs / n_per;
		// make sure we do at least one cycle
		if( n_cycles <= 0 ) n_cycles = 1;
		// PERFORM HASHING
		for( Size i = 1; i <= n_cycles; ++i )
		{
			Size const length( option[ hotspot::length ]() );
			stubset.clear();
			stubset.hotspot_length( length );
			if( target_resnum ) stubset.fill( pose, scorefxn, target_resnum, target_distance, resname, n_per );
			else stubset.fill( pose, scorefxn, resname, n_per );
			stubset.write_all( hashout_fname );
		}
	} // for each residue

	// perform a scorecut
	if ( option[ out::scorecut ].user() )
	{
		Real score_cutoff = option[ out::scorecut ]();
		std::stringstream i;
		i.str("");
		i << score_cutoff;
		stubset.clear();
		stubset.read_data( hashout_fname );
		protocols::hotspot_hashing::HotspotStubSetOP cut_stubs = stubset.subset( score_cutoff );
		std::string newfname = i.str() + "cut_" + hashout_fname;
		cut_stubs->write_all( newfname);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // end main
