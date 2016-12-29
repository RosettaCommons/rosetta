// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/** @page MembraneJumpExample
Read a membrane protein PDB with at least one jump, move a chain (this is a rigid body move; no packing),
print the PDB
Try:
"TMHTopologySampler.cc -in::file::s <pdb file> -in::file::spanfile <spanfile>
-database <database> -rigid::fragment_cycles <ncycles>  -rigid::rotation <angle> -rigid::translation>
-run::show_simulation_in_pymol <#seconds> -out::prefix <prefix>"
*/
/// @file   apps/pilot/shirst/MembraneJumpExample.cc
///
/// @brief This is to illustrate reading in of PDB file, invoking MembraneTopology, and making a Jump
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB with at least 2 residues to pack>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::spanfile <spanfile>
/// @author Stephanie DeLuca (stephanie.h.deluca@vanderbilt.edu)

//Scoring
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoreFunction.hh>

//Pose
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <utility/excn/EXCN_Base.hh>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/PoseMembraneRigidBodyMover.hh>
#include <protocols/rigid/PoseMembraneRigidBodyMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/scoring/MembranePotential.hh>

// Project Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedStubID.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <sstream>

#include <devel/init.hh>

static THREAD_LOCAL basic::Tracer tr( "apps.pilot.TMHTopologySampler" );

//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	try {
		devel::init(argc, argv);

		using namespace basic::options::OptionKeys;
		using basic::options::option;


		//pose is read in from PDB for now
		utility::vector0<std::string> pdbs;
		pdbs = option[in::file::s]();
		std::string pdb=pdbs[0];

		core::pose::Pose pose; // starts NULL, coords *never* modified!
		core::import_pose::pose_from_file(pose, pdb, core::import_pose::PDB_file);

		// Can we add the PyMOL mover here?
		{
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[OptionKeys::run::show_simulation_in_pymol].user()
					&& option[OptionKeys::run::show_simulation_in_pymol].value() > 0.0 ) {
				protocols::moves::AddPyMOLObserver(pose,
					option[OptionKeys::run::keep_pymol_simulation_history](),
					option[OptionKeys::run::show_simulation_in_pymol].value());
			}
		}
		//At this point, assert if don't have a spanfile
		std::string spanfile = option[in::file::spanfile];
		runtime_assert (option[in::file::spanfile].user());

		//set up membrane topology
		core::scoring::MembraneEmbedOP membrane_embed( new core::scoring::MembraneEmbed );
		core::scoring::MembraneTopologyOP topology( new core::scoring::MembraneTopology );
		pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topology );
		pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED, membrane_embed);
		topology->initialize(spanfile);

		//put the protein in the membrane center
		protocols::rigid::MovePoseToMembraneCenterMover center_pose;
		center_pose.apply(pose);

		core::Size const nspan = topology->tmhelix();
		core::Size const njumps = nspan-1;
		core::Size const num_cut_loops = nspan-1;
		core::Size const nres = pose.size();

		// this function will setup a fold tree to be used given a spanfile.  It basically splits up the pose into TMHs as defined
		{
			core::kinematics::FoldTree fold_tree(nres);

			tr << "setting up fold_tree with " << (njumps) << " jump(s) in a total of " << num_cut_loops << " loops to cut\n";

			// exit if spanfile not formatted correctly
			if ( nspan==0 ) {
				utility_exit_with_message("bad format for spanfile total_tmhelix==0");
			}

			//exit if nspan is not more than njumps and if # of cuttable loops != # jumps
			assert(nspan>=njumps);
			assert(njumps==num_cut_loops);

			//put all residues in spans in a one-dimensional array that is nres long, and 0=not in span, 1-3=is in span
			ObjexxFCL::FArray1D_int span(nres,0);
			for ( core::Size i=1; i<=nres; ++i ) {
				for ( core::Size j=1; j<=nspan; ++j ) {
					if ( i >= topology->span_begin(j) && i <= topology->span_end(j) ) {
						span(i)=j;
					}
				}
			}

			//Get info to make sure membrane span information is being read in correctly
			tr.Debug << "size of span array:  " << span.size() << "\tnumber of residues:  " << nres << std::endl;
			for ( core::Size i=1; i<=span.size(); ++i ) {
				tr.Debug << "residue:  " << i << " span:  " << span(i) << std::endl;
			}

			//Make an array of start and end points of loops where we can cut. Needed by FoldTree::random_tree_from_jump_points();
			// initialize some variables for setting up loops array
			core::Size span_index = 1;
			ObjexxFCL::FArray1D_int previous_span_begin((nspan-1),0);
			ObjexxFCL::FArray1D_int previous_span_end((nspan-1),0);
			ObjexxFCL::FArray1D_int span_begin((nspan-1),0);
			ObjexxFCL::FArray1D_int span_end((nspan-1),0);
			ObjexxFCL::FArray2D_int cut_loops(2,num_cut_loops);
			ObjexxFCL::FArray1D_int loop_begin(num_cut_loops,0);
			ObjexxFCL::FArray1D_int loop_end(num_cut_loops,0);

			// set up cut loops array
			for ( span_index = 1; span_index <= nspan-1; ++span_index ) {
				//need to know the beginning and end of the TMs to determine where the loops are
				previous_span_begin(span_index) = topology->span_begin(span_index);
				previous_span_end(span_index) = topology->span_end(span_index);
				span_begin(span_index) = topology->span_begin(span_index+1);
				span_end(span_index) = topology->span_end(span_index+1);

				loop_begin(span_index) = previous_span_end(span_index) + 1;
				loop_end(span_index) = span_begin(span_index) - 1;

				//if the predicted loop (that is, not span) is shorter than 3 res
				if ( loop_end(span_index) - loop_begin(span_index) < 2 ) {
					loop_begin(span_index) -= 1;
					loop_end(span_index) += 1;
				}

				assert(loop_begin(span_index) != 0);
				assert(loop_end(span_index) != 0);

				// remember that loops is a 2D array that is 2 x num_cut_loops.
				for ( core::Size cut_loop_index = 1; cut_loop_index <= num_cut_loops; ++cut_loop_index ) {
					cut_loops(1,cut_loop_index) = loop_begin(cut_loop_index);
					cut_loops(2,cut_loop_index) = loop_end(cut_loop_index);

					tr.Info << "cut_loops_array:  " << cut_loops(1,cut_loop_index) << " "
						<< cut_loops(2,cut_loop_index) << std::endl;
				}

				tr.Info << "span_index:  " << span_index << " previous_span_begin:  " << previous_span_begin(span_index) << " previous_span_end: "
					<< previous_span_end(span_index) << std::endl;
				tr.Info << "span_index:  " << span_index << " span_begin:  " << span_begin(span_index) << " span_end: "
					<< span_end(span_index) << std::endl;
			}//for span_index

			// set cut bias by checking to see if residue i is in a TMH, if it is, cut_bias=0, if so, cut_bias=1
			ObjexxFCL::FArray1D_float cut_bias(nres,0.0);
			for ( core::Size i = 1; i <= nres; ++i ) {
				if ( span(i)==0 ) {
					cut_bias(i)=1;
				}
			}

			//print out cut_bias array
			for ( core::Size i=1; i<=cut_bias.size(); ++i ) {
				tr.Debug << "cut_bias_array:  " << i << " " << cut_bias(i) << std::endl;
			}

			// Generate the fold tree
			fold_tree.random_tree_from_jump_points(nres,njumps,cut_loops,cut_bias);
			fold_tree.put_jump_stubs_intra_residue();

			// output fold tree
			tr << fold_tree;
			pose.fold_tree(fold_tree);
		}

		///  Randomly Perturb with a default rotational magnitude of 30 and a translational
		///  magnitude of 50.0 angstroms; You can change this with constructor args.
		core::Real rotation_mag;
		core::Real translation_mag;
		core::Size ncycles;
		ncycles = option[rigid::fragment_cycles]();
		rotation_mag = option[rigid::rotation]();
		translation_mag = option[rigid::translation]();

		tr << "num_move_cycles:  " << ncycles << " rotation_mag:  " << rotation_mag << " translation_mag:  " << translation_mag << std::endl;

		//store all current jumps in the pose in a vector of pairs
		assert(njumps==pose.num_jump());
		utility::vector1<std::pair<core::Size,core::kinematics::Jump> > jumps;
		for ( core::Size jump_num = 1; jump_num <= pose.num_jump(); ++jump_num ) {
			jumps.push_back(std::make_pair(jump_num,pose.jump(jump_num)));
		}

		const std::string prefix = (option[out::prefix]());

		for ( core::Size i=1; i<=ncycles; ++i ) {
			core::Size random_jump_num = static_cast<core::Size>(numeric::random::rg().random_range(1,njumps));
			tr.Debug << "random_jump_num:  " << random_jump_num << std::endl;
			protocols::rigid::RigidBodyPerturbMover mover(random_jump_num,rotation_mag,translation_mag);
			mover.apply(pose);

			// output the states after each rigid body move
			//pose.dump_pdb("cycle_"+utility::to_string(i)+".pdb","pdb");
			pose.dump_pdb("cycle_"+utility::to_string(i)+".pdb");

			//a jump in the pose stores rotational and translational information, and therefore the conformation of the pose.
			//Therefore, reverting jumps (not the jump number) will reset the pose.
			//We're not changing the FoldTree!
			assert(jumps.size() == pose.num_jump());
			utility::vector1<std::pair<core::Size,core::kinematics::Jump> >::const_iterator jump_it;
			for ( jump_it = jumps.begin(); jump_it != jumps.end(); ++jump_it ) {
				pose.set_jump(jump_it->first,jump_it->second);
			}

			pose.update_actcoords();

			// output the states after resetting jumps
			//pose.dump_pdb("cycle_"+utility::to_string(i)+"_reset.pdb","pdb");
			pose.dump_pdb("cycle_"+utility::to_string(i)+"_reset.pdb");
		}

		// Output final PDB
		const std::string output_final = (prefix+"_final.pdb");

		//pose.dump_pdb(output_final,pdb);
		pose.dump_pdb(output_final);

		tr << "Finished TMH sampling.  Outputting " << output_final << std::endl;
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return(0);
}
