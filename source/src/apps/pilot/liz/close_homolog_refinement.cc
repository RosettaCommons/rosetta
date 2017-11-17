// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
// #include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
// #include <core/chemical/VariantType.hh>
//
// #include <core/chemical/ChemicalManager.hh>

// #include <core/scoring/etable/Etable.hh>
// #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <core/pack/task/ResfileReader.hh>
#include <ObjexxFCL/format.hh>

#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>


// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/casp.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using basic::Error;
using basic::Warning;


using namespace core;


///////////////////////////////////////////////////////////////////////////////
bool
permutate( utility::vector1<int> &residues_to_randomize ) {

	int max_permutations = residues_to_randomize.size()*2;
	for ( int i=0; i < max_permutations; i++ ) {
		int index1 = std::rand() % ( residues_to_randomize.size() ) +1;
		int index2 = std::rand() % ( residues_to_randomize.size() ) +1;
		//choose two random numbers from 1-(vector_size)

		int temp_storage = residues_to_randomize[index1]; //switch the two numbers
		residues_to_randomize[index1]=residues_to_randomize[index2];
		residues_to_randomize[index2]=temp_storage;
	}
	return true;
}

int
main( int argc, char * argv [] )
{
	try {

		using namespace pose;
		using namespace core;
		using namespace scoring;
		using namespace conformation;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pack::task;
		using namespace optimization;

		// setup random numbers and options
		devel::init(argc, argv);

		// read the pose
		pose::Pose pose;
		// gets filename from -s option
		core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);

		//if optimization_radius > 100, then repack/sc_min entire structure
		core::Real optimization_radius = option[ casp::opt_radius ]();
		bool repack = option[ casp::repack ]();
		bool sidechain_min = option[ casp::sc_min ]();
		bool sequential = option[ casp::sequential ](); //whether or not mutations should be treated sequentially or all at once
		//mutations will be randomized
		core::Size iterations = option[ casp::num_iterations ]();
		//std::string weight_file = option[ casp::weight_file ]();
		//ScoreFunctionOP scorefxn(ScoreFunctionFactory::create_score_function(weight_file));
		ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		//read in file specifying which residues to consider
		std::string residue_spec_file = option[ casp::refine_res ]();
		utility::io::izstream data(residue_spec_file);

		int resnum;
		std::string line;
		utility::vector1<int> resnums_to_refine;
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open residue specification file!"+residue_spec_file );
		}
		while ( getline(data,line) ) {
			std::istringstream l(line);
			l >> resnum;
			resnums_to_refine.push_back(resnum);
		}
		//end read in file

		while ( iterations > 0 ) {
			std::cout << "current iteration: " << iterations << std::endl;
			iterations--;
			permutate(resnums_to_refine); //randomize the numbers within this vector

			for ( core::Size i = 1; i <= resnums_to_refine.size(); i++ ) {
				//create and set up packertask
				pack::task::PackerTaskOP task(pack::task::TaskFactory::create_packer_task(pose));
				task->restrict_to_repacking();
				utility::vector1<bool> residues_to_refine(pose.size(),false);

				//set up movemap
				kinematics::MoveMap movemap;
				int current_residue = resnums_to_refine[i];

				//if sequential, then go through mutations one at a time
				if ( sequential && (optimization_radius < 100) ) { //sequentially consider each mutation individually
					residues_to_refine[current_residue]=true;
					Vector const & nbr_atom( pose.residue(current_residue).nbr_atom_xyz() ); //neighbors of atom i

					if ( sidechain_min ) {
						movemap.set_chi(current_residue,true); // if sidechain minimization is wanted, set residue i to sidechain minimize
					}

					for ( Size j=1; j<= pose.size(); ++j ) {
						if ( nbr_atom.distance( pose.residue(j).nbr_atom_xyz() ) < optimization_radius ) {

							residues_to_refine[j] = true;
							//set to minimize the sidechains of this residue
							if ( sidechain_min ) {
								movemap.set_chi(j,true);
							}
						}
					}
					task->restrict_to_residues(residues_to_refine);
				} else if ( optimization_radius < 100 ) { // if ( sequential && (optimization_radius < 100) )
					//consider all residues in array resnums to refine simultaneously
					//then break loop
					for ( core::Size j=1; j <= resnums_to_refine.size(); j++ ) {
						int current_residue = resnums_to_refine[j];
						residues_to_refine[current_residue]=true;
						Vector const & nbr_atom( pose.residue(current_residue).nbr_atom_xyz() ); //neighbors of current residue
						for ( Size k=1; k<= pose.size(); ++k ) {
							if ( nbr_atom.distance( pose.residue(k).nbr_atom_xyz() ) < optimization_radius ) {
								residues_to_refine[k]=true;
							}
							if ( sidechain_min ) {
								movemap.set_chi(k,true);
							}
						}

						if ( sidechain_min ) {
							movemap.set_chi(current_residue,true);
						}
					}
					task->restrict_to_residues(residues_to_refine);

				} else {
					//all residues to be refined
					task->restrict_to_repacking();
				}

				//repack
				//this is a little silly.. why did we go through all that trouble
				//if we weren't going to repack? But this was a silly little option to easily control program flow
				if ( repack ) {
					pack::pack_rotamers(pose, (*scorefxn), task);
				}

				//sidechain minimize
				if ( sidechain_min ) {
					std::string const min_type( option[ run::min_type ]() );
					AtomTreeMinimizer().run( pose, movemap, (*scorefxn), MinimizerOptions( min_type, 0.001, true ) );
				}

				if ( !sequential ) { //at the end of refinement, break if not sequential
					break; //we've already done all our refinement
				}
			} // for ( int i = 1; i <= resnums_to_refine.size(); i++ )
		} // while ( iterations > 0 )

		std::string output_name = option[ out::file::o ]();
		pose.dump_pdb(output_name);

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
} // int main
