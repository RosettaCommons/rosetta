// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/*
 * iterative_design.cc
 *
 *  Created on: Nov 18, 2008
 *      Author: bcorreia
 */

/* parse loops file Loops
	loops::Loops loops;
	std::string filename( option[ OptionKeys::loops::loop_file ]().name() );
	loops.read_loop_file( filename );   // <== TODO: select these using density score

	check stuff about the job distributor*/

#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
// Auto-header: duplicate removed #include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopMover_CCD.hh>
// Auto-header: duplicate removed #include <protocols/loops/looprelax_protocols.hh>
#include <protocols/viewer/viewers.hh>
//#include <protocols/frags/TorsionFragment.hh>

#include <protocols/jobdist/Jobs.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include "core/scoring/packstat/compute_sasa.hh"


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>

#include <basic/Tracer.hh>
static THREAD_LOCAL basic::Tracer TR( "apps.iterative_design" );

// Utility Headers
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
using utility::string_split;
#include <numeric/xyzVector.hh>

# include <utility/file/gzip_util.hh>


// c++ headers
#include <fstream>
#include <iostream>
#include <utility>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <ObjexxFCL/string.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace chemical;
using namespace pack;
using namespace scoring;
using namespace protocols;


// Job Distributor -> there is not many advantages at this point
/*utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();


		for (core::Size jobnum = 1; jobnum <= input_jobs.size(); ++jobnum) {

				TR << "Processing " << input_jobs[jobnum]->input_tag() << "..." << std::endl;

		}
*/
// Checkpoint Stuff
// Iteration of perturbation-Design


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

		devel::init( argc, argv ); // reading options--name should be more descriptive


		TR << "Getting input filename(s)" << '\n';

		typedef vector1< utility::file::FileName > Filenames;
		Filenames pdbnames;

		typedef vector1< pose::Pose > Poses;

		Poses input_poses;


		if ( option[ in::file::l ].user() ) {
				Filenames listnames( option[ in::file::l ]().vector() );
				for ( Filenames::const_iterator filename( listnames.begin() );
				      filename != listnames.end(); ++filename ) {
					std::ifstream list( (*filename).name().c_str() );
					while ( list ) {
						std::string pdbname;
						list >> pdbname;
						pdbnames.push_back( pdbname );
					}
				}

			} else if ( option[ in::file::s ].user() ) {
				pdbnames = option[ in::file::s ]().vector();

			} else {
				std::cerr << "No files given: Use either -file:s or -file:l "
				          << "to design a single pdb or a list of pdbs"
				          << std::endl;
			}


		TR << "Reading Loop File" << '\n';

		protocols::loops::Loops loops;
		std::string filename(  protocols::loops::get_loop_file_name()  );
		loops.read_loop_file( filename );


		//Reading Natro file

		std::vector<Size> natro_res;

		if ( option[ OptionKeys::loops::keep_natro ].user() ){
			std::string filename( option[ OptionKeys::loops::keep_natro ]().name() );
			std::ifstream infile( filename.c_str() );
			TR << "Reading keep_natro file " << '\n';
			std::string res;

			while (getline(infile,res)){
				Size conv_res = (Size) atoi(res.c_str());
				TR << "Residues to keep "<< conv_res <<"\n";
				natro_res.push_back( conv_res );
				}

			}


		for ( Filenames::const_iterator filename( pdbnames.begin() );
					filename != pdbnames.end(); ++filename ) {
				std::cout << "Input pdb file " << *filename;
				if ( !utility::file::file_exists( *filename ) ) {
					std::cout << " not found, skipping" << std::endl;
					continue;
				}
				std::cout << std::endl;

				std::string pdbprefix( string_split( string_split( *filename, '/' ).back(), '.' ).front() );


				pose::Pose pose;
				core::import_pose::pose_from_file( pose, *filename , core::import_pose::PDB_file);
				input_poses.push_back( pose );


				std::cout << pdbprefix  <<std::endl;
			}


		Size refine_design_iterations;

		if ( option[ OptionKeys::loops::refine_design_iterations ].user() ){

			refine_design_iterations = option[ OptionKeys::loops::refine_design_iterations ]();
		}


		else{

			refine_design_iterations = 1;
		}


		std::string pdb_prefix;

		if ( option[ out::prefix ].user() ){
			pdb_prefix= option[ out::prefix ]()+"_"+"S";
			}

		else {

			pdb_prefix="S";

		}

		core::scoring::ScoreFunctionCOP scorefxn( get_score_function() );


		//setting a minimizer objet that will be ran after each design step

		core::optimization::AtomTreeMinimizer mzr;
		core::optimization::MinimizerOptions options("lbfgs_armijo_nonmonotone", 1e-5, true, false);
		core::kinematics::MoveMap mm;
		mm.set_bb(false);


		int pose_number = 0;

		pose::Pose init_pose;
		pose::Pose nat_pose;

		core::io::silent::SilentFileData sfd;


		for ( Poses::const_iterator iter = input_poses.begin(); iter != input_poses.end(); ++iter)
		{


			std::string filename_init (pdb_prefix+ "_" + ObjexxFCL::right_string_of(pose_number,3,'0')+"_init");


			TR << "***** Starting full-atom loop refinement protocol  ****" << std::endl;


			nat_pose = *iter;
			init_pose = *iter;

			(*scorefxn)(init_pose);


			setPoseExtraScore( init_pose, "packing", core::scoring::packstat::compute_residue_packing_score( init_pose, 2 ) );

			setPoseExtraScore( init_pose, "loop_rms", 0 );

			core::io::silent::BinarySilentStructOP pss_init (
											new core::io::silent::BinarySilentStruct ( init_pose , filename_init ));

			sfd.add_structure( pss_init );

			//(*scorefxn)(init_pose);


			//list of neighbors of the loops for a giving pose

			std::vector<Size> loop_residues;
			std::vector<Size> loops_nbr;


			for (core::Size i=1; i<=init_pose.size(); ++i ) {

				for ( int j = 1; j <= (int)loops.size(); ++j ) {
					if ( i >= core::Size( loops[j].start() ) && i <= core::Size( loops[j].stop() ) ) {
						TR << "Repacking because in loop: " << i << std::endl;
						loop_residues.push_back(i) ;
						}
					}
				}//vector of loop residues


			for ( int iter_res= 1; iter_res < (int)loop_residues.size(); ++iter_res  ){


				Vector const  & nbr_atom( init_pose.residue(loop_residues[iter_res]).nbr_atom_xyz() );

				for ( Size i=1 ; i <= init_pose.size() ; ++i){

					//bool isloop = false;

					if ( nbr_atom.distance( init_pose.residue(i).nbr_atom_xyz() ) < 5.0 ){
						loops_nbr.push_back(i);
						}

					}


				}// vector of the neighbors


			//clean redundancy on the neighbors vector

			for (std::vector<Size>::const_iterator pos=loops_nbr.begin(); pos !=loops_nbr.end(); ++pos ){

				//TR <<"Loops before filtering " << *pos << std::endl;

			}
				TR << "Residues Before Filtering  "<< loops_nbr.size()<< std::endl;

			std::sort( loops_nbr.begin(), loops_nbr.end() );

			std::vector<Size>::iterator new_end_pos;
			new_end_pos = std::unique( loops_nbr.begin(), loops_nbr.end() );

			loops_nbr.erase( new_end_pos, loops_nbr.end() );

			for (std::vector<Size>::iterator pos=loops_nbr.begin(); pos !=loops_nbr.end(); ++pos ){
				//bool redundant_res( false );

				for (std::vector<Size>::const_iterator pos_residues=loop_residues.begin(); pos_residues !=loop_residues.end(); ++pos_residues){


					if (*pos == *pos_residues && pos < loops_nbr.end()-1 ){

						loops_nbr.erase( pos );
					}


				}

				if  ( option[ OptionKeys::loops::keep_natro ].user() ){

					for (std::vector<Size>::const_iterator pos_natro= natro_res.begin(); pos_natro != natro_res.end(); ++pos_natro ){

						if (*pos == *pos_natro){
							loops_nbr.erase( pos );
						}
					}
				}//cleaning the natro res


				//TR <<"Loops after filtering " << *pos << std::endl;

			}

				TR << "Residues after filtering "<< loops_nbr.size()<< std::endl;


			if ( option [OptionKeys::loops::keep_natro ].user() ){
				for (std::vector<Size>::iterator pos_loop= loop_residues.begin(); pos_loop != loop_residues.end(); ++pos_loop ){

					for (std::vector<Size>::const_iterator pos_natro = natro_res.begin(); pos_natro != natro_res.end(); ++pos_natro){

						if (*pos_loop == *pos_natro){loop_residues.erase( pos_loop ); }

					}


				}


			}


			for (Size iteration = 1 ; iteration <= refine_design_iterations ; ++iteration ) {


				loops.auto_choose_cutpoints( init_pose );

				kinematics::FoldTree f_new, f_orig=init_pose.fold_tree();

				protocols::loops::fold_tree_from_loops( init_pose, loops, f_new );
				init_pose.fold_tree( f_new );


				protocols::loops::LoopMover_Refine_CCD refine_ccd( loops);
				refine_ccd.set_native_pose( new core::pose::Pose (init_pose) );


				refine_ccd.apply( init_pose );

				std::string filename_ref (pdb_prefix + "_" + ObjexxFCL::right_string_of(pose_number,3,'0')+"_"+ObjexxFCL::right_string_of(iteration ,3,'0') + "_ref");


				setPoseExtraScore( init_pose, "packing", core::scoring::packstat::compute_residue_packing_score( init_pose, 2 ) );

				setPoseExtraScore( init_pose, "loop_rms", protocols::loops::loop_rmsd(nat_pose, init_pose, loops ));


				core::io::silent::BinarySilentStructOP pss_ref (
									new core::io::silent::BinarySilentStruct ( init_pose , filename_ref ));

				sfd.add_structure( pss_ref );

				//design neighbors using soft_rep weights

				utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, true );

				utility::vector1< bool > residues_to_mutate( init_pose.size(), false );


				core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( init_pose ));


				for (core::Size j=0; j < loop_residues.size(); ++j ) {

						TR <<"Loop to design " << loop_residues[j] << std::endl;

						task->nonconst_residue_task(loop_residues[j]).restrict_absent_canonical_aas( allowed_aas );

						residues_to_mutate[loop_residues[j]] = true;

						mm.set_chi(loop_residues[j],true);


					}//Design loops


				for (core::Size j=0; j < loops_nbr.size(); ++j ) {

						TR <<"neighbors to design " << loops_nbr[j] << std::endl;

						task->nonconst_residue_task(loops_nbr[j]).restrict_absent_canonical_aas( allowed_aas );

						residues_to_mutate[loops_nbr[j]] = true;

						mm.set_chi(loops_nbr[j],true);

				}//Designing Neighbors


			task->restrict_to_residues( residues_to_mutate );

			task->initialize_extra_rotamer_flags_from_command_line();

			pack::pack_rotamers( init_pose, *scorefxn , task );


			//Place To minimize designed side chains

			mzr.run( init_pose, mm, *scorefxn, options );

			mm.set_chi(false);


			protocols::loops::remove_cutpoint_variants( init_pose );

			std::string filename_ref_des (pdb_prefix +"_"+ ObjexxFCL::right_string_of(pose_number,3,'0')+"_"+ObjexxFCL::right_string_of( iteration,3,'0') + "_ref_des");


			(*scorefxn)(init_pose);

			setPoseExtraScore( init_pose, "packing", core::scoring::packstat::compute_residue_packing_score( init_pose, 2 ) );

			setPoseExtraScore( init_pose, "loop_rms", 0 );


			core::io::silent::BinarySilentStructOP pss (
					new core::io::silent::BinarySilentStruct ( init_pose , filename_ref_des ));


			sfd.add_structure( pss );


			}// iteration on the refine_design


			std::string pdb_silent_file;

			if ( option[ out::file::silent ].user() ){
				pdb_silent_file= option[ out::file::silent ]();
				}

			else {
				pdb_silent_file = "test";
			}


			sfd.write_all(pdb_silent_file);


			sfd.clear();

			++pose_number;
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
