// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers


#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/PDBSilentStruct.hh>


#include <core/io/silent/silent.fwd.hh>
// Auto-header: duplicate removed #include <core/io/silent/ProteinSilentStruct.hh>
// Auto-header: duplicate removed #include <core/io/silent/SilentFileData.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/comparative_modeling/ConstraintRemodelMover.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Auto-header: duplicate removed #include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/io/mpistream.hh>


int
main( int argc, char* argv [] )
{
    try {
	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	// utility::vector1< protocols::jobdist::BasicJobOP > input_jobs;
	// int const nstruct_flag = option[ out::nstruct ];
	// int const nstruct = std::max( 1, nstruct_flag );
	//
	// protocols::jobdist::BasicJobOP job = new protocols::jobdist::BasicJob("" /*no input tag*/, "S", nstruct);
	// input_jobs.push_back( job );
	// //protocols::jobdist::PlainPdbJobDistributor< protocols::jobdist::BasicJobOP > jobdist( input_jobs );
	// protocols::jobdist::BaseJobDistributorOP jobdist;

	// jobdist->startup();
	// while ( jobdist->next_job(curr_job, curr_nstruct) ) { // loop over jobs

	// setup residue types

	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
	  rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
	  rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	//read in pose
	core::pose::Pose pose, native_pose;
	std::string input_filename = basic::options::start_file();
	core::import_pose::pose_from_file( pose, input_filename , core::import_pose::PDB_file);
	if (option[ in::file::native ].user() )
		core::import_pose::pose_from_file( native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);

	// setup ScoreFunction and constraints for minimization.
	core::scoring::constraints::ConstraintSetOP cstset
		= core::scoring::constraints::ConstraintIO::read_constraints(
		  core::scoring::constraints::get_cst_file_option(),
		  new core::scoring::constraints::ConstraintSet,pose
	);
	pose.constraint_set( cstset );


	//initialize score function
	using core::Real;
	Real init_rmsd    = core::scoring::all_atom_rmsd( native_pose, pose );
	//Real init_ca_rmsd = core::scoring::CA_rmsd( native_pose, pose );  //unused
	Real init_gdtmm   = core::scoring::CA_gdtmm( native_pose, pose );

	// protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" );

	//	core::Real initial_rms = core::scoring::all_atom_rmsd( native_pose, pose );
	protocols::comparative_modeling::ConstraintRemodelMover my_mover;
	my_mover.apply( pose );

	if ( option[ out::file::o ].user() ) {
		std::string output_name = option[ out::file::o ]();
		pose.dump_pdb(output_name);
	}

	Real final_rmsd = 0.0, final_gdtmm = 0.0;
	if ( option[ in::file::native ].user() ) {
		final_rmsd    = core::scoring::all_atom_rmsd( native_pose, pose );
		final_gdtmm   = core::scoring::CA_gdtmm( native_pose, pose );
	}

	// write score-file if specified by the user
	if ( option[ out::file::silent ].user() ) {
		core::io::silent::SilentFileData sfd;
		std::string scorefile( option[ out::file::silent ]() );
		std::string output_name = input_filename + ".min";
		core::io::silent::PDBSilentStruct pss( pose, output_name );

		pss.add_energy( "irms", init_rmsd );
		pss.add_energy( "rms" , final_rmsd );

		pss.add_energy( "iCA_rms", init_rmsd );
		pss.add_energy( "CA_rms" , final_rmsd );

		pss.add_energy( "igdtmm", init_gdtmm );
		pss.add_energy( "gdtmm" , final_gdtmm );
		sfd.write_silent_struct( pss, scorefile );
	} // if ( option[ out::file::silent ].user() )

	// jobdist->dump_pose_and_map( curr_job->output_tag(curr_nstruct), pose );    // output PDB


    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
} // int main
