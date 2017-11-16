// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <protocols/rigid/RigidBodyMover.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "apps.pilot.whitney_repack_input.main" );


/// General testing code
int
main( int argc, char * argv [] )
{

	try {


		devel::init(argc, argv);

		TR << "Starting repacking input" << std::endl;

		// create pose for wild type NGF/p75 complex and unbound pose and wild type NGF/TrkA complex and unbound
		pose::Pose NGF_pose, NT3_pose;
		core::import_pose::pose_from_file( NGF_pose, "/Users/whitneyS/Desktop/threading/NGF_alone.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( NT3_pose, "/Users/whitneyS/Desktop/threading/NT3.pdb" , core::import_pose::PDB_file);

		// create and open file
		utility::io::ozstream ddg_outstream, ddg_outstream2;

		ddg_outstream.open( "compare_NGF.out", std::ios::out );
		ddg_outstream2.open( "compare_NT3_rev.out", std::ios::out );

		// scoring function
		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		scorefxn->set_weight( core::scoring::fa_dun, 0.1 );

		// Setup packer task for NT3
		pack::task::PackerTaskOP NT3_base_packer_task( pack::task::TaskFactory::create_packer_task( NT3_pose ));
		NT3_base_packer_task->set_bump_check( false );
		NT3_base_packer_task->initialize_from_command_line();
		NT3_base_packer_task->or_include_current( true ); // allow crystallized position as possible rotamer

		TR << "Residues in NGF: " << NGF_pose.size() << std::endl;
		TR << "Residues in NT3: " << NT3_pose.size() << std::endl;

		for ( int ii = 1, resnum = NT3_pose.size(); ii <= resnum; ++ii ) {

			// allow only the current position to be changed
			pack::task::PackerTaskOP NT3_position_packer_task( NT3_base_packer_task->clone() );
			utility::vector1 <bool> allow_pos( resnum, false );
			allow_pos.at(ii) = true;
			NT3_position_packer_task->restrict_to_residues( allow_pos );

			chemical::AA const NGF_aa ( NGF_pose.residue(ii).aa() );
			pack::task::PackerTaskOP NT3_repacked_packer_task( NT3_position_packer_task->clone() );
			utility::vector1 <bool> repack_aalist( core::chemical::num_canonical_aas, false );
			repack_aalist[ NGF_aa ] = true;
			NT3_repacked_packer_task->nonconst_residue_task(ii).restrict_absent_canonical_aas( repack_aalist );
			pack::pack_rotamers( NT3_pose, *scorefxn, NT3_repacked_packer_task );

		}

		//writes out new pdb file for unbound TrkA
		std::ostringstream outPDB_name;
		outPDB_name << "NGF_seq_NT3.pdb";
		NT3_pose.dump_scored_pdb( outPDB_name.str(), *scorefxn );

		for ( int j = 1, numres = NT3_pose.size(); j <= numres; ++j ) {
			chemical::AA const NGF_res( NGF_pose.residue(j).aa() );
			std::ostringstream data_string;
			data_string << oneletter_code_from_aa(NGF_res);
			ddg_outstream << data_string.str();
		}

		ddg_outstream.close();
		ddg_outstream.clear();

		for ( int i = 1, numres = NT3_pose.size(); i <= numres; ++i ) {
			chemical::AA const NT3_res( NT3_pose.residue(i).aa() );
			std::ostringstream data_string2;
			data_string2 << oneletter_code_from_aa(NT3_res);
			ddg_outstream2 << data_string2.str();
		}

		// close file
		ddg_outstream2.close();
		ddg_outstream2.clear();


		/*
		for ( int ii = 1, resnum = NGF_pose.size(); ii <= resnum; ++ii ) {


		chemical::AA const NT3_aa( NT3_pose.residue(ii).aa() );

		std::ostringstream data_string_stream;
		data_string_stream << std::setw(9) << ii;
		data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << NGF_aa;
		data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << NT3_aa;
		ddg_outstream << data_string_stream.str() << std::endl;
		}
		*/


		TR << "Successfully finished repacking input." << std::endl;

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


