// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file - Based on fixbb.cc
/// @brief - eventually should cycle between relaxation of structure while enabling keeping
/// @        some parts fixed to designing the sequence
/// @author Avital and Nir

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/util.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <protocols/relax_protocols.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/jobdist/standard_mains.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

static THREAD_LOCAL basic::Tracer TR( "avital.stabilize" );

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
	using namespace basic;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves;
	using namespace protocols;


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	devel::init(argc, argv);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// read the pose
	 core::pose::Pose pose;
	 core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file); // gets filename from -s option

	using namespace core::pack::task;
	TaskFactoryOP main_task_factory = new TaskFactory;
	main_task_factory->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new operation::ReadResfile );
	}

	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();


	MoverOP relax_mover = new ClassicRelax ( score_fxn );

	// protocol
	SequenceMoverOP full_seq = new SequenceMover;
	full_seq->add_mover( pack_mover );
	//full_seq->add_mover( relax_mover );

	std::map < std::string, core::Real > score_map;
	//jobdist->dump_pose_and_map( curr_job->output_tag(curr_nstruct), *the_pose )
	score_map["trying"]=(*score_fxn)(pose);
	//	pack_mover->apply(pose);
	//	pose.dump_pdb("design.pdb");
	TR<<(*score_fxn)(pose)<<std::endl;
	TR << "This structure rules!!! " << std::endl;

	protocols::jobdist::main_plain_pdb_mover( *full_seq, score_fxn );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}

