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
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <protocols/moves/Mover.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.yiliu.silent" );
// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>


class SilentProlineFixMover : public protocols::moves::Mover {

public:
		SilentProlineFixMover()  {}
		~SilentProlineFixMover() {}

virtual void apply( core::pose::Pose & pose ) {

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );

	// try a repack
	core::pack::task::PackerTaskOP task(
		core::pack::task::TaskFactory::create_packer_task( pose )
	);
	task->initialize_from_command_line();
	task->restrict_to_repacking();
	core::pack::pack_rotamers(pose, (*scorefxn), task);

	// try a minimize
	core::kinematics::MoveMap mm;
	mm.set_bb(false);
	mm.set_chi(true);

	core::optimization::AtomTreeMinimizer().run(
		pose, mm, (*scorefxn), core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true)
	);

} // apply

}; // class SilentProlineFixMover


int
main( int argc, char* argv [] )
{
	try {
	// options, random initialization
	devel::init( argc, argv );

	using std::string;
	using utility::vector1;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	// configure silent-file data object
	core::io::silent::SilentFileData sfd;
	std::string infile  = *( option[ in::file::silent ]().begin() );
	utility::vector1< string > tags;
	//TR << "option[ in::file::tags ].()" <<option[ in::file::tags ]()<<std::endl;
	if ( option[ in::file::tags ].user() ) {
		tags = option[ in::file::tags ]();
	}
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			TR << "option[ in::file::tags ].(), if condition" <<option[ in::file::tags ]()<<std::endl;
			sfd.read_file( infile, option[ in::file::tags ]() );
		} else {
			TR << "option[ in::file::tags ].(), else condition" <<option[ in::file::tags ]()<<std::endl;
			sfd.read_file( infile );
		}
	} else {
		utility_exit_with_message(
			"Error: can't get any structures! Use -in::file::silent <silentfile>"
		);
	}

	core::pose::Pose pose;
	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		iter->fill_pose( pose, *rsd_set );
		pose.dump_pdb( iter->decoy_tag() + ".pdb" );
		// SilentProlineFixMover fixer_upper;
		// fixer_upper.apply( pose );
		// pose.dump_pdb( iter->decoy_tag() + ".fix_repack_min.pdb" );
	}

	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
} // main
