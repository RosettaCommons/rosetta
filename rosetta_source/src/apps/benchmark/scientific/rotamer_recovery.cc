// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/benchmark/scientific/rotamer_recovery.cc
/// @brief this determines the percent of rotamers correctly predictied during a repack
/// @author Matthew O'Meara (mattjomeara@gmail.com)



// Unit Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/RotamerRecoveryMover.hh>


// Project Headers
#include <devel/init.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using std::endl;
using std::string;

using basic::Tracer;
using core::pose::PoseOP;
using core::import_pose::pose_from_pdb;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::InitializeFromCommandline;
using core::pack::task::operation::ReadResfile;
using core::pack::task::operation::RestrictToRepacking;
using core::scoring::getScoreFunction;
using core::scoring::ScoreFunctionOP;
using devel::init;
using protocols::moves::RotamerRecoveryMover;
using protocols::moves::RotamerRecoveryMoverOP;
using protocols::jd2::JobDistributor;
using utility::excn::EXCN_Base;

static Tracer TR("apps.benchmark.scientific.rotamer_recovery");

OPT_1GRP_KEY( String, rotamer_recovery, reporter )
OPT_1GRP_KEY( String, rotamer_recovery, comparer )
//OPT_2GRP_KEY( File, rotamer_recovery, native_vs_decoy )
//OPT_2GRP_KEY( File, out, rotamer_recovery, method )
OPT_2GRP_KEY( File, out, file, rotamer_recovery )

void
register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	OPT( in::file::l );
	OPT( in::file::silent );
	OPT( in::file::native );
	OPT( packing::use_input_sc );
	OPT( packing::resfile );
	NEW_OPT( rotamer_recovery::reporter, "Rotamer Recovery Reporter component.", "");
	//	NEW_OPT( rotamer_recovery::native_vs_decoy, "Table of describing which native structures the decoys are modeling.  The first column is native pdb filename and the second column is the decoy pdb filename.", "");
	NEW_OPT( rotamer_recovery::comparer, "Rotamer Recovery Comparer component.", "");
	NEW_OPT( out::file::rotamer_recovery, "Output File  Name for reporters that write their results to a file.", "");

}

int
main( int argc, char * argv [] )
{
  register_options();

  init(argc, argv);

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace basic::options::OptionKeys::out::file;

	string reporter( option[ rotamer_recovery::reporter ].value() );
	string output_fname( option[ out::file::rotamer_recovery ].value() );
	string comparer( option[ rotamer_recovery::comparer ].value() );
	ScoreFunctionOP scfxn( getScoreFunction() );
	TaskFactoryOP task_factory( new TaskFactory );
	task_factory->push_back( new InitializeFromCommandline );
	if(option.has(OptionKeys::packing::resfile) && option[OptionKeys::packing::resfile].user()){
		task_factory->push_back( new ReadResfile );
	}
	task_factory->push_back( new RestrictToRepacking );

  RotamerRecoveryMoverOP rr(
		new RotamerRecoveryMover(
			reporter,
			output_fname,
			comparer,
			scfxn,
			task_factory
		)
	);


  try{
    JobDistributor::get_instance()->go( rr );
		rr->show();
  } catch (EXCN_Base & excn ) {
    TR.Error << "Exception: " << endl;
    excn.show( TR.Error );

    TR << "Exception: " << endl;
    excn.show( TR ); //so its also seen in a >LOG file
  }

  return 0;

}








