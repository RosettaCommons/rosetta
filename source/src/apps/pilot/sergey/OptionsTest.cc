// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief
/// @author Sergey Lyskov

#include <core/init/init.hh>
#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>
#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <basic/options/option.hh>

#include <numeric/random/random.hh>


#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/database/open.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <ostream>
#include <fstream>

//#include <boost/libs/thread/thread.cpp>
//#include <boost/libs/thread/tss.cpp>
//#include <boost/libs/thread/mutex.cpp>
//#include <boost/libs/thread/exceptions.cpp>
//#include <boost/libs/thread/once.cpp>

#include <iostream>


// option key includes

#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers


static basic::Tracer TM( "TMemory" );

static basic::Tracer TR_( "global" );


void TracerDiskSpaceTest(void)
{
	for ( int i=0; i<500000; ++i ) {
		TM << "Writing line: " << i << " _______________________________________________________________________________________"<< std::endl;
	}

	TM << "Normal exit..." << std::endl;
}


class OutputMover : public protocols::moves::Mover
{
public:
	void apply( core::pose::Pose & )
	{
		TracerDiskSpaceTest();
	}

	std::string get_name() const { return ""; }
};

typedef utility::pointer::owning_ptr< OutputMover > OutputMoverOP;
typedef utility::pointer::owning_ptr< OutputMover const > OutputMoverCOP;


using basic::Error;
using basic::Warning;

//basic::Tracer TR("TTest", basic::t_info, true);
//basic::Tracer TR("core.TTest");
//basic::Tracer TR2("core.TTest.T2");


using namespace core;
using namespace basic;


#include <utility/tools/make_vector.hh>

int main( int argc, char * argv [] )
{
	try {

		basic::Tracer TR( "main" );

		using namespace core;
		using namespace basic::options::OptionKeys;

		basic::options::option.add_relevant(in::path::database);

		devel::init(argc, argv);

		{
			using namespace basic::options;
			using namespace basic::options::OptionKeys;


			if ( option[ out::levels ].active() ) T("Levels:") << option[ out::levels ]() <<  std::endl;

			//option[ out::levels ]( utility::tools::make_vector(std::string("some.namespace:info"), std::string("some.other.namespace:debug")) );

			option[ out::levels ]("some.namespace:info");
			option[ out::levels ]("some.other.namespace:debug");

			if ( option[ out::levels ].active() ) T("Levels now:") << option[ out::levels ]() <<  std::endl;

			{
				core::pose::Pose pose;
				core::import_pose::pose_from_file(pose, "test_in.pdb", core::import_pose::PDB_file);

				core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );
				T("Score:") << scorefxn->score(pose)  << std::endl;
				pose.energies().residue_total_energies(1);
				T("Scoring done!") << "---------------------" << std::endl;
			}
			{
				T("Testing pose_from_sequence...") << std::endl;
				std::string sequence(1000, 'V');
				core::pose::Pose pose;
				core::pose::make_pose_from_sequence(pose, sequence, core::chemical::FA_STANDARD );
				T("Testing pose_from_sequence... Done!") << std::endl;
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
