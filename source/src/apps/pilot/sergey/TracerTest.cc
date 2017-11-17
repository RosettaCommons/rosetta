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

#include <devel/init.hh>
#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>
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

// Utility headers
#include <utility/exit.hh>


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

void show(std::ostream & tr)
{
	tr << "Something \n";
	tr << 1 << 2 << 3 << "\n";
}

class AA {};

//typedef std::ostream& ENDL ( std::ostream& os );

std::ostream& operator <<(std::ostream &tr, AA)
{
	tr << "Class A \n";
	return tr;
}


void test_Tracer()
{
	basic::Tracer TR( "core.TTest" );
	basic::Tracer TR2( "core.TTest.T2" );

	TR.Fatal   << "Some Fatal Error here... " << "000 " << std::endl;
	TR.Fatal   << "2222 Fatal Error here... " << "___ " << std::endl;
	TR << "... regular output ";

	TR.Error   << "Some Error here... " << "100" << std::endl;
	TR << "... regular output " << std::endl;
	TR.Warning << "Warning... " << "200" << std::endl;
	TR << "... regular output " << std::endl;
	TR.Info    << "Info... " << "300" << std::endl;
	TR << "... regular output " << std::endl;
	TR.Debug   << "Debug... " << "400" << std::endl;
	TR << "... regular output " << std::endl;
	TR.Trace   << "Trace... " << "500" << std::endl;
	TR << "... regular output " << std::endl;

	TR2 << "TR2 regular" << std::endl;
	TR2.Debug << "TR2.Debug " << std::endl;
	TR2.Trace << "TR2.Trace " << std::endl;


	TR.Fatal   << "Setting default level to 500 -----------------------------------" << std::endl;
	TR(300) << "500!" << std::endl;
	TR.Fatal   << "Some Fatal Error here... unflushed...\n";// << std::endl;
	TR << "... regular output " << std::endl;
	TR.Error   << "Some Error here..." << std::endl;
	//TR.flush();

	TR << "... regular output \n" << std::endl;
	TR.Warning << "Warning...\n" << std::endl;
	TR << "... regular output \n" << std::endl;
	TR.Info    << "Info...\n" << std::endl;
	TR << "... regular output \n" << std::endl;
	TR.Debug   << "Debug...\n" << std::endl;
	TR << "... regular output \n" << std::endl;
	TR.Trace   << "Trace... unflushed...\n";// << std::endl;
	TR.Fatal   << "Fatal... unflushed...\n";// << std::endl;
	TR << "... regular output \n" << std::endl;

	show( TR );
	//show( std::cout );

	AA a;
	TR << a;

	//OStreamTracer(std::cout) << a;

	//std::fstream file1("123");
	//std::fstream f2(file1);

	std::ostringstream osstr;
	//std::ostringstream osstr2 = osstr;
	osstr << std::endl;

	//TR << "123..." << std::endl;


	//show( OStreamTracer(std::cout) );

	//TR << "DB:"  << basic::options::option[ in::path::database ]() << "\n";

	//chemical::ResidueTypeSetCAP residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );


	//TR << "Specified()=" << basic::options::option[ run::benchmark_scale ].specified() << "\n";
	//TR << "Legal" << basic::options::option[ run::benchmark_scale ].legal() << "\n";
	//TR << "Active:" << basic::options::option[ run::benchmark_scale ].active() << "\n";
	//TR << "native:"  << basic::options::option[ james::native ]() << "\n";
	//TR << "DB:"  << basic::options::option[ in::path::database ]() << "\n";

	//int scale = basic::options::option[ run::benchmark_scale ]();

	//TR << "Mini Benchmark started! Scale factor: " << scale << " -------------\n";

}


void test_Random(void)
{
	//devel::init_random_generators(1000, numeric::random::_RND_TestRun_, "ran3");
	devel::init_random_generators(1000, numeric::random::_RND_TestRun_, "mt19937");
	for ( int i=0; i<100; i++ ) {
		double r = numeric::random::uniform();
		TR_.precision(25);
		TR_ << i << " " << r << std::endl;
	}
}


int main( int argc, char * argv [] )
{

	try {

		basic::Tracer TR( "main" );

		using namespace core;
		using namespace basic::options::OptionKeys;

		basic::options::option.add_relevant(in::path::database);

		devel::init(argc, argv);

		{
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

		test_Tracer();
		std::cout << "Done !-------------------------------" << std::endl;
		//return 0;


		TR << "Some unflushed output 1...\n";
		TR << "Some unflushed output 2...\n";

		//std::cout << "Calling basic::Tracer::flush_all_tracers manually..." << std::endl;
		//basic::Tracer::flush_all_tracers();
		//std::cout << "Calling basic::Tracer::flush_all_tracers manually... Done." << std::endl;
		utility_exit_with_message("\nExiting without flushing the tracers...");


		// Testing std::IO errors
		std::cerr << "Point 1" << std::endl;
		OutputMoverOP om = new OutputMover();
		protocols::jd2::JobDistributor::get_instance()->go(om);
		std::cerr << "Point 2" << std::endl;


		std::cout << "-------------------------------" << std::endl;
		{
			//#boost::thread thrd(&hello);

			//boost::thread_specific_ptr< int > * mint;
			//mint = new boost::thread_specific_ptr< int >;


		}
		std::cout << "Done !-------------------------------" << std::endl;
		return 0;


		//TR << "clock:" << clock() << "\n";


		test_Tracer();
		//test_Random();

		TR << "TTest ended. #@!  --------------------------------" << std::endl;
		//TR.flush();
		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
