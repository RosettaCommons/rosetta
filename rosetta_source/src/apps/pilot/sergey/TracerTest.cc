// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///
/// @brief
/// @author Sergey Lyskov

#include <core/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>
#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <basic/options/option.hh>

#include <numeric/random/random.hh>



#include <ostream>
#include <fstream>


basic::Tracer TM("TTest");

void TracerDiskSpaceTest(void)
{
	for(int i=0; i<500000; ++i) {
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


using basic::T;
using basic::Error;
using basic::Warning;

basic::Tracer TR("TTest", basic::t_info, true);

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
	//core::init_random_generators(1000, numeric::random::_RND_TestRun_, "ran3");
	core::init_random_generators(1000, numeric::random::_RND_TestRun_, "mt19937");
	for(int i=0; i<100; i++) {
		double r = numeric::random::uniform();
		TR.precision(25);
		TR << i << " " << r << std::endl;
	}
}




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



int main( int argc, char * argv [] )
{
	using namespace core;
	using namespace basic::options::OptionKeys;

	basic::options::option.add_relevant(run::benchmark_scale);
	basic::options::option.add_relevant(in::path::database);

	devel::init(argc, argv);

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
}

