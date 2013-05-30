// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   scratch1.cc
/// @brief  Rosetta exploration code
/// @author Leif Arthur Halvorsen (lah435@nyu.edu)
//
//
//////////::::
//#include <protocols/jobdist/standard_mains.hh>
//#include <core/scoring/constraints/util.hh>
//#include <basic/options/option.hh>
//#include <basic/options/after_opts.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/pose/Pose.fwd.hh>
//
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/init/init.hh>
//#include <devel/init.hh>
//#include <protocols/relax/RelaxProtocolBase.hh>
//#include <protocols/relax/ClassicRelax.hh>
//#include <protocols/relax/util.hh>
//
//// C++ headers
//#include <cstdlib>
//#include <fstream>
//#include <iostream>
//#include <string>
//
////silly using/typedef
//
//#include <protocols/jd2/JobDistributor.hh>
//#include <core/init/init.hh>
//
//#include <protocols/moves/PackRotamersMover.hh>
//#include <protocols/moves/MoverContainer.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//
//#include <core/pack/task/operation/TaskOperations.hh>
//#include <core/pack/task/TaskFactory.hh>
//
//#include <basic/Tracer.hh>
//
//// option key includes
//
//#include <basic/options/keys/relax.OptionKeys.gen.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/packing.OptionKeys.gen.hh>


////////::::

// Project headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

// Prococol headers

// C++ headers
#include <iostream>
#include <string>

// Options headers
//#include <basic/options/option.hh>
//
//// Utility headers
//#include <basic/Tracer.hh>
//#include <utility/vector1.hh>

//static basic::Tracer TR("MyRosetta");


//using utility::vector1;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

int
//main (int argc, const char * argv[])
main (int argc, char * argv[])
{

	//devel::init(argc, argv);
	//TR << "This is fucking tracer output!!!" << std::endl;
	
	std::cout << "Hello, World!\n" << std::endl;
    return 0;
}



















//// Project headers
//
//// Core headers
//#include <core/pose/Pose.hh>
//#include <core/init/init.hh>
//
//// Prococol headers
//
//// C++ headers
//#include <iostream>
//#include <string>
//
//// Options headers
////#include <basic/options/option.hh>
//
//// Utility headers
//#include <basic/Tracer.hh>
//#include <utility/vector1.hh>
//
//static basic::Tracer TR("MyRosetta");
//
//
//using utility::vector1;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

//main (int argc, char * argv[])
//{
//	// initialize the options system
//	core::init::init(argc, argv);
//
//	TR << "This is fucking tracer output!!!" << std::endl;
//	
//	// Make a vector of ints and add some numbers
//	utility::vector1< int > my_vector;
//	my_vector.push_back(3);
//	my_vector.push_back(2);
//	my_vector.push_back(1);
//	
//	// Now let's print those numbers
//	std::cout << my_vector.front() << std::endl;
//	std::cout << my_vector.at(2) << std::endl;
//	std::cout << my_vector.back() << std::endl;
//	
//	// Can use overloaded operator instead of 'at' accessor
//	std::cout << my_vector[2] << std::endl;
//	std::cout << my_vector[4] << std::endl;
//	
//	// Get the capacity
//	std::cout << "Capacity: " << my_vector.size() << std::endl;
//	
//	
//	
//	// Get an instance of the job distributor and pass it 
//	
//	std::cout << "Hello, World!\n" << std::endl;	
//    return 0;
//}
