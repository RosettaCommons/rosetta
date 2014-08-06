// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/utility/true_false.cc
/// @brief  demo of True/False query class which handles basic boolean logic
/// @author Will Sheffler (willsheffler@gmail.com)
/// @date   Thu Aug  9 19:49:23 2007
///

#include <cstdlib>
#include <iostream>
#include <vector>

#include <utility/query/types.hh>
#include <utility/query/Metric.hh>
#include <utility/query/Filter.hh>

#include <utility/pointer/owning_ptr.hh>
#include <ObjexxFCL/string.functions.hh>


using utility::pointer::owning_ptr;
using namespace utility::query;
using namespace std;
using ObjexxFCL::string_of;

void print_passed( Filter<int> cond, int min=0, int max=15 )
{
	std::cerr << "Condition: " << std::endl << cond.description() << std::endl << std::endl;
	for( int i = min; i < max; ++i)
		if( cond( i ) ) std::cerr << "condition holds for " << i << std::endl;
	std::cerr << std::endl;
}

void test_int() {
	std::cerr << "============== int tests =============" << std::endl;

	Metric0Param<int> I;

	print_passed( I==8 | I==7 );
	print_passed( 3 < I & I < 7 );

}


#include <Vec.hh>

typedef utility::pointer::owning_ptr< FilterBase<Vec> > TF_Vec_OP;
Filter<Vec> __cplusplus_is_stupid__;

Vec rand_Vec_condition( Filter<Vec> cond, int len=-1, int int_max=10, int MAX_TRY=9999 ) {
	for( int it = 0; it < MAX_TRY; ++it) {
		Vec v = rand_Vec(len,int_max); // from Vec.hh
		if( cond(v) ) return v;
	}
	std::cerr << "COULDN'T SATISFY CONDITION: " << std::endl;
	std::cerr << cond.description() << std::endl;
	return Vec();
}

void test_Vec() {
	std::cerr << "============== Vec tests =============" << std::endl;

	Filter1ParamGenerator<Vec,int> contains( new VecContains );

	Metric0Param<Vec>  min( new VecMinElement );
	Metric0Param<Vec>  max( new VecMaxElement );
	Metric0Param<Vec>  len( new VecLen );
	Metric0Param<Vec> mean( new VecMean );

	Filter<Vec> cond = contains(7) & ( len < 3 | 5 < len & mean > 3.7 & max < 8 & min >= 2);
	std::cerr << "Condition: " << std::endl << cond.description() << std::endl << std::endl;

	for(int i=0; i<10; ++i) {
		std::cerr <<  rand_Vec_condition( cond ) << std::endl;
	}
}



int main (int argc, char const* argv[])
{
	argc++;
	char const** shutup = argv;
	shutup++;

	test_int();

	test_Vec();

	std::cerr << "DONE!" << std::endl;

	return 0;

}














