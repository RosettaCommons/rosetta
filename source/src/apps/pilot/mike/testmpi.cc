// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


#include <protocols/wum/MPI_Relax.hh>
#include <cstdio>

#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/FArray4D.hh>

#include <unistd.h>
#include <ctime>
#include <sys/time.h>


int
main( int argc, char * argv [] )
{
    try {
	// initialize core
	devel::init(argc, argv);
	//protocols::wum::MPI_Relax wum;
	//wum.go();

	//ObjexxFCL::FArray4D< double > monster_table( 100, 100, 100, 10 );


	// Unique string ident:

//	core::Real cl = clock()/CLOCKS_PER_SEC;
//	long       ti = time(NULL);
//	timeval a;
//  gettimeofday(&a, 0);
//	timespec tp;
//	clock_gettime(CLOCK_REALTIME, &tp);
//	std::cout << cl << "  " << ti << a.tv_sec << " " << a.tv_usec << "  " << tp.tv_sec << "  " << tp.tv_nsec << std::endl;

	exit(0);


	utility::vector1< float > myvec;

	std::cout << myvec.capacity() << std::endl;

	for(int i =0; i < 1500; i++){
		myvec.resize(i);
		myvec.reserve( myvec.size()+1 );
		std::cout << myvec.size() << "  " << myvec.capacity() << std::endl;

	}


///	std::string test;
///
///	test = "LLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdb";
///
///	std::cout << test.capacity() << std::endl;
///
///	test.clear();
///
///	std::cout << test.capacity() << std::endl;
///
///	test.reserve(2);
///
///	std::cout << test.capacity() << std::endl;
///
///	test = "LLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAdb";
///
///	std::cout << test.capacity() << std::endl;
///
///	core::Real a=1e-20;
///	for( int inf=0; inf < 1; inf=inf) a = a*1.00000001;
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}


