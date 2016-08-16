// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/delucasl/roc_test.cc
/// @author Sam DeLuca

#include <iostream>

#include <devel/init.hh>
#include <numeric/roc_curve.hh>


int main(int argc, char* argv[])
{
    try {
	devel::init(argc,argv);

	numeric::RocCurve roc;

	roc.insert_point(false,true,"a",0.2);
	roc.insert_point(true,true,"b",0.1);
	roc.insert_point(true,false,"c",0.0);

	roc.generate_roc_curve();
	roc.print_roc_curve();
	std::cout <<roc.calculate_auc() <<std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
       return 0;
}
