// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/brunette/cluster_alns.cc
///
/// @brief  Divide input alns into clusters based on gdtmm comparison of partial models.
/// @author TJ Brunette


#include <devel/init.hh>
#include <protocols/comparative_modeling/AlignmentClustering.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


int main( int argc, char * argv [] ) {
	try {

		using namespace protocols::comparative_modeling;
		devel::init(argc, argv);
		AlignmentClusteringOP cluster( new AlignmentClustering() );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
