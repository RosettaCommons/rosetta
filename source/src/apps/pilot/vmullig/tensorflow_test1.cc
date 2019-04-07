// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file tensorflow_test1.cc
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief Copied verbatim from the Tensorflow website.  This is intended to be
/// a quick test of whether we're properly linking the Tensorflow libraries.

#include <tensorflow/c/c_api.h>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

static basic::Tracer TR( "apps.pilot.vmullig.tensorflow_test1" );

int
main( int argc, char * argv [] ) {
	try {
		devel::init(argc, argv);
#ifdef USE_TENSORFLOW
		TR << "Hello from TensorFlow C library version " << TF_Version() << std::endl;
#else
		utility_exit_with_message("This application requires tensorflow!\n");
#endif
	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}
