// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/scenarios/FloppyTail.cc
/// @brief This app was initially intended for modeling the binding of a long unstructured C-terminal tail to some other part of a protein.  It now works for N-terminal, C-terminal, and internal flexible regions.  It works best as a method for sampling space to see what is possible, preferably in conjunction with extensive experimental constraints.  It is not meant to produce ab-initio style models of folded complexes.
/// @author Steven Lewis

// Unit Headers
#include <protocols/floppy_tail/FloppyTailMover.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

static thread_local basic::Tracer TR( "apps.public.scenarios.FloppyTail" );

int main( int argc, char* argv[] )
{
	try {
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new protocols::floppy_tail::FloppyTailMover ));

		TR << "************************d**o**n**e**************************************" << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
