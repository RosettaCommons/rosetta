// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <devel/init.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "apps.pilot.kale.examples.Tracer" );

int main(int argc, char* argv[]) {

	// This function is necessary.  The tracer will not work until the rosetta 
	// libraries have been initialized.
	
	devel::init(argc, argv);

	TR << "Use the '-out:level' argument to control how " << std::endl;
	TR << "the following messages are displayed:\n" << std::endl;

	TR.Trace   << "Trace:   Priority 500" << std::endl;
	TR.Debug   << "Debug:   Priority 400" << std::endl;
	TR.Info    << "Info:    Priority 300" << std::endl;
	TR.Warning << "Warning: Priority 200" << std::endl;
	TR.Error   << "Error:   Priority 100" << std::endl;
	TR.Fatal   << "Fatal:   Priority 0"   << std::endl;
}


