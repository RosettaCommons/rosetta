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
/// @brief  Demo for Tracer Color output
/// @author Sergey Lyskov

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tools/make_vector.hh>

#include <iostream>

#include <unistd.h>

basic::Tracer TR("ColorDemo");

int main( int argc, char * argv [] )
{
	try {
		using namespace basic;

		devel::init(argc, argv);

		TR << "isatty:" << isatty(fileno(stdout)) << '\n' << std::endl << std::endl;

		TR << TR.bgWhite << TR.Black << "Example" << TR.Reset << std::endl;
		TR << TR.Red << "Red text " << TR.Bold << "Red-bold-or-bright " << TR.Reset << TR.Black << TR.bgRed << "Black-on-red" << TR.Reset << std::endl;
		TR << TR.Green << TR.Underline << " Green with underline\n" << TR.Reset << std::endl;

		std::vector< basic::CSI_Sequence > all_codes = utility::tools::make_vector< basic::CSI_Sequence > (
			Tracer::Black, Tracer::Red, Tracer::Green, Tracer::Yellow, Tracer::Blue, Tracer::Magenta, Tracer::Cyan, Tracer::White,
			Tracer::bgBlack, Tracer::bgRed, Tracer::bgGreen, Tracer::bgYellow, Tracer::bgBlue, Tracer::bgMagenta, Tracer::bgCyan, Tracer::bgWhite,
			Tracer::Bold, Tracer::Underline);

		for(int i=0; i<all_codes.size(); ++i) {
			for(int j=0; j<all_codes.size(); ++j) TR << all_codes[i] << all_codes[j] << " X " << TR.Reset;
			TR << std::endl;
		}

		TR << "Regular text...\n" << std::endl;

		utility_exit_with_message("\nExample of utility_exit_with_message...");
		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}
