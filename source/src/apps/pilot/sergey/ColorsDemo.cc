// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief  Demo for Tracer Color output
/// @author Sergey Lyskov

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tools/make_vector.hh>
#include <utility/CSI_Sequence.hh>

#include <iostream>

#include <unistd.h>

static basic::Tracer TR( "ColorDemo" );

using namespace utility;

static basic::Tracer     Blue("blue",       CSI_Blue());
static basic::Tracer    Green("green",      CSI_Green());
static basic::Tracer      Red("red",        CSI_Red());
static basic::Tracer VeryBlue("very.blue",  CSI_Blue()+CSI_Bold(), CSI_Blue());


int main( int argc, char * argv [] )
{
	try {
		using namespace basic;

		devel::init(argc, argv);

		TR << "isatty:" << isatty(fileno(stdout)) << '\n' << std::endl << std::endl;

		TR << TR.bgWhite << TR.Black << "Example" << TR.Reset << std::endl;
		TR << TR.Red << "Red text " << TR.Bold << "Red-bold-or-bright " << TR.Reset << TR.Black << TR.bgRed << "Black-on-red" << TR.Reset << std::endl;
		TR << TR.Green << TR.Underline << "Green with underline\n" << TR.Reset << std::endl;

		std::vector< CSI_Sequence > all_codes = utility::tools::make_vector< CSI_Sequence > (
			Tracer::Black, Tracer::Red, Tracer::Green, Tracer::Yellow, Tracer::Blue, Tracer::Magenta, Tracer::Cyan, Tracer::White,
			Tracer::bgBlack, Tracer::bgRed, Tracer::bgGreen, Tracer::bgYellow, Tracer::bgBlue, Tracer::bgMagenta, Tracer::bgCyan, Tracer::bgWhite,
			Tracer::Bold, Tracer::Underline);

		for ( unsigned int i=0; i<all_codes.size(); ++i ) {
			for ( unsigned int j=0; j<all_codes.size(); ++j ) TR << all_codes[i] << all_codes[j] << " X " << TR.Reset;
			TR << std::endl;
		}

		TR << "Regular text...\n" << std::endl;

		Blue << "Output with default color for Blue Tracer..." << std::endl;
		Green << "Output with default color for Green Tracer..." << std::endl;
		Red << "Output with default color for Red Tracer...\n" << std::endl;

		VeryBlue << "Output with default color for VeryBlue Tracer..." << std::endl;
		VeryBlue << CSI_Green() << "We" << CSI_Yellow() << " can still use " << CSI_Red() << "other" << CSI_Bold() << " colors in Tracers " << CSI_bgBlue() << CSI_Yellow() << "with default color set!" << std::endl;
		VeryBlue << "..." << std::endl;
		VeryBlue << CSI_Reset() << "After reset on VeryBlue Tracer...\n" << std::endl;

		TR.Warning << "Example of Warning on default tracer" << std::endl;
		TR.Error   << "Example of Error on default tracer" << std::endl;
		TR.Fatal   << "Example of Fatal on default tracer" << std::endl;

		VeryBlue.Warning << "Example of Warning on VeryBlue tracer" << std::endl;
		VeryBlue.Error   << "Example of Error on VeryBlue tracer" << std::endl;
		VeryBlue.Fatal   << "Example of Fatal on VeryBlue tracer" << std::endl;

		std::string empty_string("");

		//assert( "Plain Assert Example" == empty_string );
		//try {
		// debug_assert( "Debug Assert Example" == empty_string );
		//} catch (utility::excn::Exception const & e ) {}
		//try {
		// runtime_assert( "Runtime Assert Example" == empty_string );
		//} catch (utility::excn::Exception const & e ) {}

		utility_exit_with_message("\nExample of utility_exit_with_message...");
		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
