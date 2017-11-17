// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/fragments_to_ss
/// @brief  outputs the fragments as if they were a psipred secondary structure prediction. This will be used to calibrate
//          rosetta abinitio
/// @author TJ Brunette

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FragData.hh>

#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <ObjexxFCL/format.hh>

#include <utility/io/ozstream.hh>

class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}
OPT_KEY( File, fragment_file_in )
OPT_KEY( File, fragment_file_out )

using namespace core;
using utility::vector1;

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( fragment_file_in, "fragment file in", "" );
  NEW_OPT( fragment_file_out, "fragment file out", "" );
}

int main( int argc, char * argv [] ) {
	try {
		using namespace basic;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace fragment;
		using namespace ObjexxFCL::format;
		ThisApplication::register_options();
		devel::init(argc, argv);
    FragSetOP orig_frags = FragmentIO().read_data( option[ fragment_file_in ]() );
		FragmentIO().write_data(option[fragment_file_out](),*orig_frags);
 	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

