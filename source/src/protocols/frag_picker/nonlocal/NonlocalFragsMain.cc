// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/frag_picker/nonlocal/NonlocalFragsMain.cc
/// @author David Kim

/*
The purpose of this protocol is to get interacting fragment pairs non-local in sequence from
a given PDB. This protocol will be used to generate a VALL database of non-local fragment pairs.
*/


// Project headers
#include <protocols/frag_picker/nonlocal/NonlocalFragsMain.hh>
#include <protocols/frag_picker/nonlocal/NonlocalFrags.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
//#include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// C/C++ headers
#include <iostream>
//#include <string>

namespace protocols  {
namespace frag_picker {
namespace nonlocal {

static basic::Tracer TR( "protocols.frag_picker.nonlocal.NonlocalFragsMain" );

using namespace std;

void NonlocalFrags_main() {
	//using namespace basic::options;
	//using namespace basic::options::OptionKeys;
	using protocols::jd2::JobDistributor;

	// necessary for outputting intermediate structures
	//option[ OptionKeys::run::intermediate_structures ].value(true);

	NonlocalFragsOP mover;

	mover = NonlocalFragsOP( new NonlocalFrags() );

	try {
		JobDistributor::get_instance()->go(mover);
	} catch (utility::excn::Exception& e) {
		cerr << "Exception: " << endl;
		e.show(cerr);
	}

}

}  // namespace nonlocal
}  // namespace frag_picker
}  // namespace protocols
