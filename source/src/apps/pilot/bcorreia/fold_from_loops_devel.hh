// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers

#include <core/types.hh>
#include <utility/vector1.hh> // there is no forward declaration possible for const_iterator(?)
#include <basic/options/option.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <string>
#include <list>
#include <iosfwd>
#include <iostream>
#include <sstream>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
static thread_local basic::Tracer TRA( "fold_from_loops", basic::t_info );



void write_checkpoint( Size iter )
{
	if ( ! option[ OptionKeys::fold_from_loops::checkpoint ].user() ) return;
	std::string fileroot( option[ OptionKeys::fold_from_loops::checkpoint ]() );

	TRA << "writing  iteration checkpoint " << '\n';

	// write checkpoint file
	std::string checkpointname( fileroot + ".checkpoint" );
	utility::io::ozstream out( checkpointname.c_str() );

	if ( !out ) {
		std::cerr << "trouble opening file " << checkpointname
		          << " for writing... skipping checkpoint" << std::endl;
		runtime_assert( false ); // die here in debug mode
		return;
	}

	// here iter should refer to the last complete iteration
	out << "Iteration " << iter << '\n';
	out.close();
}


void
load_checkpoint( Size & iter )
{
	if ( ! option[ OptionKeys::fold_from_loops::checkpoint ].user() ) return;
	std::string fileroot( option[ OptionKeys::fold_from_loops::checkpoint ]() );

	utility::io::izstream file;
	std::string filename( fileroot + ".checkpoint" );
	file.open( filename.c_str() );
	if ( !file ) return;

	TRA << "Reading nstruct checkpoint info from " << filename << '\n';

	std::string line, word ;
	// get iteration
	Size last_iter;
	file >> word >> last_iter >> utility::io::skip; // first line
	if ( ( word != "Iteration" ) ) return;
	file.close();


	// here iter should refer to the last complete iteration
	iter = last_iter + 1 ;

}
