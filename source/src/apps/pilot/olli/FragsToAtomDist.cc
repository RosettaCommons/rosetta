// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#include <string>
#include <core/types.hh>
#include <ostream>
#include <stdlib.h>
#include <devel/init.hh>
#include <utility>
#include <protocols/noesy_assign/FragsToAtomDist.fwd.hh>
#include <protocols/noesy_assign/FragsToAtomDist.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <utility/excn/EXCN_Base.hh>

static basic::Tracer tr( "main" );

OPT_KEY( Integer, cycles )
OPT_KEY( Integer, freq )
OPT_KEY( File, fragments )
OPT_KEY( Boolean, no_r6_averaging )
OPT_KEY( File, distances )

int main( int argc, char * argv [])
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core;
		using namespace std;
		using namespace protocols;
		OPT( in::file::fasta );
		OPT( score::weights );
		OPT( score::patch );
		NEW_OPT( cycles, "number of cycles", 1000 );
		NEW_OPT( freq, "calculate distance every nth cycle", 10 );
		NEW_OPT( fragments, "input fragment file", "frags.dat.gz" );
		//NEW_OPT( no_r6_averaging, "whether use r6 average or not", false );
		NEW_OPT( distances, "distance database file", "dist.dat" );

		devel::init( argc, argv );
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		std::string const full_sequence = fasta_sequence->sequence();
		protocols::noesy_assign::FragsToAtomDist Atomdist1;
		Atomdist1.generate_from_frag_file( option[ fragments ](), full_sequence, option[ cycles ](), option[ freq ]() );
		Atomdist1.write_hist_to_file( option[ OptionKeys::distances ]() );
	} catch ( utility::excn::EXCN_Base& excn ) {
		excn.show( std::cerr );
		std::exit( 1 );
	}

	return 0;
}
