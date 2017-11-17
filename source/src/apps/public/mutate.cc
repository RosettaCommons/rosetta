// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Mutate a residue
/// @author JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/conformation/Residue.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

// utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/mutate.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::simple_moves;
using namespace utility;
using namespace core::scoring;
using namespace core::pack::task;

static basic::Tracer TR( "apps.public.mutate" );

////////////////////////////////////////////////////////////////////////////////

Pose read_pose() {

	// cry if PDB not given
	if ( ! option[OptionKeys::in::file::s].user() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide PDB file!");
	}

	// read in pose
	Pose pose;
	core::import_pose::pose_from_file( pose, option[OptionKeys::in::file::s].value_string() , core::import_pose::PDB_file);
	TR.Debug << "got pose of length " << pose.size() << std::endl;

	return pose;

}// read pose

////////////////////////////////////////////////////////////////////////////////

/// @brief Convert one to three letter code
std::string one2three( std::string one ) {

	// error checking
	if ( one == "B" || one == "J" || one == "O" || one == "U" || one == "X" || one == "Z" ) {
		utility_exit_with_message( "One letter code doesn't belong to any of the 20 natural amino acids!" );
	}

	// create the vectors
	utility::vector1< std::string > olc;
	olc.push_back( "A" );
	olc.push_back( "C" );
	olc.push_back( "D" );
	olc.push_back( "E" );
	olc.push_back( "F" );
	olc.push_back( "G" );
	olc.push_back( "H" );
	olc.push_back( "I" );
	olc.push_back( "K" );
	olc.push_back( "L" );
	olc.push_back( "M" );
	olc.push_back( "N" );
	olc.push_back( "P" );
	olc.push_back( "Q" );
	olc.push_back( "R" );
	olc.push_back( "S" );
	olc.push_back( "T" );
	olc.push_back( "V" );
	olc.push_back( "W" );
	olc.push_back( "Y" );

	utility::vector1< std::string > tlc;
	tlc.push_back( "ALA" );
	tlc.push_back( "CYS" );
	tlc.push_back( "ASP" );
	tlc.push_back( "GLU" );
	tlc.push_back( "PHE" );
	tlc.push_back( "GLY" );
	tlc.push_back( "HIS" );
	tlc.push_back( "ILE" );
	tlc.push_back( "LYS" );
	tlc.push_back( "LEU" );
	tlc.push_back( "MET" );
	tlc.push_back( "ASN" );
	tlc.push_back( "PRO" );
	tlc.push_back( "GLN" );
	tlc.push_back( "ARG" );
	tlc.push_back( "SER" );
	tlc.push_back( "THR" );
	tlc.push_back( "VAL" );
	tlc.push_back( "TRP" );
	tlc.push_back( "TYR" );

	// do the conversion
	for ( Size i = 1; i <= 20; ++i ) {
		if ( olc[ i ] == one ) {
			return tlc[ i ];
		}
	}

	return "this is wrong";

} // one to three letter code


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		// initialize option system, RNG, and all factory-registrators
		devel::init(argc, argv);

		// read spanfile
		Pose pose = read_pose();

		// read option for mutation
		// get mutation
		if ( ! option[ OptionKeys::mutate::mutation ].user() ) {

			utility_exit_with_message("No mutations given. Use flag -mutate:mutation. Quitting.");
		} else {
			// input format A163F, can give multiple mutations
			std::string mutation = option[ OptionKeys::mutate::mutation ]();
			TR << "Looking at mutation " << mutation << std::endl;

			// initialize single point variables
			std::string wt, mut;
			Size seqid;

			// get wt and mutant from vector: split string by character
			wt = mutation[ 0 ];
			mut = mutation[ mutation.size()-1 ];

			// get residue number
			utility::vector1< std::string > tmp;
			for ( Size i = 1; i <= mutation.size()-2; ++i ) {
				tmp.push_back( to_string( mutation[ i ] ) );
			}

			seqid = string2Size( join( tmp, "" ) );
			TR << "wt " << wt << ", seqid " << seqid << ", mut " << mut << std::endl;

			// make mutation
			MutateResidueOP mutate( new MutateResidue( seqid, one2three( mut ) ) );
			mutate->apply( pose );

			// create scorefunction or get full-atom one
			//   ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "talaris2014.wts" );
			ScoreFunctionOP sfxn = get_score_function( true );

			// repack the mutated residue
			utility::vector1< bool > repack_residues( pose.total_residue(), false );
			repack_residues[ seqid ] = true;
			PackerTaskOP repack = TaskFactory::create_packer_task( pose );
			//   repack->pack_residue( static_cast< int >( seqid ) );
			repack->restrict_to_repacking();
			repack->restrict_to_residues( repack_residues );
			core::pack::pack_rotamers( pose, *sfxn, repack );

			// output name
			std::string protein = option[ OptionKeys::in::file::s ](1);
			const std::string tmp1( file_basename( protein ) );
			std::string output = trim( tmp1, ".pdb");

			// add mutation to output filename
			output += "_" + wt + to_string( seqid ) + mut + ".pdb";

			// dump pdb
			pose.dump_pdb( output );
		}
	}
catch (utility::excn::Exception const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
