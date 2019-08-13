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
//#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>

// utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/mutate.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/relax/membrane/util.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

#include <utility/io/util.hh>

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
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		using namespace utility;
		using namespace utility::io;
		using namespace protocols::relax::membrane;

		// initialize option system, RNG, and all factory-registrators
		devel::init(argc, argv);

		// read pose
		Pose input_pose = read_pose();

		// simple input checking
		if ( ! option[ OptionKeys::mutate::mutation ].user()
				&& ! option[ OptionKeys::mutate::mutant_file ].user() ) {
			utility_exit_with_message("No mutations given. Use flag -mutate:mutation or -mutate:mutant_file. Quitting.");
		}

		if ( option[ OptionKeys::mutate::mutant_file ].user() &&
				option[ OptionKeys::mutate::mutation ].user() ) {

			utility_exit_with_message( "Too many inputs: You must EITHER specify the mutant file with -mutate:mutant_file OR specify a single mutation with -mutate:mutation! Quitting..." );
		}

		// one line per input file: several mutants within a single construct
		utility::vector1< std::string > mutations;

		// initialize vectors for wildtype residue, mutant residue, and sequence ID
		// Outer vector is different constructs, inner vector
		//   is multiple residues per construct
		utility::vector1< utility::vector1< char > > wt_res, mut_res;
		utility::vector1< utility::vector1< core::Size > > resid;

		// mutants given in command line
		// input format A163F in pose numbering, can give multiple mutations
		// split multiple mutations by whitespace
		if ( option[ OptionKeys::mutate::mutation ].user() ) {
			mutations = option[ OptionKeys::mutate::mutation ]();
		}
		// mutants given in mutant file
		if ( option[ OptionKeys::mutate::mutant_file ].user() ) {
			mutations = get_lines_from_file_data( option[ OptionKeys::mutate::mutant_file ]() );
		}

		TR << "mutations " << mutations[1] << std::endl;

		// create scorefunction or get full-atom one
		//   ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "talaris2014.wts" );
		ScoreFunctionOP sfxn = get_score_function( true );

		// loop over constructs
		for ( core::Size i = 1; i <= mutations.size(); ++i ) {

			// add mutants in that line to vectors
			add_mutant_to_vectors( input_pose, mutations[i], wt_res, resid, mut_res );

		}

		// check that the wt residues exist in the pose
		if ( check_mutants_ok( input_pose, wt_res, resid ) == false ) {
			utility_exit_with_message( "Residue identity in input file doesn't match the pose!" );
		}

		// loop over constructs
		for ( core::Size c = 1; c <= wt_res.size(); ++c ) {

			Pose pose = input_pose;

			// get base output name
			std::string const protein = option[ OptionKeys::in::file::s ](1);
			const std::string tmp1( file_basename( protein ) );
			std::string output;
			if ( option[ OptionKeys::out::path::pdb ].user() ) {
				output += to_string( option[ OptionKeys::out::path::pdb ]() ) + "/";
			}
			output += trim( tmp1, ".pdb");

			// loop over mutants within that construct
			for ( core::Size m = 1; m <= wt_res[c].size(); ++m ) {

				TR << "Looking at mutation " << std::endl;
				TR << "wt " << wt_res[c][m] << ", seqid " << resid[c][m] << ", mut " << mut_res[c][m] << std::endl;

				// make mutation
				using namespace core::chemical;
				std::string aa3 = name_from_aa( aa_from_oneletter_code( mut_res[c][m] ) );
				MutateResidueOP mutate( new MutateResidue( resid[c][m], aa3 ) );
				mutate->apply( pose );

				// add mutation to output filename
				output += "_" + to_string( wt_res[c][m] ) + to_string( resid[c][m] ) + to_string( mut_res[c][m] );

				// repack the mutated residue
				utility::vector1< bool > repack_residues( pose.total_residue(), false );
				repack_residues[ resid[c][m] ] = true;
				PackerTaskOP repack = TaskFactory::create_packer_task( pose );
				//   repack->pack_residue( static_cast< int >( seqid ) );
				repack->restrict_to_repacking();
				repack->restrict_to_residues( repack_residues );
				core::pack::pack_rotamers( pose, *sfxn, repack );

			}// mutations

			// add extension to output filename
			output += ".pdb";

			// dump pdb
			if ( utility::file::file_exists( output ) ) {
				utility_exit_with_message( "File " + output + " exists. Please delete it before rerunning this app." );
			}
			pose.dump_pdb( output );

			TR << "Wrote " << output << std::endl;
		}// constructs
	}
catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}
}
