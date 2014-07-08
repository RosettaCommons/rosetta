// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file extract_pdbs.cc
/// @brief simple application for extracting PDBs from a silent-file.
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>

#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

#include <core/pose/symmetry/util.hh>

// C++ headers
#include <string>

#include <basic/options/option_macros.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <boost/algorithm/string/erase.hpp>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char* argv [] ) {
	try {

	// define relevant options
	OPT(in::path::database);
	OPT(in::file::silent);
	OPT(in::file::tags);
	OPT(in::file::silent_struct_type);
	OPT(in::file::silent_renumber);
	OPT(out::file::residue_type_set);
	OPT(out::prefix);
	OPT(in::path::database);
	OPT(in::file::rescore);
	OPT(score::weights);
	OPT(score::patch);

	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// options, random initialization
	devel::init( argc, argv );

	std::string usage("");
	usage += "\n\nusage:  extract_pdbs [options] -in::file::silent <silent_files>\n";
	usage += "\tTo see a list of other valid options, use the option -help.\n";

	if ( !option[ in::file::silent ].user() ) {
		std::cerr << usage << std::endl;
		std::exit(1);
	}

	basic::Tracer tr( "extract_pdbs" );

	core::chemical::ResidueTypeSetCAP rsd_set( NULL );
	if ( 	option[ out::file::residue_type_set ].user() ) {
		rsd_set = ChemicalManager::get_instance()->residue_type_set(
					option[ out::file::residue_type_set ]()
		);
	}
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
			);
		} else if( option[ in::file::tagfile ].user() ) {

			utility::vector1< std::string > input_tags;
			utility::io::izstream tag_file( option[ in::file::tagfile ]() );
			std::copy( std::istream_iterator< std::string >( tag_file ), std::istream_iterator< std::string >(),
								 std::back_inserter( input_tags ) );
			input = new SilentFilePoseInputStream( option[ in::file::silent ](), input_tags );

		} else {
			input = new SilentFilePoseInputStream( option[ in::file::silent ]() );
		}
	}

	core::scoring::ScoreFunctionOP scorefxn;

	if ( option[ in::file::rescore ]() ){
		scorefxn = core::scoring::get_score_function();
		tr.Debug << "scoring using ScoreFunction with weights: " << std::endl;
		scorefxn->show( tr.Debug );
		tr.flush();
	}

	core::pose::Pose pose;
	std::string out_prefix = option[ out::prefix ]();
	while ( input->has_another_pose() ) {
		//fpd  'fill_pose' assumes an asymmetric pose is input
    if ( core::pose::symmetry::is_symmetric( pose ) ) {
      core::pose::symmetry::make_asymmetric_pose( pose );
    }
		if ( rsd_set ) {
			input->fill_pose( pose, *rsd_set );
		} else {
			input->fill_pose( pose );
		}

		// make sure the parent/template name gets included in the REMARK line of the PDB
		using std::map;
		using std::string;
		map< string, string > score_line_strings(	core::pose::get_all_score_line_strings( pose ) );
		for ( map< string, string >::const_iterator it = score_line_strings.begin(),
									end = score_line_strings.end();
									it != end; ++it )
		{
					if ( it->first != "aln_id" ) continue;
					core::pose::add_comment( pose, "parents", it->second );
		}
		std::string tag( tag_from_pose( pose ) );
		std::string fn( out_prefix + tag + ".pdb" );

		tr << "extracting Pose with tag " << tag << " into PDB file " << fn
			<< std::endl;

		if ( option[ in::file::rescore ]() ) {
			tr.Debug << "rescoring Pose with tag " << tag << std::endl;
			scorefxn->show( tr.Debug );
			(*scorefxn)(pose);

			scorefxn->show( tr.Debug, pose );
			pose.dump_scored_pdb( fn, *scorefxn );
		} else {
			pose.dump_pdb( fn );
		}

		//fpd

		tr.flush();
	}

	exit( 0 );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // main
