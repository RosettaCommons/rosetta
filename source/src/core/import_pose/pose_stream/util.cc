// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pose_stream/util.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_core_pose_stream_util_HH
#define INCLUDED_core_pose_stream_util_HH


#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/ExtendedPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/LazySilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <string>

#include <utility/io/izstream.hh>

#include <utility/vector1.hh>
#include <boost/bind.hpp>


namespace core {
namespace import_pose {
namespace pose_stream {


utility::vector1< utility::file::FileName > filenames_from_list_file(
	utility::vector1< utility::file::FileName > const & list_fns
) {

	using std::string;
	using utility::vector1;
	using utility::file::FileName;

	vector1< FileName > fns_to_return;

	for ( vector1< FileName >::const_iterator it = list_fns.begin(),
			it_end = list_fns.end();
			it != it_end; ++it
			) {
		std::string filename( it->name() );
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}

		std::string line;
		while ( getline(data, line) ) {
			fns_to_return.push_back( utility::file::FileName(line) );
		}
		data.close();
	}

	return fns_to_return;
}

/// @brief Get all input streams based on command-line input.
/// @details If do_renumber_decoys is true, silent file decoys are sorted in alphabetical order of tags.
MetaPoseInputStream streams_from_cmd_line( bool const do_renumber_decoys ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	MetaPoseInputStream input;
	if ( option[ in::file::s ].user() ) {
		PDBPoseInputStreamOP pdb_input( new PDBPoseInputStream( option[ in::file::s ]() ) );
		input.add_pose_input_stream( pdb_input );
	}
	if ( option[ in::file::l ].user() ) {
		PDBPoseInputStreamOP pdb_list_input( new PDBPoseInputStream( filenames_from_list_file( option[ in::file::l ]() ) ) );
		input.add_pose_input_stream( pdb_list_input );
	}

	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::lazy_silent ]() ) {
			LazySilentFilePoseInputStreamOP silent_input( new LazySilentFilePoseInputStream( option[ in::file::silent ]() ) );
			input.add_pose_input_stream(silent_input);

			if ( option[ in::file::silent_list ].user() ) {
				LazySilentFilePoseInputStreamOP silent_list_input( new LazySilentFilePoseInputStream( filenames_from_list_file(
					option[ in::file::silent_list ]()
					) ) );
				input.add_pose_input_stream(silent_list_input);
			}
		} else {
			SilentFilePoseInputStreamOP silent_input;
			if ( option[ in::file::tags ].user() ) {
				silent_input = SilentFilePoseInputStreamOP( new SilentFilePoseInputStream(
					option[ in::file::silent ](),
					option[ in::file::tags ](),
					option[ in::file::silent_energy_cut ]()
					) );
			} else {
				silent_input = SilentFilePoseInputStreamOP( new SilentFilePoseInputStream(
					option[ in::file::silent ](), option[ in::file::silent_energy_cut ]()
					) );
			}
			silent_input->renumber_decoys( do_renumber_decoys && option[ in::file::silent_renumber ]() );
			input.add_pose_input_stream( silent_input );
		}
	}

	if ( option[ in::file::extended_pose ].user() ) {
		Size const ntimes = option[ in::file::extended_pose ]();

		using utility::vector1;
		using utility::file::FileName;

		vector1< FileName > fasta_files( option[ in::file::fasta ]() );
		for ( vector1< FileName >::const_iterator it = fasta_files.begin(),
				end = fasta_files.end(); it != end; ++it
				) {

			std::string sequence
				= core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence();

			PoseInputStreamOP extended_protein_input( new ExtendedPoseInputStream( sequence, ntimes ) );
			input.add_pose_input_stream( extended_protein_input );
		}
	}

	return input;
} // streams_from_cmd_line

/// @brief Get all input streams based on command-line input, sorting silent file decoys in alphabetical order of tags.
///
MetaPoseInputStream streams_from_cmd_line () {
	return streams_from_cmd_line(true);
}


} // pose_stream
} // import_pose
} // core

#endif
