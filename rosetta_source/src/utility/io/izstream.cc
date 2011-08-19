// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/izstream.cc
/// @brief  Input file stream wrapper for uncompressed and compressed files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author David Kim (dekim@u.washington.edu)


// Unit header
#include <utility/io/izstream.hh>

// Project headers
#include <utility/file/file_sys_util.hh>


namespace utility {
namespace io {


	/// @brief Open a file
	void
	izstream::open(
		std::string const & filename_a,
		std::ios_base::openmode open_mode // Ignored for gzip files
	)
	{
		using std::ios_base;
		using utility::file::file_extension;
		using utility::file::trytry_ifstream_open;
		using zlib_stream::zip_istream;

		// Close the file if open and reset the state
		close();

		// Open the ifstream and set the compression state and file name
		if ( ( open_mode & ios_base::ate ) ||
		 ( open_mode & ios_base::app ) ||
		 ( ( open_mode & ios_base::in ) &&
		 ( open_mode & ios_base::out ) ) ) { // Unsupported for gzip files: Use ifstream
			trytry_ifstream_open( if_stream_, filename_a, open_mode );
			compression_ = UNCOMPRESSED;
			filename_ = filename_a;
		} else if ( file_extension( filename_a ) == "gz" ) { // gzip file
			trytry_ifstream_open( if_stream_, filename_a, ios_base::in|ios_base::binary );
			if ( if_stream_ ) { // Open succeeded
				compression_ = GZIP;
			} else { // Leave stream state so that failure can be detected
				compression_ = NONE;
			}
			filename_ = filename_a;
		} else { // Try with .gz extension added
			trytry_ifstream_open( if_stream_, filename_a + ".gz", ios_base::in|ios_base::binary );
			if ( if_stream_ ) { // Found/opened with .gz added
				compression_  = GZIP;
				filename_ = filename_a + ".gz";
			} else { // Try as an uncompressed file
				if_stream_.clear();
				trytry_ifstream_open( if_stream_, filename_a, open_mode );
				if ( if_stream_ ) { // Open succeeded
					compression_ = UNCOMPRESSED;
				} else { // Leave stream state so that failure can be detected
					compression_ = NONE;
				}
				filename_ = filename_a;
			}
		}

		// Attach zip_istream to ifstream if gzip file
		if ( compression_ == GZIP ) {
			// zip_stream_p_ deleted by close() above so don't have to here
			zip_stream_p_ = new zip_istream( if_stream_ );
			if ( ( !zip_stream_p_ ) || ( !( *zip_stream_p_ ) ) || ( !zip_stream_p_->is_gzip() ) ) {
				// zip_stream not in good state
				delete zip_stream_p_; zip_stream_p_ = 0;
				if_stream_.close();
				if_stream_.setstate( ios_base::failbit ); // set failbit so failure can be detected
			}
		}
	}


} // namespace io
} // namespace utility
