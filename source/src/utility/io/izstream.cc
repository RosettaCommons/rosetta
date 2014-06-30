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

#if defined( USE_FILE_PROVIDER )
#include <utility/inline_file_provider.hh>
#endif

namespace utility {
namespace io {

	// Initialize private static data
	vector1< std::string > izstream::alternative_search_paths_;

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

		//std::cout << "FILELIST: Trying to open file: '" << filename_a << "' " << std::endl;

#if defined( USE_FILE_PROVIDER )
	utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
	if( (open_mode & std::ios_base::out ) ){
		throw( "Cannot open output file in istream inline file provider " );
	}
	
	if(!provider->get_istream( filename_a , &file_provider_stream )){
		 std::cerr << "Cannot find inline file: " << filename_a << std::endl;	
		 file_provider_stream = &bad_stream;
		 file_provider_stream->setstate( ios_base::failbit | ios_base::badbit );
	}

#ifdef __native_client__
  // in native client mode there is no recourse - no actauly file access is allowed
  // so ut this point we can only hand back the stream no matter if it is good() or not.
  return;
#endif

  if (file_provider_stream->good() ){ 
    return;
  }
#endif

		// Close the file if open and reset the state
		close();

		// Open the ifstream and set the compression state and file name
		if ( ( open_mode & ios_base::ate ) ||
		 ( open_mode & ios_base::app ) ||
		 ( ( open_mode & ios_base::in ) &&
		 ( open_mode & ios_base::out ) ) ) { // Unsupported for gzip files: Use ifstream
			open_ifstream(filename_a, open_mode);
			compression_ = UNCOMPRESSED;
		} else if ( file_extension( filename_a ) == "gz" ) { // gzip file
			open_ifstream(filename_a, ios_base::in|ios_base::binary );
			if( if_stream_ ){ // Open succeeded
				compression_ = GZIP;
			} else { // Leave stream state so that failure can be detected
				compression_ = NONE;
			}
		} else { // Try with .gz extension added
			open_ifstream( filename_a + ".gz", ios_base::in|ios_base::binary );
			if ( if_stream_ ) { // Found/opened with .gz added
				compression_  = GZIP;
			} else { // Try as an uncompressed file
				if_stream_.clear();
				open_ifstream( filename_a, open_mode );
				if ( if_stream_ ) { // Open succeeded
					compression_ = UNCOMPRESSED;
				} else { // Leave stream state so that failure can be detected
					compression_ = NONE;
				}
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


	void
	izstream::open_ifstream(
		std::string const & name,
		std::ios_base::openmode open_mode
	){
		using std::ios_base;
		using utility::file::trytry_ifstream_open;
		trytry_ifstream_open( if_stream_, name, open_mode );
		if ( if_stream_ ) { // Open succeeded
			filename_ = name;
		} else { // Try opening file in alternative search paths
			vector1<std::string>::const_iterator i(alternative_search_paths_.begin());
			vector1<std::string>::const_iterator const ie(alternative_search_paths_.end());
			for( ; i!=ie; ++i ){
				trytry_ifstream_open( if_stream_,
					(*i + platform::file::PATH_SEPARATOR) + name,
					ios_base::in|ios_base::binary );
				if(if_stream_){
					filename_ = (*i + platform::file::PATH_SEPARATOR) + name;
					break;
				}
			}
		}
	}


} // namespace io
} // namespace utility
