// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/ozstream.cc
/// @brief  Output file stream wrapper for uncompressed and compressed files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author David Kim (dekim@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


// Unit header
#include <utility/io/ozstream.hh>

// Project headers
#include <utility/file/file_sys_util.hh>
#if defined( USE_FILE_PROVIDER )
#include <utility/inline_file_provider.hh>
#endif
#include <utility/exit.hh>

// C++ headers
#include <cstdlib>


namespace utility {
namespace io {

/// @detail special "atomic" method to open a stream and print a header if it is new. This is necessary to avoid process competition when
/// checking for existance to decided whether or not to write a header.
void
ozstream::open_append_if_existed( std::string const& filename_a, std::stringstream& preprinted_header ) {
	// Close the file if open and reset the state
	close();
	filename_ = filename_a;
#ifdef USEMPI
		//		std::cout << "MPI_Reroute " << (bMPI_reroute_stream_ ? " active " : " not-active ") << std::endl;
	if ( bMPI_reroute_stream_ ) { // this is switched via call to static function enable_MPI_reroute()
		mpi_stream_p_ = new mpi_stream::mpi_ostream( filename_a, mpi_FileBuf_master_rank_, preprinted_header, true );
		if ( ( !mpi_stream_p_ ) || ( !( *mpi_stream_p_ ) ) ) {
			compression_ = NONE;
		} else 	compression_ = UNCOMPRESSED;
		return;
	}
#endif

	if ( !utility::file::file_exists( filename_a ) ) {
		open( filename_a );
		if ( good() ) ( *this) << preprinted_header.str();
	} else {
		open_append( filename_a );
	}
}

/// @brief Open a file
void
ozstream::open(
	std::string const & filename_a,
	std::ios_base::openmode open_mode // Ignored for gzip files
)
{
	using std::ios;
	using std::ios_base;
	using utility::file::file_extension;
	using utility::file::trytry_ofstream_open;
	using zlib_stream::zip_ostream;


	// #if defined( USE_FILE_PROVIDER )
	//  utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
	//
	//  if(!provider->get_ostream( filename_a , &file_provider_stream )){
	//    std::cerr << "Cannot find inline file: " << filename_a << std::endl;
	//    file_provider_stream = &bad_stream;
	//    file_provider_stream->setstate( ios_base::failbit | ios_base::badbit );
	//  }
	//
	//    if (file_provider_stream->good() ){
	//      return;
	//    }
	// #endif

	// Close the file if open and reset the state
	close();

	// Open the ofstream and set the compression state and file name
	filename_ = filename_a;

#ifdef USEMPI
		//		std::cout << "MPI_Reroute " << (bMPI_reroute_stream_ ? " active " : " not-active ") << std::endl;
		if ( bMPI_reroute_stream_ ) {
			std::stringstream no_header;
			mpi_stream_p_ = new mpi_stream::mpi_ostream( filename_, mpi_FileBuf_master_rank_, no_header, open_mode & ios::app );
			if ( ( !mpi_stream_p_ ) || ( !( *mpi_stream_p_ ) ) ) {
				compression_ = NONE;
			} else {
				compression_ = UNCOMPRESSED;
			}
			return;
		}
#endif
	if ( ( open_mode & ios::ate ) || ( open_mode & ios::app )
			|| ( ( open_mode & ios::in ) && ( open_mode & ios::out ) ) ) {

		// prepare new character buffer -- must do this before file is opened
		allocate_assign_char_buffer();

		// Unsupported for gzip files: Use ofstream
		trytry_ofstream_open( of_stream_, filename_a, open_mode );

		compression_ = UNCOMPRESSED;

	} else if ( file_extension( filename_a ) == "gz" ) { // gzip file

		trytry_ofstream_open( of_stream_, filename_a, ios_base::out|ios_base::binary );
		if ( of_stream_ ) { // Open succeeded
			compression_ = GZIP;
		} else { // Leave stream state so that failure can be detected
			compression_ = NONE;
		}

	} else { // Uncompressed file

		// prepare new character buffer -- must do this before file is opened
		allocate_assign_char_buffer();

		trytry_ofstream_open( of_stream_, filename_a, ios_base::out );
		if ( of_stream_ ) { // Open succeeded
			compression_ = UNCOMPRESSED;
		} else { // Leave stream state so that failure can be detected
			compression_ = NONE;
		}

	}

	// Attach zip_ostream to ofstream if gzip file
	if ( compression_ == GZIP ) {
		// zip_stream_p_ deleted by close() above so don't have to here
		zip_stream_p_ = new zip_ostream( of_stream_, true, static_cast< size_t >( Z_DEFAULT_COMPRESSION ), zlib_stream::DefaultStrategy, 15, 8, buffer_size_ );
		if ( ( !zip_stream_p_ ) || ( !( *zip_stream_p_ ) ) ||
				( !zip_stream_p_->is_gzip() ) ) { // zip_stream not in good state
			if ( zip_stream_p_ ) delete zip_stream_p_;
			zip_stream_p_ = nullptr;
			of_stream_.close();
			// Set failbit so failure can be detected
			of_stream_.setstate( ios_base::failbit );
		}
	}
}


/// @brief Open a text file or gzip'd file for appending
void
ozstream::open_append( std::string const & filename_a )
{
	using std::cout;
	using std::endl;
	using std::exit;
	using std::ios_base;
	using utility::file::file_extension;
	using utility::file::trytry_ofstream_open;
	using zlib_stream::zip_ostream;

	// Close the file if open and reset the state
	close();

	// Open the ofstream and set the compression state and file name
	filename_ = filename_a;

#ifdef USEMPI
		//		std::cout << "MPI_Reroute " << (bMPI_reroute_stream_ ? " active " : " not-active ") << std::endl;
		if ( bMPI_reroute_stream_ ) {
			std::stringstream no_header;
			mpi_stream_p_ = new mpi_stream::mpi_ostream( filename_, mpi_FileBuf_master_rank_, no_header, true );
			if ( ( !mpi_stream_p_ ) || ( !( *mpi_stream_p_ ) ) ) {
				compression_ = NONE;
			} else compression_ = UNCOMPRESSED;
			return;
		}
#endif

	if ( file_extension( filename_a ) == "gz" ) { // gzip file

		trytry_ofstream_open( of_stream_, filename_a, ios_base::out|ios_base::binary|ios_base::app );
		if ( of_stream_ ) { // Open succeeded
			compression_ = GZIP;
		} else { // Leave stream state so that failure can be detected
			compression_ = NONE;
		}

	} else { // Uncompressed file

		// prepare new character buffer -- must do this before file is opened
		allocate_assign_char_buffer();

		trytry_ofstream_open( of_stream_, filename_a, ios_base::out|ios_base::app );
		if ( of_stream_ ) { // Open succeeded
			compression_ = UNCOMPRESSED;
		} else { // Leave stream state so that failure can be detected
			compression_ = NONE;
		}

	}

	// Attach zip_ostream to ofstream if gzip file
	if ( compression_ == GZIP ) {
		// zip_stream_p_ deleted by close() above so don't have to here
		zip_stream_p_ = new zip_ostream( of_stream_, true, static_cast< size_t >( Z_DEFAULT_COMPRESSION ), zlib_stream::DefaultStrategy, 15, 8, buffer_size_ );
		if ( ( !zip_stream_p_ ) || ( !( *zip_stream_p_ ) ) ||
				( !zip_stream_p_->is_gzip() ) ) { // zip_stream not in good state
			delete zip_stream_p_; zip_stream_p_ = nullptr;
			of_stream_.close();
			// Set failbit so failure can be detected
			of_stream_.setstate( ios_base::failbit );
		}
	}
}

#ifndef USEMPI
void ozstream::enable_MPI_reroute( int, int ) {
	utility_exit_with_message( "enable_MPI_reroute called in non-mpi version of mini ");
}
#endif

#ifdef USEMPI
void ozstream::enable_MPI_reroute( int min_rank, int master_rank ) {
	int my_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);/* get current process id */
	if ( my_rank >= min_rank ) bMPI_reroute_stream_ = true;
	mpi_FileBuf_master_rank_ = master_rank;
}
#endif

bool ozstream::bMPI_reroute_stream_( false );
int  ozstream::mpi_FileBuf_master_rank_( 0 );

} // namespace io
} // namespace utility
