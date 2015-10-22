// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/gzip_util.cc
/// @brief  gzip utility functions
/// @author David Kim (dekim@u.washington.edu)

// Unit header
#include <utility/file/gzip_util.hh>

// Project headers
// This has to come before boinc_util.hh or we get this error on VC++
// '_read' : is not a member of 'std::basic_istream<_Elem,_Traits>'
#include <utility/io/izstream.hh>

// Boinc headers
// This has to come before some/all other headers or we get this error on Mac:
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:34: error: storage class specifiers invalid in parameter declarations
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:34: error: storage class specified for parameter 'parameter'
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:55: error: storage class specifiers invalid in parameter declarations
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:55: error: storage class specified for parameter 'parameter'
#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#endif // BOINC

#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>


namespace utility {
namespace file {


/// @brief gzip: file compression
long
gzip(
	std::string const & uncompressedfile,
	bool overwrite
)
{
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::ifstream;
	using std::ios_base;
	using std::string;
	using utility::io::ozstream;

	string compressed = uncompressedfile + ".gz";
	string resolved_compressed = compressed;

#ifdef BOINC
	// Files that are not temporary need to have resolved names.
	boinc::resolve_filename( resolved_compressed );
#endif // BOINC

	if ( file_extension( uncompressedfile ) == "gz" ) {
		cout << "WARNING! attempt to gzip file " << uncompressedfile << " failed: already has .gz suffix -- unchanged." << endl;
		cerr << "WARNING! attempt to gzip file " << uncompressedfile << " failed: already has .gz suffix -- unchanged." << endl;
		return 0;
	}

	// Check if it exists alread
	if ( !overwrite && file_exists( resolved_compressed ) ) {
		cout << "WARNING! attempt to gzip file " << uncompressedfile << " failed: file " << resolved_compressed << " already exists." << endl;
		cerr << "WARNING! attempt to gzip file " << uncompressedfile << " failed: file " << resolved_compressed << " already exists." << endl;
		return 0;
	}

	// Check if uncompressed file exists
	if ( !file_exists( uncompressedfile ) ) {
		cout << "WARNING! attempt to gzip file " << uncompressedfile << " failed: file does not exist." << endl;
		cerr << "WARNING! attempt to gzip file " << uncompressedfile << " failed: file does not exist." << endl;
		return 0;
	}

	// Get uncompressed ifstream
	ifstream f_stream;
	trytry_ifstream_open( f_stream, uncompressedfile, ios_base::in|ios_base::binary );
	if ( !f_stream ) {
		cout << "WARNING! attempt to gzip file " << uncompressedfile << " failed: cannot read file." << endl;
		cerr << "WARNING! attempt to gzip file " << uncompressedfile << " failed: cannot read file." << endl;
		f_stream.close();
		return 0;
	}

	// Write compressed file
	// Use compressed and not resolved_compressed since ozstream
	// resolves filenames for BOINC
	ozstream zf_stream( compressed );
	if ( !zf_stream || zf_stream.uncompressed() ) {
		cout << "WARNING! attempt to create gzipped file " << resolved_compressed << " failed." << endl;
		cerr << "WARNING! attempt to create gzipped file " << resolved_compressed << " failed." << endl;
		f_stream.close();
		zf_stream.close();
		return 0;
	}

	// Write to ozstream
	zf_stream << f_stream.rdbuf();
	if ( zf_stream().fail() ) {
		cout << "WARNING! cannot write gzipped stream to file " << resolved_compressed << endl;
		cerr << "WARNING! cannot write gzipped stream to file " << resolved_compressed << endl;
		f_stream.close();
		zf_stream.close();
		return 0;
	}
	// Flush all buffers (must do for gzip)
	zf_stream.zflush();

	// Check for success and then delete uncompressed file
	long in_size = zf_stream.get_in_size();
	long out_size = zf_stream.get_out_size();
	long crc = zf_stream.get_crc();
	long uncompressedfile_size = file_size( uncompressedfile );
	f_stream.close();
	zf_stream.close();
	if ( in_size && out_size && crc && in_size == uncompressedfile_size ) {
		// gzipped file created so remove uncompressed file
		if ( file_delete( uncompressedfile ) == -1 ) {
			cout << "WARNING! error deleting file " << uncompressedfile << endl;
			cerr << "WARNING! error deleting file " << uncompressedfile << endl;
		}
		return out_size;
	}
	cout << "WARNING! gzip failed to create file " << resolved_compressed << endl;
	return 0;
}


/// @brief gunzip: file decompression
long
gunzip(
	std::string const & compressedfile,
	bool overwrite
)
{
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::ios_base;
	using std::ofstream;
	using std::string;
	using utility::io::izstream;

	// Check if compressedfile ends with .gz
	if ( file_extension( compressedfile ) != "gz" ) {
		cout << "WARNING! attempt to gunzip file " << compressedfile << " failed: unknown suffix -- ignored." << endl;
		cerr << "WARNING! attempt to gunzip file " << compressedfile << " failed: unknown suffix -- ignored." << endl;
		return 0;
	}

	// Get uncompressed file name
	string uncompressed( compressedfile );
	uncompressed.replace( uncompressed.rfind(".gz", uncompressed.length()), 3, "" );
	string resolved_uncompressed( uncompressed );

#ifdef BOINC
	// Files that are not temporary need to have resolved names.
	boinc::resolve_filename( resolved_uncompressed );
#endif // BOINC

	// Check if uncompressed file exists alread
	if ( !overwrite && file_exists( resolved_uncompressed ) ) {
		cout << "WARNING! attempt to gunzip file " << compressedfile << " failed: file " << resolved_uncompressed << " already exists." << std::endl;
		cerr << "WARNING! attempt to gunzip file " << compressedfile << " failed: file " << resolved_uncompressed << " already exists." << std::endl;
		return 0;
	}

	// Check if compressed file exists
	if ( !file_exists( compressedfile ) ) {
		cout << "WARNING! attempt to gzip file " << compressedfile << " failed: file does not exist." << endl;
		cerr << "WARNING! attempt to gzip file " << compressedfile << " failed: file does not exist." << endl;
		return 0;
	}

	// Get compressed zipstream
	izstream zf_stream( compressedfile );
	if ( !zf_stream ) {
		cout << "WARNING! attempt to gunzip file " << compressedfile << " failed: cannot read file." << endl;
		cerr << "WARNING! attempt to gunzip file " << compressedfile << " failed: cannot read file." << endl;
		zf_stream.close();
		return 0;
	}
	if ( !zf_stream.is_gzip() ) {
		cout << "WARNING! attempt to gunzip file " << compressedfile << " failed: unknown file type." << endl;
		cerr << "WARNING! attempt to gunzip file " << compressedfile << " failed: unknown file type." << endl;
		zf_stream.close();
		return 0;
	}
	// Write uncompressed file
	// Use uncompressed and not resolved_uncompressed since trytry_ofstream_open
	// Resolves filenames for BOINC
	ofstream o_stream;
	trytry_ofstream_open( o_stream, uncompressed, ios_base::out|ios_base::binary );
	if ( !o_stream ) {
		cout << "WARNING! attempt to create unzipped file " << uncompressed << " failed." << endl;
		cerr << "WARNING! attempt to create unzipped file " << uncompressed << " failed." << endl;
		o_stream.close();
		return 0;
	}
	o_stream << zf_stream.rdbuf(); // write to ozstream
	if ( o_stream.fail() ) {
		cout << "WARNING! cannot write gunzipped stream to file " << resolved_uncompressed << endl;
		cerr << "WARNING! cannot write gunzipped stream to file " << resolved_uncompressed << endl;
		o_stream.close();
		zf_stream.close();
		return 0;
	}

	long in_size = zf_stream.get_in_size();
	long out_size = zf_stream.get_out_size();
	long crc = zf_stream.get_crc();
	o_stream.close();
	zf_stream.close();
	if ( in_size && out_size && crc ) {
		// gunzipped file created so remove compressed file
		if ( file_delete( compressedfile ) == -1 ) {
			cout << "WARNING! error deleting file " << compressedfile << endl;
			cerr << "WARNING! error deleting file " << compressedfile << endl;
		}
		return out_size;
	}
	cout << "WARNING! gunzip failed to create file " << resolved_uncompressed << endl;
	cerr << "WARNING! gunzip failed to create file " << resolved_uncompressed << endl;
	return 0;
}


} // namespace file
} // namespace utility
