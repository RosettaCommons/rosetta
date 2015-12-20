// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/file_sys_util.cc
/// @brief  Platform independent operations on files (except I/O)
/// @author David Kim (dekim@u.washington.edu)
/// @author Ion Yannopoulos (ion@rosettacommons.org)
/// @todo   Break out platform-specific code.
/// @todo   trytry_* functions create mutual dependency between
///         utility/file and utility/io.  Resolve this.

// Boinc headers
// This has to come before some/all other headers or we get this error on Mac:
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:34: error: storage class specifiers invalid in parameter declarations
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:34: error: storage class specified for parameter 'parameter'
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:55: error: storage class specifiers invalid in parameter declarations
//   /System/Library/Frameworks/CoreFoundation.framework/Headers/CFMessagePort.h:55: error: storage class specified for parameter 'parameter'
#ifdef USEMPI
#include <mpi.h>
#endif

#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#endif // BOINC

// Unit headers
#include <utility/file/file_sys_util.hh>
#include <utility/file/PathName.hh>

// Project headers
#include <utility/basic_sys_util.hh> //  for rand_sleep?
#include <utility/exit.hh>

// C++ headers
#include <iostream>

// Platforms headers
#include <sys/stat.h>
#ifndef _WIN32
#include <dirent.h>
#endif

// only needed for linux?
#include <errno.h>

#ifdef CXX11
#include <random>
#endif

// Platform headers - Win32
#ifndef _WIN32
// #include <unistd.h>
#endif // _WIN32

// C POSIX library header for directory manipulation
#if (defined WIN32) //&& (!defined WIN_PYROSETTA)
#include <windows.h>
#include <direct.h>
#include <io.h>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <io.h>
#include <sys/types.h>
#endif

#ifdef WIN_PYROSETTA
#include <direct.h>
#define mkdir _mkdir
#endif

// Windows compatibility shim (*nixes should have it already defined.)
#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif


namespace utility {
namespace file {

/// @brief Does File Exist?
bool
file_exists( std::string const & path )
{
	// NOTE: this is not entirely reliable, stat may fail also when a file does exist
#if defined USE_FILE_PROVIDER
	utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
	// if file exists in file_provider then this will be the one used. Thus return true - file exists.
	if ( provider->file_exists( path ) ) return true;
	//otherwise try and locate the file on the file system
#endif

#ifdef _WIN32
	/*bool res = false;
	std::fstream f;
	f.open( path.c_str(), std::ios::in );
	if( f.is_open() ) res = true;
	f.close();
	return res; */

	if ( access(path.c_str(), 0) ) return false;
	else return true;

#else
	struct stat buf;
	return !stat( path.c_str(), &buf ); // stat() returns zero on success
#endif
}


/// @brief Is the given path a directory?
bool
is_directory( std::string const & path ) {
	// Does this need to be specialized for Windows builds, etc.?
	// As far as I can tell, not necessarily.
	struct stat buf;
	if ( stat(path.c_str(), &buf) ) { // stat() returns non-zero on failure
		return false; // A non-existant file is not a directory.
	}
	return S_ISDIR(buf.st_mode);
}

/// @brief Delete File
int
file_delete( std::string const & path )
{
	if ( !file_exists( path ) ) {
		return 0;
	}

	int retval;
#ifdef _WIN32
	for ( int i = 0; i < 5; ++i ) {
		retval = remove( path.c_str() );
		if ( !retval ) break;
		rand_sleep(); // avoid lockstep
	}
#else
	retval = remove( path.c_str() );
#endif // _WIN32
	return retval;
}


/// @brief Extension of a File Name
std::string
file_extension( std::string const & filename )
{
	using std::string;

	string::size_type const iver = filename.find_last_of( '~' );
	if ( iver == string::npos ) { // No version
		string::size_type const iext = filename.find_last_of( '.' );
		if ( iext == string::npos ) { // No extension
			return string();
		} else { // Return the extension (without the '.')
			return filename.substr( iext + 1 );
		}
	} else { // Find extension before version
		string::size_type const iext = filename.find_last_of( '.', iver );
		if ( iext == string::npos ) { // No extension
			return string();
		} else { // Return the extension (without the '.')
			return filename.substr( iext + 1, iver - iext - 1 );
		}
	}
}


/// @brief Basename of a File Name
std::string
file_basename( std::string const & filename )
{
	using std::string;

	string::size_type const dot_location = filename.find_last_of( '.' );
	if ( dot_location == string::npos ) { // No extension
		return filename;
	} else { // Return the filename without the extension
		return filename.substr( 0, dot_location );
	}
}


/// @brief Platform independent way of getting file size
long
file_size( std::string const & filename )
{
	using std::cerr;
	using std::endl;
	using std::ifstream;
	using std::ios_base;

	// Skip if file does not exist
	if ( !file_exists( filename ) ) return -1;
	long l, m;
	ifstream file_;
	if ( !trytry_ifstream_open( file_, filename, ios_base::in|ios_base::binary ) ) {
		cerr << "WARNING! cannot get file size for " << filename << ": could not open file." << endl;
		return -1;
	}
	l = file_.tellg();
	file_.seekg( 0, ios_base::end );
	m = file_.tellg();
	file_.close();
	return m - l;
}


/// @brief Create a blank file if it doesn't already exist
bool
create_blank_file( std::string const & blank_file )
{
	//utility::io::ozstream blank_out_stream;
	//blank_out_stream.open( blank_file );
	std::ofstream blank_out_stream;
	trytry_ofstream_open( blank_out_stream, blank_file, std::ios::out );
	if ( !blank_out_stream ) {
		std::cout << "Open failed for file: " << blank_file << std::endl;
		std::cerr << "Open failed for file: " << blank_file << std::endl;
		utility_exit();
	}
	blank_out_stream << std::endl;
	blank_out_stream.close();
	blank_out_stream.clear();
	return true;
}

/// @brief Find an unused random tempfile name with a given prefix (which may include a directory)
std::string
create_temp_filename( std::string const & dir, std::string const & prefix ) {
	// Ugly C stuff...
#ifdef WIN32
	char * tmpname_output = _tempnam( dir.c_str(), prefix.c_str() );
	std::string tempfilename = tmpname_output;
	free( tmpname_output );
#elif __CYGWIN__
	char * tmpname_output = tempnam( dir.c_str(), prefix.c_str() );
	std::string tempfilename = tmpname_output;
	free( tmpname_output );
#else
	std::string dirname = dir;
	if ( *dirname.rbegin() != '/' ) {
		dirname += "/";
	}
	bool exists = false;
	std::string const charset("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	std::string str = "";
	std::string tempfilename;
#ifdef CXX11
	std::random_device rd; // Gets unseeded random numbers from OS.
	std::mt19937 mt(rd()); // Start the mt with the random device.
	std::uniform_int_distribution< int > rand_int(0,charset.size()-1); // Inclusive range
#else
	srand(time(NULL));
#endif
	do {
		// We're using the C++ random number generator because we don't want
		// to pull extra draws on the Rosetta one - besides,
		// in utility we don't have access to it
#ifdef CXX11
		str += charset[rand_int(mt)];
#else
		str += charset[rand()%charset.size()];
#endif
		tempfilename  = dirname + str + prefix;
		std::ifstream file(tempfilename.c_str());
		if ( file ) {
			file.close();
			exists = true;
		} else {
			exists = false;
		}
	} while (exists);
#endif
	return tempfilename;
}

/// @details
/// Code is based on a discussion on this web site:
///
/// http://www.brainbell.com/tutorials/C++/Creating_A_Directory.htm
///
/// On non-Windows systems, the permissions are solely determined by the umask.
/// If directory creation is successful, the function returns true.
///
/// *** Currently untested on Windows! If you find it works, remove this warning. Thanks. ***
bool
create_directory(
	std::string const & dir_path
)
{
#if !defined( __native_client__ )
#if (defined WIN32)
	// Windows code
	return !mkdir(dir_path.c_str());
#else
	// non-Windows code
	return !mkdir(dir_path.c_str(), 0777);
#endif // _WIN32
#endif
	return false;
}


/// @details
/// Based on the above create_directory() function.
/// If directory creation is successful or the directory already exists, the function returns true.
///
/// *** Currently untested on Windows! If you find it works, remove this warning. Thanks. ***
bool
create_directory_recursive(
	std::string const & dir_path
)
{
	// This code works on Linux and Mac, and it *appears* that Windows uses the same error codes:
	// http://msdn.microsoft.com/en-us/library/2fkk4dzw(VS.80).aspx
	bool success = create_directory(dir_path);
	if ( success ) return true;
	else if ( errno == EEXIST ) return true; // directory already exists
	else if ( errno == ENOENT ) {
		std::string const parent = PathName(dir_path).parent().name();
		if ( parent == dir_path ) return false;
		else return create_directory_recursive(parent) && create_directory(dir_path);
	} else return false;
}


/// @brief Try to open file a few times just in case it is locked (from BOINC LIB)
bool
trytry_ifstream_open(
	std::ifstream & ifstream_,
	std::string const & name,
	std::ios_base::openmode open_mode
)
{
	using std::ios;
	using std::ios_base;
	using std::string;

	ifstream_.close();
	ifstream_.clear();

	string resolved_name( name );

#ifdef BOINC
	// Files that are not temporary need to have resolved names.
	// Resolve them here since all file input should use this function.
	// Be sure to resolve file names used in other functions
	boinc::resolve_filename( resolved_name );
#endif // BOINC

	// If opening for read, and file isn't there,
	// leave now (avoid 5-second delay!!)
	if ( open_mode & ios::in ) {
		if ( !file_exists( resolved_name ) ) {
			ifstream_.setstate( ios_base::failbit ); // set ios state to failbit
			return false;
		}
		// On some platforms, "opening" a directory as a file for reading succeeds and gives an empty "file".
		// Catch this and error out instead.
		if ( is_directory( resolved_name ) ) {
			ifstream_.setstate( ios_base::failbit ); // set ios state to failbit
			return false;
		}
	}

	ifstream_.open( resolved_name.c_str(), open_mode );
	if ( ifstream_ && ifstream_.is_open() ) return true;

#ifdef _WIN32
	// On Windows: if open fails, try again for 5 seconds
	// (since the file might be open by FastFind, Diskeeper etc.)
	if ( !ifstream_ ) {
		for ( int i = 0; i < 5; ++i ) {
			rand_sleep();
			ifstream_.open( resolved_name.c_str(), open_mode );
			if ( ifstream_ && ifstream_.is_open() ) return true;
		}
	}
#else
	// Unix - if call was interrupted, retry a few times
	if ( !ifstream_ ) {
		for ( int i = 0; i < 5; ++i ) {
			if ( errno != EINTR ) break;
			ifstream_.open( resolved_name.c_str(), open_mode );
			if ( ifstream_ && ifstream_.is_open() ) return true;
		}
	}
#endif // _WIN32
	ifstream_.setstate( ios_base::failbit ); // set ios state to failbit
	return false;
}


/// @brief Try to open file a few times just in case it is locked (from BOINC LIB)
bool
trytry_ofstream_open(
	std::ofstream & ofstream_,
	std::string const & name,
	std::ios_base::openmode open_mode
)
{
	using std::ios;
	using std::ios_base;
	using std::string;

	ofstream_.close();
	ofstream_.clear();

	string resolved_name( name );

#ifdef BOINC
	// Files that are not temporary need to have resolved names.
	// Resolve them here since all file output should use this function.
	// Be sure to resolve file names used in other functions
	boinc::resolve_filename( resolved_name );
#endif // BOINC

	// If opening for read, and file isn't there,
	// leave now (avoid 5-second delay!!)
	if ( open_mode & ios::in ) {
		if ( !file_exists( resolved_name ) ) {
			ofstream_.setstate( ios_base::failbit ); // set ios state to failbit
			return false;
		}
	}

	ofstream_.open( resolved_name.c_str(), open_mode );
	if ( ofstream_ && ofstream_.is_open() ) return true;

#ifdef _WIN32
	// On Windows: if open fails, try again for 5 seconds
	// (since the file might be open by FastFind, Diskeeper etc.)
	if ( !ofstream_ ) {
		for ( int i = 0; i < 5; ++i ) {
			rand_sleep();
			ofstream_.open( resolved_name.c_str(), open_mode );
			if ( ofstream_ && ofstream_.is_open() ) return true;
		}
	}
#else
	// Unix - if call was interrupted, retry a few times
	if ( !ofstream_ ) {
		for ( int i = 0; i < 5; ++i ) {
			if ( errno != EINTR ) break;
			ofstream_.open( resolved_name.c_str(), open_mode );
			if ( ofstream_ && ofstream_.is_open() ) return true;
		}
	}
#endif // _WIN32
	ofstream_.setstate( ios_base::failbit ); // set ios state to failbit
	return false;
}

int list_dir (std::string dir, utility::vector1<std::string> & files)
{
	//#ifndef WIN_PYROSETTA
#if (defined WIN32) //&& (!defined PYROSETTA)
	std::string file_name = dir + "\\*";

	WIN32_FIND_DATA find_file_data;
	HANDLE h_find;

	h_find = FindFirstFile(file_name.c_str(), &find_file_data);
	if ( h_find == INVALID_HANDLE_VALUE ) {
		return -1;
	}

	do {
		files.push_back(find_file_data.cFileName);
	} while(FindNextFile(h_find, &find_file_data) != 0);
	FindClose(h_find);
#else
			DIR *dp;
			struct dirent *dirp;
			if((dp  = opendir(dir.c_str())) == NULL) {
			 //  cout << "Error(" << errno << ") opening " << dir << endl;
				return errno;
			}

			while ((dirp = readdir(dp)) != NULL) {
				files.push_back(std::string(dirp->d_name));
			}
			closedir(dp);
#endif
	//#endif
	return 0;
}

FileName combine_names(utility::vector1<std::string> file_name_strings){
	std::vector<FileName> file_names;
	std::vector<std::string>::const_iterator begin= file_name_strings.begin();
	for ( ; begin != file_name_strings.end(); ++begin ) {
		file_names.push_back(FileName(*begin));
	}
	return FileName(file_names);
}


} // namespace file
} // namespace utility
