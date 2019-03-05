// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/io/GeneralFileManager.hh
/// @brief A singleton class for managing arbitrary files to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_utility_io_GeneralFileManager_hh
#define INCLUDED_utility_io_GeneralFileManager_hh

// Unit headers
#include <utility/io/GeneralFileManager.fwd.hh>

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <string>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

namespace utility {
namespace io {

class GeneralFileContents : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor is explicitly deleted.
	GeneralFileContents() = delete;

	/// @brief File contents constructor.
	GeneralFileContents( std::string const & filename );

	/// @brief Destructor.
	~GeneralFileContents();

	/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
	GeneralFileContentsOP clone() const;

	std::string const & get_file_contents() const {
		return file_contents_;
	}

private:

	std::string file_contents_;

};

class GeneralFileManager : public utility::SingletonBase< GeneralFileManager > {
	friend class utility::SingletonBase< GeneralFileManager >;

public:

	/// @brief Get a weights file.  Load it from disk if it has not already been loaded.
	/// @details Threadsafe and lazily loaded.
	std::string const & get_file_contents( std::string const & filename ) const;

private:

	/// @brief Empty constructor.
	GeneralFileManager();

	/// @brief Explicitly deleted copy constructor.
	GeneralFileManager(GeneralFileManager const & ) = delete;

	/// @brief Explicitly deleted assignment operator.
	GeneralFileManager operator=(GeneralFileManager const & ) = delete;

	/// @brief Create an instance of a GeneralFileContents object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of GeneralFileManager.
	static GeneralFileContentsOP create_instance( std::string const & script_file );

private:  // Private data /////////////////////////////////////////////////////

	/// @brief A map of filename to file contents.
	mutable std::map < std::string, GeneralFileContentsOP > filename_to_filecontents_map_;

#ifdef MULTI_THREADED
	/// @brief Mutex for accessing the filename_to_filecontents_map_ object.
	mutable utility::thread::ReadWriteMutex io_script_mutex_;
#endif //MULTI_THREADED


};

} //utility
} //io

#endif //INCLUDED_utility/io_GeneralFileManager_fwd_hh



