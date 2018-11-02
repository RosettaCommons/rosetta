// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/relax/RelaxScriptManager
/// @brief A singleton class for managing relax scripts, to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_relax_RelaxScriptManager_hh
#define INCLUDED_protocols_relax_RelaxScriptManager_hh

// Unit headers
#include <protocols/relax/RelaxScriptManager.fwd.hh>

// Protocols headers
#include <protocols/relax/FastRelax.hh>

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ header
#include <map>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

namespace protocols {
namespace relax {

/// @brief A simple wrapper class to store a vector of file contents.
/// @details Used because owning pointers to vectors behave in a wonky way.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class RelaxScriptFileContents : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor is explicitly deleted.
	RelaxScriptFileContents() = delete;

	/// @brief File contents constructor.
	RelaxScriptFileContents( utility::vector1< std::string > const & file_lines_in );

	/// @brief Destructor.
	~RelaxScriptFileContents();

	/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
	RelaxScriptFileContentsOP clone() const;

	inline utility::vector1< std::string > const & get_file_lines() const { return file_lines_; }

private:

	/// @brief Lines of the relax script file.
	utility::vector1< std::string > file_lines_;

};

/// @brief A singleton class for managing relax scripts, to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class RelaxScriptManager : public utility::SingletonBase< RelaxScriptManager > {
	friend class utility::SingletonBase< RelaxScriptManager >;

public: // Public methods //////////////////////////////////////////////////

	/// @brief Get a relax script.  Load it from disk if it has not already been loaded.
	/// @details Threadsafe and lazily loaded.
	RelaxScriptFileContents const & get_relax_script( std::string const &filename ) const;

private:  // Private methods //////////////////////////////////////////////////

	/// @brief Empty constructor.
	RelaxScriptManager();

	/// @brief Explicitly deleted copy constructor.
	RelaxScriptManager(RelaxScriptManager const & ) = delete;

	/// @brief Explicitly deleted assignment operator.
	RelaxScriptManager operator=(RelaxScriptManager const & ) = delete;

	/// @brief Create an instance of a RelaxScriptFileContents object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of RelaxScriptManager.
	static RelaxScriptFileContentsOP create_relax_script_instance( std::string const & script_file );

private:  // Private data /////////////////////////////////////////////////////

	/// @brief A map of filename to file contents.
	mutable std::map < std::string, RelaxScriptFileContentsOP > filename_to_filecontents_map_;

#ifdef MULTI_THREADED
	/// @brief Mutex for accessing the filename_to_filecontents_map_ object.
	mutable utility::thread::ReadWriteMutex relax_script_mutex_;
#endif //MULTI_THREADED


};

} //protocols
} //relax

#endif //INCLUDED_protocols/relax_RelaxScriptManager_fwd_hh



