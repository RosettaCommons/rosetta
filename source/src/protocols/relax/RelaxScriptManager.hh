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

// Core headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.hh>

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ header
#include <map>
#include <tuple>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

namespace protocols {
namespace relax {

/// @brief Represents a single line of database/sampling/relax_scripts/index.dat
struct DatabaseRelaxScript {
	std::string title;
	bool sfxn_specialized;
	utility::vector1< std::string > sfxns_supported;

	DatabaseRelaxScript() = delete;

	DatabaseRelaxScript(
		std::string const & title_arg,
		bool sfxn_specialized_arg
	) :
		title( title_arg ),
		sfxn_specialized( sfxn_specialized_arg ),
		sfxns_supported()
	{}

	DatabaseRelaxScript(
		std::string const & title_arg,
		bool sfxn_specialized_arg,
		std::string const & sfxns_arg
	) :
		title( title_arg ),
		sfxn_specialized( sfxn_specialized_arg ),
		sfxns_supported( utility::string_split( sfxns_arg, ':' ) )
	{}
};

/// @brief A class for storing energy maps, since they're not created by owning pointer
/// by default.

class EnergyMapContainer : public utility::pointer::ReferenceCount {
public:

	EnergyMapContainer() = delete;

	/// @brief Constructor -- READS FROM DISK!
	EnergyMapContainer( std::string const & filename );
	EnergyMapContainer( EnergyMapContainer const & src ) = default;
	~EnergyMapContainer() = default;

	/// @brief Access the energy map.
	inline core::scoring::EnergyMap const & get_energy_map() const { return energy_map_; }

private:

	core::scoring::EnergyMap energy_map_;
};

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

	utility::vector1< std::string > const &
	get_file_lines() const {
		return file_lines_;
	}

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
	/// @details Threadsafe and lazily loaded.  Requires the FastRelax or FastDesign mover's scorefunction to be provided.
	RelaxScriptFileContents const & get_relax_script(
		std::string const & filename,
		core::scoring::ScoreFunctionCOP const & mover_sfxn,
		bool const dualspace ) const;

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
	RelaxScriptFileContentsOP create_relax_script_instance( std::string const & script_file, std::string const & nearest_sfxn_name, bool const dualspace ) const;

	/// @brief Check whether (a) there is a period in the script file name, and, if so, (b) whether the base name exists in the
	/// database.  If (a) is true, throw an error message, where the text depends on whether (b) is true.
	/// @details PERFORMS DISK READ.
	static void handle_script_file_with_period( std::string const & script_file, std::string::size_type position_of_first_period );

	/// @brief Get the full script file path from the script name.
	/// @details PEFRORMS DISK READ.
	std::string script_file_path( std::string const & script_file, std::string const & nearest_sfxn_name, bool const dualspace ) const;

	/// @brief Get the database scorefunction that is closest to the current one, if the relax script is a
	/// database relax script that specifies a scorefunction.
	/// @details Performs onetime read from disk.
	std::string get_nearest_sfxn_if_in_database( std::string const & script_file, core::scoring::ScoreFunctionCOP mover_sfxn ) const;

	/// @brief Read the list of relax scripts from the database.
	/// @details Performs onetime read from disk.
	void initialize_relax_scripts_in_database() const;


	/// @brief Get the difference between two EnergyMaps.
	static core::Real sfxn_distance_squared( core::scoring::EnergyMap const & target_weights, core::scoring::EnergyMap const & candidate_weights );

	/// @brief Create an instance of an EnergyMapContainer.  Needed for threadsafe lazy loading.
	/// @details Peforms onetime read from disk.
	static EnergyMapContainerCOP create_energy_map_instance( std::string const & filename );

	/// @brief Given a scorefunction and a list of candidate scorefunctions, determine the closest one.
	/// @details Performs onetime read from disk.
	std::string determine_closest_scorefunction( core::scoring::ScoreFunction const & target, utility::vector1< std::string > const & names_of_candidates ) const;

private:  // Private data /////////////////////////////////////////////////////

	/// @brief A map of filename to file contents.
	mutable std::map < std::tuple< std::string, std::string, bool>, RelaxScriptFileContentsOP > filename_to_filecontents_map_;

	/// @brief cached data from database/sampling/relax_scripts/index.dat
	mutable utility::vector1< DatabaseRelaxScript > relax_scripts_in_database_;

	mutable std::map < std::string, EnergyMapContainerCOP > energy_maps_map_;

#ifdef MULTI_THREADED
	/// @brief Mutex for accessing the filename_to_filecontents_map_ object.
	mutable utility::thread::ReadWriteMutex relax_script_mutex_;

	/// @brief Mutex for accessing the list of relax scripts in the database.
	mutable utility::thread::ReadWriteMutex relax_scripts_in_database_mutex_;

	/// @brief Mutex for accessing the energy_maps_map_ object.
	mutable utility::thread::ReadWriteMutex energy_maps_mutex_;
#endif //MULTI_THREADED


};

} //protocols
} //relax

#endif //INCLUDED_protocols/relax_RelaxScriptManager_fwd_hh



