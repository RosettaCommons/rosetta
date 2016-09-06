// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/antibody/grafting/RegExManager.hh
/// @brief   Method declarations for RegExManager.
/// @author  Brian D. Weitzner <brian.weitzner@gmail.com>


#ifndef INCLUDED_protocols_antibody_grafting_regex_manager_HH
#define INCLUDED_protocols_antibody_grafting_regex_manager_HH

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

// Unit header
#include <protocols/antibody/grafting/regex_manager.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// C++ header
#include <string>


namespace protocols {
namespace antibody {
namespace grafting {

/// @details  This class is a singleton and manages CarbohydratesInfo data that should only be read from the database
/// one time and shared among all instances of CarbohydrateInfo.
class RegExManager : public utility::SingletonBase< RegExManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< RegExManager >;

public:  // Static constant data access ///////////////////////////////////////
	std::string H1_pattern() const;
	std::string H3_pattern() const;
	
	std::string L1_pattern() const;
	std::string L3_pattern() const;

private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	RegExManager();

	// Singleton-creation function for use with utility::thread::threadsafe_singleton
	static RegExManager * create_singleton_instance();
	
	void load_regex_from_db();

private:  // Private data /////////////////////////////////////////////////////
	std::string H1_pattern_;
	std::string H3_pattern_;
	
	std::string L1_pattern_;
	std::string L3_pattern_;
	
};

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif  // INCLUDED_protocols_antibody_grafting_regex_manager_HH
