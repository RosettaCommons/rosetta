// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/distributions.hh
/// @brief  enumerate some distributions and map them to strings
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_DistributionMap_hh
#define INCLUDED_protocols_ligand_docking_DistributionMap_hh

// C++ headers
#include <map>
#include <string>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace ligand_docking {

enum Distribution{
	Uniform,
	Gaussian
};

/// A singleton class that returns a map of strings to enum types
class DistributionMap{
public:
	Distribution operator[](std::string distribution);
	static DistributionMap* get_instance();

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:
	DistributionMap(); // private constructor

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static DistributionMap * create_singleton_instance();

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< DistributionMap * > instance_;
#else
	static DistributionMap * instance_;
#endif

	std::map< std::string, Distribution > distribution_map_;


};

Distribution get_distribution(std::string distribution_str);

} //namespace ligand_docking
} //namespace protocols

#endif
