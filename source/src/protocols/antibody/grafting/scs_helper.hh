// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_helper.hh
/// @brief Singleton to store and load data for Structural Component Selector (SCS) filters
/// @author Jeliazko Jeliazkov

#ifndef INCLUDED_protocols_antibody_grafting_scs_helper_hh
#define INCLUDED_protocols_antibody_grafting_scs_helper_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_helper.fwd.hh>

#include <utility/SingletonBase.hh>
#include <core/types.hh>

#include <map>
#include <string>

namespace protocols {
namespace antibody {
namespace grafting {

/// @details This class is a singleton and manages data that is necessary for SCS filters and should be read only once
class SCS_Helper : public utility::SingletonBase< SCS_Helper >
{
public:
	friend class utility::SingletonBase< SCS_Helper >;
	
public: // Static constant data access
	/// @brief Return a map of pdb to map of CDR to bool (true if bfactor cutoff is NOT met)
	static std::map< std::string, std::map<std::string, bool> > get_ab_cdr_bfactors();
  
  /// @brief Return a map of pdb to map of pdb to a float (OCD distance between the two pdbs)
  static std::map< std::string, std::map<std::string, core::Real> > get_ab_OCDs();
  
  /// @brief Return a map of pdb to map of CDR to bool (true if region is outlier)
	static std::map< std::string, std::map<std::string, bool> > get_ab_region_outliers();
  
private: // Private methods
	// Empty constructor
	SCS_Helper();
	
	static SCS_Helper * create_singleton_instance();
	
	// Function to load bfactor data into map
	std::map< std::string, std::map<std::string, bool> > parse_bfactor_data();

  // Function to load OCD data into map
  std::map< std::string, std::map<std::string, core::Real> > parse_OCD_data();
  
  // Function to load outlier data into map
  std::map< std::string, std::map<std::string, bool> > parse_outlier_data();

private: // Private data
	std::map< std::string, std::map<std::string, bool> > bfactor_data_;
  std::map< std::string, std::map<std::string, core::Real> > OCD_data_;
  std::map< std::string, std::map<std::string, bool> > outlier_data_;
};

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_scs_functor_hh
