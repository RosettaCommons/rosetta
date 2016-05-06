// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_helper.cc
/// @brief Singleton to store and load data for Structural Component Selector (SCS) filters
/// @author Jeliazko Jeliazkov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_helper.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <utility/vector0.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <core/types.hh>
#include <basic/database/open.hh>

// C++ headers
#include <map>
#include <string>
#include <fstream>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif

#include <basic/Tracer.hh>
static THREAD_LOCAL basic::Tracer TR("protocols.antibody.SCS_Helper");

// Singleton setup
namespace utility {
#if defined MULTI_THREADED && defined CXX11
	template <> std::mutex utility::SingletonBase< protocols::antibody::grafting::SCS_Helper > ::singleton_mutex_;
	template <> std::atomic< protocols::antibody::grafting::SCS_Helper * > utility::SingletonBase< protocols::antibody::grafting::SCS_Helper >::instance_( 0 );
#else
	template <> protocols::antibody::grafting::SCS_Helper * utility::SingletonBase< protocols::antibody::grafting::SCS_Helper >::instance_( 0 );
#endif
}

namespace protocols {
namespace antibody {
namespace grafting {

// Public methods ////////////////////////////////////////////////////////////
// Return bfactor map
std::map< std::string, std::map<std::string, bool> >
SCS_Helper::get_ab_cdr_bfactors()
{
	return get_instance()->bfactor_data_;
}

// Return OCD map
std::map< std::string, std::map<std::string, core::Real> >
SCS_Helper::get_ab_OCDs()
{
  return get_instance()->OCD_data_;
}

// Return outlier map
std::map< std::string, std::map<std::string, bool> >
SCS_Helper::get_ab_region_outliers()
{
  return get_instance()->outlier_data_;
}

// Private methods ////////////////////////////////////////////////////////////
// Constructor
SCS_Helper::SCS_Helper()
{
	// load data
	bfactor_data_ = parse_bfactor_data();
  OCD_data_ = parse_OCD_data();
  outlier_data_ = parse_outlier_data();
}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
SCS_Helper *
SCS_Helper::create_singleton_instance()
{
	return new SCS_Helper;
}


/// @details Parse bfactor data (true implies b-factor criterion is not met)
///          Expected input:
//						# ab	 l1    l2    l3    h1    h2    h3    // NEEDS TO BE ADDED TO FILE
//						12e8 False False False False False False
//						15c8 False False False False False True
//						1a0q False False False False True True
//						1a14 False False False False True False
//						1a2y False False False False False False
//						1a3l False False False False False False
//						1a3r False False False False False False
/// @return std::map< std::string pdb, std::map<std::string region, bool is_bfactor_50> >
std::map< std::string, std::map<std::string, bool> >
SCS_Helper::parse_bfactor_data()
{
	std::map< std::string, std::map<std::string, bool> > results;

	std::string dir = "protocol_data/antibody/";
	std::string filename = "list_bfactor50";
	std::string path = basic::database::find_database_path(dir, filename);

	// should be a vector zero of string to string maps
	auto lines( parse_plain_text_with_columns( path ) );

	// loop over pdb id and cdrs
	for(auto fields : lines ) {
		// Converting text based results to map
		std::map<std::string, bool> cdr_bfactors;

		std::string cdrs[] = { "l1", "l2", "l3", "h1", "h2", "h3" }; // is it ok to just use a list like this?
		for (auto cdr: cdrs) {
			assert( fields.at(cdr) == "True" | fields.at(cdr) == "False");
			cdr_bfactors[cdr] = fields.at(cdr) == "True";
		}

		results[fields.at("ab")] = cdr_bfactors;
	}

	return results;
}

/// @details Parse OCD data
/// input data is a symmetric matrix of orientational coordinate distances (OCD) between
/// all antibody structures in our template database as calculated by Nick Marze's packing_angle app
/// inital orientional coordinate calculations are stored in angles.sc and then comparisons.txt is generated
/// by the make_comparisons.py script
///          Expected input:
//            NOTE: LEGEND PREFIX '# ' NEEDS TO BE MANUALLY ADDED to be read by parse_plain_text_with_columns
//            pdb_code pdb12e8_chothia.pdb pdb15c8_chothia.pdb ...
//						pdb12e8_chothia.pdb 0.0 14.723
//            pdb15c8_chothia.pdb 14.723 0.0
/// @return std::map< std::string pdb1, std::map<std::string pdb2, core::Real OCD> >
std::map< std::string, std::map<std::string, core::Real> >
SCS_Helper::parse_OCD_data()
{
  std::map< std::string, std::map<std::string, core::Real> > results;

	std::string dir = "protocol_data/antibody/";
	std::string filename = "comparisons.txt";
	std::string path = basic::database::find_database_path(dir, filename);

  // should be a vector zero of string to string maps
  auto lines( parse_plain_text_with_columns( path ) );

  // loop over pdb id and cdrs
  for(auto fields : lines ) {
    // Converting text based results to map
    std::map<std::string, core::Real> OCDs;

    for (auto pairs : fields) {
      std::string key = pairs.first; // get column name
			if (key.compare("pdb_code")!=0) { // skip first column
				OCDs[key.substr(3,4)] = core::Real(std::stof(fields[key])); // OCDs[12e8] = 0.0;
			}
    }

    results[fields.at("pdb_code").substr(3,4)] = OCDs;
  }

  return results;
}

/// @details Parse outlier data
///          Expected input:
//            # pdb cdr outlier // NEEDS TO BE ADDED TO FILE
//            pdb12e8_chothia.pdb FRL false
//            pdb12e8_chothia.pdb FRH false
//            pdb12e8_chothia.pdb L1 false
//            pdb12e8_chothia.pdb L2 false
//            pdb12e8_chothia.pdb L3 false
//            pdb12e8_chothia.pdb H1 false
//            pdb12e8_chothia.pdb H2 false
//            pdb15c8_chothia.pdb FRL false
//            pdb15c8_chothia.pdb FRH false
//            pdb15c8_chothia.pdb L1 false
/// @return std::map< std::string pdb, std::map<std::string region, bool is_outlier> >
std::map< std::string, std::map<std::string, bool> >
SCS_Helper::parse_outlier_data()
{
	// e.g. 12e8_FRL
  std::map< std::string, std::map<std::string, bool> > results;
	string temp_cdr;

	// can we find the file?
	std::string dir = "protocol_data/antibody/";
	std::string filename = "outlier_list";
	std::string path = basic::database::find_database_path(dir, filename);

  // should be a vector zero of string to string maps
	// utility::vector0< std::map<string, string> >
  auto lines( parse_plain_text_with_columns( path ) );

	// set last stored pdb id, used for storing data later
	std::string last_pdb = lines[0]["pdb"];
	std::map<std::string, bool> outliers;

  // loop over lines (pdb id and cdrs)
  for(auto fields : lines ) {
		// store all regions for a single pdb in a single map
		// then store that single map to results map
		std::string current_pdb = fields["pdb"];

		if( current_pdb.compare(last_pdb)!=0 ) {
			// we've stored all outliers for all regions for a single pdb
			// map PDB - CDR - outlier
			results[last_pdb.substr(3,4)] = outliers;
			outliers.clear(); // re-initialize map for next pdb
			last_pdb = current_pdb;
		}

		// check for true/false value in outlier field
		assert( fields["outlier"] == "true" | fields["outlier"] == "false");

		// compare outlier field, store true if true, else false
		outliers[fields["cdr"]] = fields["outlier"] == "true";

  }

  return results;
}


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
