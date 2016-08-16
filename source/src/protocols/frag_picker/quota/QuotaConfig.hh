// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/QuotaSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_QuotaConfig_hh
#define INCLUDED_protocols_frag_picker_quota_QuotaConfig_hh

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace frag_picker {
namespace quota {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

/// @brief read a config file for quota selector
class QuotaConfig {
public:

	/// @brief  Constructor reads a config file
	QuotaConfig(std::string config_file_name);

	/// @brief  Constructor used by derived classes
	QuotaConfig() {}

	/// @brief how many pools have been defined in a config file
	inline Size count_pools() { return pool_names_.size(); }

	/// @brief returns a fraction for a given pool
	inline Real get_fraction(Size pool_id) { return pool_weights_[pool_id]; }

	/// @brief returns a fraction for a given pool
	inline void set_fraction(Size pool_id,Real fraction) { pool_weights_[pool_id] = fraction; }

	/// @brief returns a fraction for a given pool
	/// @details if the given string is not a valid name of a quota pool,
	/// the method returns 0
	inline Real get_fraction(std::string pool_name) {

		for ( Size i=1; i<=pool_names_.size(); ++i ) {
			if ( pool_names_[i].compare(pool_name) == 0 ) {
				return pool_weights_[i];
			}
		}

		return 0;
	}

	/// @brief returns true if a config file defined a given pool name
	bool is_valid_quota_pool_name(std::string & pool_name) {

		for ( Size i=1; i<=pool_names_.size(); ++i ) {
			if ( pool_names_[i].compare(pool_name) == 0 ) {
				return true;
			}
		}

		return false;
	}

	/// @brief return a string id (name) assigned to a given pool
	inline std::string & get_pool_name(Size pool_id) { return pool_names_[pool_id]; }

protected:
	utility::vector1<Real> pool_weights_;
	utility::vector1<std::string> pool_names_;
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_QuotaConfig_HH */
