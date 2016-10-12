// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/QuotaConfig.hh
/// @brief Reads a config file for quota selector
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/quota/QuotaConfig.hh>

// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace quota {

using namespace core;

static THREAD_LOCAL basic::Tracer trQuotaConfig(
	"protocols.frag_picker.quota.QuotaConfig");

QuotaConfig::QuotaConfig(std::string file_name) {

	utility::io::izstream data(file_name.c_str());
	trQuotaConfig.Info << "reading quota config from " << file_name << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open quota config file: "
			+ file_name);
	}
	std::string line;
	core::Size pool_id;
	core::Real fraction;
	std::string pool_name;
	while ( data ) {
		getline(data, line);
		if ( line[0]=='#' ) continue;
		if ( line.length() > 7 ) {
			std::istringstream line_stream(line);
			line_stream >> pool_id >> pool_name >> fraction;
			pool_names_.push_back( pool_name );
			pool_weights_.push_back( fraction );
		}
	}
	data.close();
}

} // quota
} // frag_picker
} // protocols

