// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/util.hh
/// @brief  Utility functions for the distance-dependent dielectric potential.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit Headers
#include <core/scoring/elec/FA_ElecEnergy.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "core.scoring.elec.util" );

namespace core {
namespace scoring {
namespace elec {

/// @brief Read the CP tables from the database and return an owning pointer to the
/// new object created in memory.
/// @details Called by the ScoringManager to allow these data to be read in once and
/// only once.
CPRepMapTypeOP
read_cp_tables_from_db(
	std::string const & filename
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	CPRepMapTypeOP cp_rep_map_byname( new CPRepMapType );

	// search in the local directory first
	utility::io::izstream iunit;
	iunit.open( filename );

	if ( !iunit.good() ) {
		iunit.close();
		if ( !basic::database::open( iunit, filename ) ) {
			std::stringstream err_msg;
			err_msg << "Unable to open fa_elec countpair groups '" << filename << "'.";
			utility_exit_with_message(err_msg.str());
		}
	}

	std::string line;
	std::string name3, atom1, atom2;
	while ( iunit ) {
		getline( iunit, line );
		if ( line[0] == '#' ) continue;

		std::istringstream linestream(line);

		linestream >> name3 >> atom1 >> atom2;  // format is NAME REPRESENTATIVE_ATOM TARGET_ATOM

		if ( option[ score::elec_representative_cp_flip ]() ) {
			(*cp_rep_map_byname)[name3].insert(std::make_pair( atom1, atom2 )) ;
		} else {
			(*cp_rep_map_byname)[name3].insert(std::make_pair( atom2, atom1 )) ;
		}
	}

	TR << "Read " << cp_rep_map_byname->size() << " countpair representative atoms" << std::endl;
	return cp_rep_map_byname;
}

} // namespace elec
} // namespace scoring
} // namespace core
