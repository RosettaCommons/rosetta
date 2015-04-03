// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/file_util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_file_util_HH
#define INCLUDED_protocols_stepwise_modeler_file_util_HH

#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace protocols {
namespace stepwise {
namespace modeler {

	std::string
	get_file_name( std::string const & silent_file, std::string const & tag );

	void
	remove_silent_file_if_it_exists( std::string const & silent_file);

} //modeler
} //stepwise
} //protocols

#endif
