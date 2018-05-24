// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/segment_file_generator.fwd.hh
/// @brief a generator of smart sewing segment files
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_apps_pilot_frankdt_segment_file_generator_fwd_hh
#define INCLUDED_apps_pilot_frankdt_segment_file_generator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace apps {
namespace pilot {
namespace frankdt {

class segment_file_generator;

typedef utility::pointer::shared_ptr< segment_file_generator > segment_file_generatorOP;
typedef utility::pointer::shared_ptr< segment_file_generator const > segment_file_generatorCOP;



} //apps
} //pilot
} //frankdt


#endif //INCLUDED_apps_pilot_frankdt_segment_file_generator_fwd_hh





