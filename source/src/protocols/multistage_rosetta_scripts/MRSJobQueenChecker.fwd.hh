// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/MRSJobQueenChecker.fwd.hh
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_protocols_multistage_rosetta_scripts_MRSJobQueenChecker_FWD_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_MRSJobQueenChecker_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace multistage_rosetta_scripts {

class MRSJobQueenChecker;
typedef utility::pointer::shared_ptr< MRSJobQueenChecker > MRSJobQueenCheckerOP;
typedef utility::pointer::shared_ptr< MRSJobQueenChecker const > MRSJobQueenCheckerCOP;

} //multistage_rosetta_scripts
} //protocols

#endif
