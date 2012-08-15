// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/Outputter.fwd.hh
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a file
/// @author Ken Jung

#ifndef INCLUDED_protocols_outputter_Outputter_fwd_hh
#define INCLUDED_protocols_outputter_Outputter_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace outputter {

class Outputter;
typedef boost::shared_ptr< Outputter > OutputterSP;

} // outputter
} // protocols

#endif //INCLUDED_protocols_outputter_Outputter_fwd_hh
