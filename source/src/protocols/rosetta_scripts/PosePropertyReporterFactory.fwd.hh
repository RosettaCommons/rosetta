// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rosetta_scripts/PosePropertyReporterFactory.fwd.hh
/// @brief	Factory for PosePropertyReporters forward definitions
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_moves_PosePropertyReporterFactory_fwd_hh
#define INCLUDED_protocols_moves_PosePropertyReporterFactory_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rosetta_scripts {

class PosePropertyReporterFactory;

typedef utility::pointer::owning_ptr< PosePropertyReporterFactory > PosePropertyReporterFactoryOP;
typedef utility::pointer::owning_ptr< PosePropertyReporterFactory const > PosePropertyReporterFactoryCOP;

} //namespace rosetta_scripts
} //namespace protocols

#endif
