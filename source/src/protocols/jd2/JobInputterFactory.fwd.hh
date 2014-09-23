// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobInputterFactory.fwd.hh
/// @brief  JobInputterFactory forward header
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_jd2_JobInputterFactory_fwd_hh
#define INCLUDED_protocols_jd2_JobInputterFactory_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class JobInputterFactory;

typedef utility::pointer::shared_ptr< JobInputterFactory > JobInputterFactoryOP;
typedef utility::pointer::shared_ptr< JobInputterFactory const > JobInputterFactoryCOP;

} //namespace jd2
} //namespace protocols

#endif //INCLUDED_protocols_jd2_JobInputterFactory_fwd_hh
