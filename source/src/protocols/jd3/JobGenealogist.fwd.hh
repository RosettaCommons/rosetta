// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobGenealogist.fwd.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_jd3_JobGenealogist_fwd_hh
#define INCLUDED_protocols_jd3_JobGenealogist_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace jd3 {

class JGResultNode;
using JGResultNodeOP = utility::pointer::shared_ptr< JGResultNode >;
using JGResultNodeCOP = utility::pointer::shared_ptr< JGResultNode const >;
using JGResultNodeAP = utility::pointer::weak_ptr< JGResultNode >;
using JGResultNodeCAP = utility::pointer::weak_ptr< JGResultNode const >;

class JGJobNode;
using JGJobNodeOP = utility::pointer::shared_ptr< JGJobNode >;
using JGJobNodeCOP = utility::pointer::shared_ptr< JGJobNode const >;
using JGJobNodeAP = utility::pointer::weak_ptr< JGJobNode >;
using JGJobNodeCAP = utility::pointer::weak_ptr< JGJobNode const >;

class JobGenealogist;
using JobGenealogistOP = utility::pointer::shared_ptr< JobGenealogist >;
using JobGenealogistCOP = utility::pointer::shared_ptr< JobGenealogist const >;

} // namespace jd3
} // namespace protocols

#endif
