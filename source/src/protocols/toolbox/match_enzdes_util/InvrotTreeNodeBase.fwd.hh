// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.fwd.hh
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNodeBase_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNodeBase_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class InvrotTreeNodeBase;

typedef utility::pointer::shared_ptr< InvrotTreeNodeBase > InvrotTreeNodeBaseOP;
typedef utility::pointer::shared_ptr< InvrotTreeNodeBase const > InvrotTreeNodeBaseCOP;

typedef utility::pointer::weak_ptr< InvrotTreeNodeBase > InvrotTreeNodeBaseAP;
typedef utility::pointer::weak_ptr< InvrotTreeNodeBase const > InvrotTreeNodeBaseCAP;

class  AllowedSeqposForGeomCst;
typedef utility::pointer::shared_ptr< AllowedSeqposForGeomCst >AllowedSeqposForGeomCstOP;
typedef utility::pointer::shared_ptr< AllowedSeqposForGeomCst const > AllowedSeqposForGeomCstCOP;

class InvrotCollector;
typedef utility::pointer::shared_ptr< InvrotCollector >InvrotCollectorOP;
typedef utility::pointer::shared_ptr< InvrotCollector const > InvrotCollectorCOP;


}
}
}

#endif
