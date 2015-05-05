// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/WorkUnit.fwd.hh
/// @brief  Foward decls for base work unit (abstract) for wum2 and some commonly used derived work units
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_WorkUnit_fwd_hh
#define INCLUDED_protocols_wum2_WorkUnit_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace wum2 {

// ABSTRACT
class WorkUnit;
typedef boost::shared_ptr< WorkUnit > WorkUnitSP;

class WorkUnit_Wait;
typedef boost::shared_ptr< WorkUnit_Wait > WorkUnit_WaitSP;

class WorkUnit_ElScripts;
typedef boost::shared_ptr< WorkUnit_ElScripts > WorkUnit_ElScriptsSP;

} //namespace wum2
} //namespace protocols

#endif

