// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ABPSWrapper.fwd.hh
/// @brief  APBSWrapper forward delcaration
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_core_scoring_APBSWrapper_FWD_HH
#define INCLUDED_core_scoring_APBSWrapper_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
class APBSWrapper;
typedef utility::pointer::shared_ptr< APBSWrapper > APBSWrapperOP;
typedef utility::pointer::shared_ptr< APBSWrapper const > APBSWrapperCOP;
typedef utility::pointer::weak_ptr< APBSWrapper > APBSWrapperAP;
typedef utility::pointer::weak_ptr< APBSWrapper const > APBSWrapperCAP;

class APBSConfig;
typedef utility::pointer::shared_ptr< APBSConfig > APBSConfigOP;
typedef utility::pointer::shared_ptr< APBSConfig const > APBSConfigCOP;
typedef utility::pointer::weak_ptr< APBSConfig > APBSConfigAP;
typedef utility::pointer::weak_ptr< APBSConfig const > APBSConfigCAP;

class APBSResult;
typedef utility::pointer::shared_ptr< APBSResult > APBSResultOP;
typedef utility::pointer::shared_ptr< APBSResult const > APBSResultCOP;
typedef utility::pointer::weak_ptr< APBSResult > APBSResultAP;
typedef utility::pointer::weak_ptr< APBSResult const > APBSResultCAP;

class PQR;
typedef utility::pointer::shared_ptr< PQR > PQROP;
typedef utility::pointer::shared_ptr< PQR const > PQRCOP;
typedef utility::pointer::weak_ptr< PQR > PQRAP;
typedef utility::pointer::weak_ptr< PQR const > PQRCAP;
}
}
#endif
