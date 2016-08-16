// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/PackDaemon.hh
/// @brief  declaration for class PackDaemon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_pack_daemon_PackDaemon_fwd_hh
#define INCLUDED_protocols_pack_daemon_PackDaemon_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace pack_daemon {

// #define APL_MEASURE_MSD_LOAD_BALANCE

class PackDaemon;

typedef utility::pointer::shared_ptr< PackDaemon > PackDaemonOP;
typedef utility::pointer::shared_ptr< PackDaemon const > PackDaemonCOP;

class DaemonSet;

typedef utility::pointer::shared_ptr< DaemonSet > DaemonSetOP;
typedef utility::pointer::shared_ptr< DaemonSet const > DaemonSetCOP;

class NPDPropCalculator;
typedef utility::pointer::shared_ptr< NPDPropCalculator > NPDPropCalculatorOP;
typedef utility::pointer::shared_ptr< NPDPropCalculator const > NPDPropCalculatorCOP;

class NPDPropCalculatorCreator;
typedef utility::pointer::shared_ptr< NPDPropCalculatorCreator > NPDPropCalculatorCreatorOP;
typedef utility::pointer::shared_ptr< NPDPropCalculatorCreator const > NPDPropCalculatorCreatorCOP;

class QuickRepacker;

typedef utility::pointer::shared_ptr< QuickRepacker > QuickRepackerOP;
typedef utility::pointer::shared_ptr< QuickRepacker const > QuickRepackerCOP;

class BasicSimAnnealerRepacker;

typedef utility::pointer::shared_ptr< BasicSimAnnealerRepacker > BasicSimAnnealerRepackerOP;
typedef utility::pointer::shared_ptr< BasicSimAnnealerRepacker const > BasicSimAnnealerRepackerCOP;

class DenseIGRepacker;

typedef utility::pointer::shared_ptr< DenseIGRepacker > DenseIGRepackerOP;
typedef utility::pointer::shared_ptr< DenseIGRepacker const > DenseIGRepackerCOP;

class DoubleDenseIGRepacker;

typedef utility::pointer::shared_ptr< DoubleDenseIGRepacker > DoubleDenseIGRepackerOP;
typedef utility::pointer::shared_ptr< DoubleDenseIGRepacker const > DoubleDenseIGRepackerCOP;

}
}

#endif
