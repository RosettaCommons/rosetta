// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file fwd.hh file for movers that mess around with additional ligand rigid body conformations
/// stored in the enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu, oct 09

#ifndef INCLUDED_protocols_enzdes_ModifyStoredLigandRBConfsMovers_fwd_hh
#define INCLUDED_protocols_enzdes_ModifyStoredLigandRBConfsMovers_fwd_hh

#include <utility/pointer/owning_ptr.hh>

//#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace enzdes {

class ModifyStoredRBConfs;
typedef utility::pointer::shared_ptr< ModifyStoredRBConfs > ModifyStoredRBConfsOP;
typedef utility::pointer::shared_ptr< ModifyStoredRBConfs const > ModifyStoredRBConfsCOP;

class GenerateStoredRBConfs;
typedef utility::pointer::shared_ptr< GenerateStoredRBConfs > GenerateStoredRBConfsOP;
typedef utility::pointer::shared_ptr< GenerateStoredRBConfs const > GenerateStoredRBConfsCOP;

class ApplyRandomStoredRBConf;
typedef utility::pointer::shared_ptr< ApplyRandomStoredRBConf > ApplyRandomStoredRBConfOP;
typedef utility::pointer::shared_ptr< ApplyRandomStoredRBConf const > ApplyRandomStoredRBConfCOP;

class MinimizeStoredRBConfs;
typedef utility::pointer::shared_ptr< MinimizeStoredRBConfs > MinimizeStoredRBConfsOP;
typedef utility::pointer::shared_ptr< MinimizeStoredRBConfs const > MinimizeStoredRBConfsCOP;

class DiversifyStoredRBConfs;
typedef utility::pointer::shared_ptr< DiversifyStoredRBConfs > DiversifyStoredRBConfsOP;
typedef utility::pointer::shared_ptr< DiversifyStoredRBConfs const > DiversifyStoredRBConfsCOP;

} // enzdes
} //protocols


#endif
