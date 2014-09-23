// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/DetectProteinLigandInterface.fwd.hh
///
/// @brief
/// @author Sinisa Bjelic

#ifndef INCLUDED_protocols_enzdes_EnzdesTaskOperations_fwd_hh
#define INCLUDED_protocols_enzdes_EnzdesTaskOperations_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace enzdes {

class SetCatalyticResPackBehavior;
typedef utility::pointer::shared_ptr< SetCatalyticResPackBehavior > SetCatalyticResPackBehaviorOP;
typedef utility::pointer::shared_ptr< SetCatalyticResPackBehavior const > SetCatalyticResPackBehaviorCOP;

class DetectProteinLigandInterface;
typedef utility::pointer::shared_ptr< DetectProteinLigandInterface > DetectProteinLigandInterfaceOP;
typedef utility::pointer::shared_ptr< DetectProteinLigandInterface const > DetectProteinLigandInterfaceCOP;

typedef utility::pointer::weak_ptr< DetectProteinLigandInterface > DetectProteinLigandInterfaceAP;
typedef utility::pointer::weak_ptr< DetectProteinLigandInterface const > DetectProteinLigandInterfaceCAP;

class ProteinLigandInterfaceUpweighter;
typedef utility::pointer::shared_ptr< ProteinLigandInterfaceUpweighter > ProteinLigandInterfaceUpweighterOP;
typedef utility::pointer::shared_ptr< ProteinLigandInterfaceUpweighter const > ProteinLigandInterfaceUpweighterCOP;

class AddRigidBodyLigandConfs;
typedef utility::pointer::shared_ptr< AddRigidBodyLigandConfs > AddRigidBodyLigandConfsOP;
typedef utility::pointer::shared_ptr< AddRigidBodyLigandConfs const > AddRigidBodyLigandConfsCOP;

class AddLigandMotifRotamers;
typedef utility::pointer::shared_ptr< AddLigandMotifRotamers > AddLigandMotifRotamersOP;
typedef utility::pointer::shared_ptr< AddLigandMotifRotamers const > AddLigandMotifRotamersCOP;

}  //namespace enzdes
}  //namespace protocols

#endif
