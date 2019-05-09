// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/conformation/membrane/MembraneInfo.fwd.hh
/// @brief    Data describing the relationship between protein(s) and a membrane environment
///
/// @details  MembraneInfo is a container object that describes membrane-protein relationships
///             1. Coordinates of the membrane
///             2. A pointer to MEM which describes relative orientation (Residue)
///             3. Topology of the transmembrane spans (SpanningTopology)
///             4. Physical and chemical properties of the implicit lipid membrane (ImplicitLipidInfo)
///
/// @note     This object is a member of Conformation and should only be accessed using
///           pose.conformation().membrane_info(). Do not access MEM outside of the framework!
///
/// @author   Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {

/// @brief MembraneInfo describes the membrane bilayer and its relationship with the protein
class MembraneInfo;
typedef utility::pointer::shared_ptr< MembraneInfo > MembraneInfoOP;
typedef utility::pointer::shared_ptr< MembraneInfo const > MembraneInfoCOP;

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

